# ============================================================
# 4-NODE GENE CIRCUIT SIMULATION PIPELINE
# Description: Executes Sobol parameter sampling and runs a 
#              two-pulse ODE experiment.
# Matrix Convention: A[i, j] -> interaction from node j to i
# Dependencies: DifferentialEquations, QuasiMonteCarlo, 
#               DataFrames, Statistics, CSV, Dates
# ============================================================

using DifferentialEquations
using QuasiMonteCarlo
using DataFrames
using Statistics
using CSV
using Dates

# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

"""
    log_sample(u::Real, min_val::Real, max_val::Real) -> Float64

Transforms a uniformly distributed variable `u` ∈ [0, 1] into a logarithmically 
scaled value between `min_val` and `max_val`.
"""
function log_sample(u::Real, min_val::Real, max_val::Real)::Float64
    return 10^(log10(min_val) + u * (log10(max_val) - log10(min_val)))
end

"""
    get_stimulus(t, T1, T3, duration, strength) -> Real

Evaluates the square pulse stimulus at time `t`. Generates a pulse of `strength` 
at start times `T1` and `T3`, each lasting for `duration`.
"""
function get_stimulus(t, T1, T3, duration, strength)
    pulse1 = (T1 <= t < T1 + duration)
    pulse2 = (T3 <= t < T3 + duration)
    return (pulse1 || pulse2) ? strength : 0.0
end

"""
    model!(du, u, p, t)

In-place ODE function for the 4-node gene circuit.
Indexing Convention:
- `A[i, j]` represents the directed edge FROM node `j` (source) TO node `i` (target).
- `du[i]` (target) is modified by iterating over `j` (source) using `A[i, j]`.
"""
function model!(du, u, p, t)
    A     = p.A
    alpha = p.alpha
    kappa = p.kappa
    act   = p.act
    inh   = p.inh
    stim  = p.stim

    S = get_stimulus(t, stim.T1, stim.T3, stim.duration, stim.strength)

    for i in 1:4 # i is the TARGET node
        ui = u[i]
        gate = 1.0

        if i == 1
            gate = S
        elseif i == 4
            gate = u[1]
        end

        activation = 0.0

        # Basal production (only applied to Y1 and Y2)
        if i != 1 && i != 4
            activation += act[i] * (1.0 - ui)
        end

        deactivation = 4 * inh[i] * ui

        # Stimulus direct activation (only applied to X)
        if i == 1
            activation += stim.alpha * S * (1.0 - ui) / (stim.kappa + (1.0 - ui))
        end

        # Network interactions
        for j in 1:4 # j is the SOURCE node
            # Skip X (1) -> Z (4) direct interaction entirely
            if j == 1 && i == 4
                continue
            end

            # A[i, j] is the edge from source j to target i
            if A[i, j] == 1.0
                activation += gate * alpha[i, j] * u[j] * (1.0 - ui) / (kappa[i, j] + (1.0 - ui))
            elseif A[i, j] == -1.0
                deactivation += alpha[i, j] * u[j] * ui / (kappa[i, j] + ui)
            end
        end

        du[i] = activation - deactivation
    end
    return nothing
end


"""
    find_T3(sol::ODESolution, eta_t3::Float64, T1::Float64) -> Float64

Determines the start of the second pulse (T3) based on the decay of both X and Z. 
Restricts the search to find the peaks of X and Z strictly at or after T1. T3 is 
identified as the first time point after both peaks where both X and Z have 
decayed below their respective thresholds.
"""
function find_T3(sol::ODESolution, eta_t3::Float64, T1::Float64)::Float64
    # Assuming X is at index 1 and Z is at index 4 based on standard motif layouts.
    # Adjust these indices if your state variables are ordered differently.
    X = sol[1, :] 
    Z = sol[4, :]
    t = sol.t

    # Find the first index where t >= T1 (O(log N) fast search)
    start_idx = searchsortedfirst(t, T1)
    if start_idx > length(t)
        return t[end]
    end

    # Isolate the window of interest for X and Z
    X_window = @view X[start_idx:end]
    Z_window = @view Z[start_idx:end]
    
    # Find peak indices relative to the original arrays
    peak_idx_X = start_idx + argmax(X_window) - 1
    peak_idx_Z = start_idx + argmax(Z_window) - 1
    
    # Calculate individual thresholds based on their respective peaks
    threshold_X = eta_t3 * X[peak_idx_X]
    threshold_Z = eta_t3 * Z[peak_idx_Z]

    # The search for decay must begin only after BOTH variables have reached their peaks.
    # This prevents false positives if one variable dips while the other is still rising.
    decay_search_start = max(peak_idx_X, peak_idx_Z)

    # Iterate to find when BOTH variables drop below their thresholds
    for i in decay_search_start:(length(t)-1)
        if X[i+1] <= threshold_X && Z[i+1] <= threshold_Z
            return t[i+1]
        end
    end

    # Return the final time point if the condition is never met
    return t[end]
end



"""
    compute_auc(t::Vector{Float64}, Z::Vector{Float64}, t_start::Float64, t_end::Float64) -> Float64

Computes the Area Under the Curve (AUC) for Z using the trapezoidal rule.
"""
function compute_auc(t::Vector{Float64}, Z::Vector{Float64}, t_start::Float64, t_end::Float64)::Float64
    t_stop = min(t_end, t[end])
    if t_stop <= t_start
        return 0.0
    end

    auc = 0.0
    for i in 1:(length(t)-1)
        if t[i] >= t_start && t[i+1] <= t_stop
            auc += (Z[i] + Z[i+1]) / 2.0 * (t[i+1] - t[i])
        end
    end
    return auc
end

"""
    compute_peak(t::Vector{Float64}, Z::Vector{Float64}, t_start::Float64, t_end::Float64) -> Float64

Finds the maximum amplitude of Z within the specified time window.
"""
function compute_peak(t::Vector{Float64}, Z::Vector{Float64}, t_start::Float64, t_end::Float64)::Float64
    t_stop = min(t_end, t[end])
    peak = 0.0
    for i in 1:length(t)
        if t_start <= t[i] <= t_stop
            peak = max(peak, Z[i])
        end
    end
    return peak
end


# ============================================================
# 2. MOTIF SIMULATION
# ============================================================

"""
    simulate_motif(row::DataFrameRow) -> NamedTuple

Runs the parameter sampling and dual-pulse simulation for a single motif row.
Constructs the adjacency matrix such that A[tgt, src] is the edge from src to tgt.
CSV format Cij means edge FROM i TO j, so it maps to A[j, i].
"""
function simulate_motif(row::DataFrameRow)
    n_samples = Int(row.n_sampling)
    T1 = Float64(row.T1)

    # A[tgt, src] mapping
    # Since Cij is Source i -> Target j, Cij must map to A[j, i]
    A = zeros(Float64, 4, 4)
    A[1, :] .= [row.C11, row.C21, row.C31, row.C41] # Target 1 edges from Sources 1, 2, 3, 4
    A[2, :] .= [row.C12, row.C22, row.C32, row.C42] # Target 2 edges from Sources 1, 2, 3, 4
    A[3, :] .= [row.C13, row.C23, row.C33, row.C43] # Target 3 edges from Sources 1, 2, 3, 4
    A[4, :] .= [row.C14, row.C24, row.C34, row.C44] # Target 4 edges from Sources 1, 2, 3, 4

    basal = Float64(row.basal_e)
    act = fill(basal, 4)
    inh = fill(basal, 4)

    # Determine valid target <- source edges
    valid_edges = Tuple{Int, Int}[]
    for tgt in 1:4, src in 1:4
        if A[tgt, src] != 0.0 && !(src == 1 && tgt == 4)
            push!(valid_edges, (tgt, src))
        end
    end
    
    n_params = 2 * length(valid_edges) + 2
    sobol_samples = QuasiMonteCarlo.sample(n_samples, n_params, SobolSample())

    # Pre-allocate arrays
    SI_values  = zeros(Float64, n_samples)
    PRR_values = zeros(Float64, n_samples)
    T3_values  = zeros(Float64, n_samples)
    count_both_gt1 = 0
    count_both_lt1 = 0
    count_SI_gt1_PRR_lt1 = 0
    count_SI_lt1_PRR_gt1 = 0
    count_SI_gt1_PRR_et1 = 0
    count_SI_lt1_PRR_et1 = 0

    u0 = Float64[row.x_o, row.y1_o, row.y2_o, row.z_o]

    for k in 1:n_samples
        alpha = zeros(Float64, 4, 4)
        kappa = zeros(Float64, 4, 4)

        @views sobol_vec = sobol_samples[:, k]
        param_idx = 1
        
        for (tgt, src) in valid_edges
            alpha[tgt, src] = log_sample(sobol_vec[param_idx], row.min_alpha, row.max_alpha)
            kappa[tgt, src] = log_sample(sobol_vec[param_idx + 1], row.min_kappa, row.max_kappa)
            param_idx += 2
        end

        stim_alpha = log_sample(sobol_vec[param_idx], row.min_alpha, row.max_alpha)
        stim_kappa = log_sample(sobol_vec[param_idx + 1], row.min_kappa, row.max_kappa)

        # Pulse 1 Specification
        stim_first = (
            alpha = stim_alpha, kappa = stim_kappa,
            strength = Float64(row.max_s), duration = Float64(row.t_duration),
            T1 = T1, T3 = Inf
        )

        p_first = (A=A, alpha=alpha, kappa=kappa, act=act, inh=inh, stim=stim_first)
        
        stops_1 = [T1, T1 + Float64(row.t_duration)]
        sol_first = solve(ODEProblem(model!, u0, (0.0, Float64(row.T)), p_first), 
                          Rodas5P(); abstol=1e-8, reltol=1e-6, tstops=stops_1)
        
        T3 = find_T3(sol_first, Float64(row.eta_t3), T1)
            
        # Pulse 2 Specification
        stim_second = (
            alpha = stim_alpha, kappa = stim_kappa,
            strength = Float64(row.max_s), duration = Float64(row.t_duration),
            T1 = T1, T3 = T3
        )

        p_second = (A=A, alpha=alpha, kappa=kappa, act=act, inh=inh, stim=stim_second)
        
        stops_2 = [T1, T1 + Float64(row.t_duration), T3, T3 + Float64(row.t_duration)]
        sol = solve(ODEProblem(model!, u0, (0.0, Float64(row.T)), p_second), 
                    Rodas5P(); abstol=1e-8, reltol=1e-6, tstops=stops_2)

        # Metrics
        t = sol.t
        Z = sol[4, :]

        delta_T = max(T3 - T1, 0.0)

        auc_pre  = compute_auc(t, Z, T1, T3)
        auc_post = compute_auc(t, Z, T3, T3 + delta_T)
        peak_pre  = compute_peak(t, Z, T1, T3)
        peak_post = compute_peak(t, Z, T3, T3 + delta_T)

        SI  = auc_post / (auc_pre)
        PRR = peak_post / (peak_pre)

        SI_values[k] = SI
        PRR_values[k] = PRR
        T3_values[k] = T3

        if SI <= 1.0 && PRR < 1.0; count_both_lt1 += 1;
        elseif SI > 1.0 && PRR < 1.0; count_SI_gt1_PRR_lt1 += 1;
        elseif SI <= 1.0 && PRR > 1.0; count_SI_lt1_PRR_gt1 += 1;
        elseif SI > 1.0 && PRR > 1.0; count_both_gt1 += 1; 
        elseif SI > 1.0 && PRR == 1.0; count_SI_gt1_PRR_et1 += 1;
        elseif SI < 1.0 && PRR == 1.0; count_SI_lt1_PRR_et1 += 1;    
        end
        
    end

    return (
        Motif = Int(row[Symbol("Motif Number")]),
        
        # SI Statistics
        
        SI_min = minimum(SI_values),
        SI_q25 = quantile(SI_values, 0.25),
        SI_median = median(SI_values),
        SI_q75 = quantile(SI_values, 0.75),
        SI_max = maximum(SI_values),
        SI_mean = mean(SI_values),
        SI_std = std(SI_values),
        
        
        # PRR Statistics
        
        PRR_min = minimum(PRR_values),
        PRR_q25 = quantile(PRR_values, 0.25),
        PRR_median = median(PRR_values),
        PRR_q75 = quantile(PRR_values, 0.75),
        PRR_max = maximum(PRR_values),
        PRR_mean = mean(PRR_values),
        PRR_std = std(PRR_values),

        
        # T3 Statistics
        
        T3_min = minimum(T3_values),
        T3_q25 = quantile(T3_values, 0.25),
        T3_median = median(T3_values),
        T3_q75 = quantile(T3_values, 0.75),
        T3_max = maximum(T3_values),
        T3_mean = mean(T3_values),
        T3_std = std(T3_values),
        
        
        successful_runs = n_samples,
        count_both_gt1 = count_both_gt1,
        count_both_lt1 = count_both_lt1,
        count_SI_gt1_PRR_lt1 = count_SI_gt1_PRR_lt1,
        count_SI_lt1_PRR_gt1 = count_SI_lt1_PRR_gt1,
        count_SI_gt1_PRR_et1 = count_SI_gt1_PRR_et1,
        count_SI_lt1_PRR_et1 = count_SI_lt1_PRR_et1
    )
end

# ============================================================
# 3. EXECUTION PIPELINE
# ============================================================

const INPUT_FILE = "input_params.csv"

# 1. Generate a unique timestamp (e.g., 2026-04-06_20-13-12)
timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")

# 2. Append the timestamp to ensure unique output files
const OUTPUT_FILE = "simulation_results_$timestamp.csv"
const LOG_FILE = "simulation_log_$timestamp.txt"

println("Loading dataset from: ", INPUT_FILE)
data = CSV.read(INPUT_FILE, DataFrame)
rename!(data, Symbol.(strip.(string.(names(data)))))

results_array = NamedTuple[]

for row in eachrow(data)
    # Handle slight variations in the CSV column name
    motif_col = hasproperty(row, Symbol("Motif Number")) ? Symbol("Motif Number") : Symbol("Motif Num")
    println("Processing Motif: ", row[motif_col])
    
    motif_result = simulate_motif(row)
    push!(results_array, motif_result)
end

results_df = DataFrame(results_array)

# 3. Stitch the input parameters and results together horizontally
final_df = hcat(data, results_df[:, 2:end])

# 4. ADD THE TIMESTAMP AS A COLUMN TO THE DATA
# This guarantees you know exactly when the simulation ran, even if the file is renamed!
final_df[!, :Run_Timestamp] .= timestamp

println("Writing results to: ", OUTPUT_FILE)
CSV.write(OUTPUT_FILE, final_df)

open(LOG_FILE, "w") do io
    println(io, "Simulation Timestamp: ", Dates.now())
    println(io, "Input File Used: ", INPUT_FILE)
    println(io, "Total Motifs Processed: ", nrow(data))
    
    for res in results_array
        println(io, "\nMotif: ", res.Motif)
        println(io, "  Successful runs: ", res.successful_runs)
        println(io, "  T3 -> mean: ", res.T3_mean, ", median: ", res.T3_median)
        println(io, "  SI mean: ", res.SI_mean)
        println(io, "  PRR mean: ", res.PRR_mean)
    end
end

println("Simulation complete. Results saved to $OUTPUT_FILE")