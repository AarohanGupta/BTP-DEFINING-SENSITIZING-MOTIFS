4-NODE GENE CIRCUIT ODE SIMULATION PIPELINE
===========================================

------------------------------------------------------------
OVERVIEW
------------------------------------------------------------
This repository provides a high-performance Julia pipeline 
for simulating a 4-node gene regulatory circuit using 
Ordinary Differential Equations (ODEs).

The pipeline:
- Uses Sobol sequence sampling (Quasi-Monte Carlo)
- Simulates dual-pulse stimuli experiments
- Evaluates motif behavior using SI and PRR metrics

The entire workflow is implemented in a single optimized script.


------------------------------------------------------------
KEY FEATURES
------------------------------------------------------------
- Deterministic Sobol sampling (fully reproducible)
- Stiff ODE solving using Rodas5P
- Dual-pulse simulation
- Automated feature extraction (SI, PRR, T3)
- Scalable to large motif datasets
- Thread-safe for parallel execution


------------------------------------------------------------
REPRODUCIBILITY
------------------------------------------------------------
Sobol sampling is completely deterministic:
- No randomness across runs
- Identical outputs for identical inputs
- Consistent motif comparison


------------------------------------------------------------
MODEL DESCRIPTION
------------------------------------------------------------

Nodes:
- X (1): Input
- Y1 (2): Intermediate
- Y2 (3): Intermediate
- Z (4): Output


------------------------------------------------------------
MATRIX CONVENTION
------------------------------------------------------------

Internal (ODE Model):
- A[i, j] represents interaction FROM node j TO node i
- Rows → Targets
- Columns → Sources

Input CSV:
- Cij represents interaction FROM node i TO node j
- Converted internally as:
  A[j, i] = Cij


------------------------------------------------------------
INTERACTION VALUES
------------------------------------------------------------
  1  : Activation
 -1  : Inhibition
  0  : No interaction


------------------------------------------------------------
PROJECT STRUCTURE
------------------------------------------------------------

final_single_simulation_pipeline.jl

Contains:
- Helper functions
- ODE model (model!)
- Sobol sampling logic
- Simulation loop
- Output generation

Run this file to execute the full pipeline.


------------------------------------------------------------
DEPENDENCIES
------------------------------------------------------------

Install required packages:

using Pkg
Pkg.add([
    "DifferentialEquations",
    "CSV",
    "DataFrames",
    "QuasiMonteCarlo",
    "StatsBase"
])

Standard libraries:
- Statistics
- Dates


------------------------------------------------------------
INPUT FORMAT
------------------------------------------------------------

File: input_params.csv

Place in the same directory as the script.

Required columns:

Motif Information:
- Motif Number

Adjacency Matrix:
- C11 to C44

Initial Conditions:
- x_o, y1_o, y2_o, z_o

Global Parameters:
- basal_e
- n_sampling
- T
- T1
- t_duration
- max_s
- min_alpha, max_alpha
- min_kappa, max_kappa
- eta_t3


------------------------------------------------------------
IMPORTANT ASSUMPTION
------------------------------------------------------------
Basal deactivation rate = 4 × basal activation rate
(for Y1 and Y2 nodes)


------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------

1. simulation_results_<timestamp>.csv

Contains:
- SI, PRR, T3 statistics:
  - Min, Max, Mean, Median, Quartiles, Std
- Counts:
  - SI > 1
  - PRR > 1
  - Both > 1
  - Mixed regimes


2. simulation_log_<timestamp>.txt

Contains:
- Execution timestamp
- Input file used
- Per-motif summaries:
  - SI mean
  - PRR mean
  - T3 statistics
  - Successful runs


------------------------------------------------------------
HOW TO RUN
------------------------------------------------------------

1. Open terminal
2. Navigate to project directory
3. Ensure input file is correctly named
4. Run:

julia pipeline.jl


------------------------------------------------------------
SIMULATION DETAILS
------------------------------------------------------------

Solver:
- Rodas5P (stiff Rosenbrock method)

Stimulus:
- Two square pulses:
  - First at T1
  - Second at dynamically computed T3

T3 Logic:
- Determined based on decay of X and Z
- Ensures both signals peak before decay threshold check


------------------------------------------------------------
METRICS
------------------------------------------------------------

SI (Sensitivity Index):
- Ratio of post-pulse to pre-pulse AUC

PRR (Peak Response Ratio):
- Ratio of post-pulse to pre-pulse peak

T3:
- Start time of second pulse


------------------------------------------------------------
PERFORMANCE OPTIMIZATION
------------------------------------------------------------
- Local Sobol sampling per motif
- Memory-efficient for large sampling sizes
- Pre-allocation for speed

------------------------------------------------------------
NOTES
------------------------------------------------------------
- Direct interaction X → Z is excluded
- Designed for:
  - Large-scale motif screening
  - Systems biology analysis
  - Synthetic circuit evaluation