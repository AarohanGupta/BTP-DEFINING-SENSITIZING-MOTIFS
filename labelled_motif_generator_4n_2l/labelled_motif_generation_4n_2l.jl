############################################################
# 4-Node Sensitizing Motif Enumeration (CORRECTED)
############################################################

using DelimitedFiles
using Combinatorics


############################################################
# MODULE 1: Matrix Utilities
############################################################

function vec_to_mat(v)
    reshape(collect(v), 4, 4)'
end

function mat_to_tuple(A)
    Tuple(vec(A'))
end


############################################################
# MODULE 2: Canonical Representation
############################################################

# Remove redundancy due to Y1 ↔ Y2 symmetry
function canonical_form(A)

    perms = [[1,2,3,4], [1,3,2,4]]  # swap Y1 ↔ Y2

    forms = NTuple{16,Int}[]

    for p in perms
        Ap = A[p, p]
        push!(forms, mat_to_tuple(Ap))
    end

    return minimum(forms)
end


############################################################
# MODULE 3: Structural Conditions
############################################################

# Condition 1: Direct positive X → Z
# X = 1, Z = 4
function direct_positive_XZ(A)
    return A[4,1] == 1
end

function feedback_reachable_from_x(A)

    # Direct from X
    if A[2,1] != 0 || A[3,1] != 0
        return true
    end

    # Indirect via Z (since X → Z exists)
    if A[2,4] != 0 || A[3,4] != 0
        return true
    end

    return false
end

# Condition 2: Positive feedback between Y1 and Y2
# Y1 = 2, Y2 = 3
# --Should the Feedback be exclusively +ve?--
#Will revise when finalized
function positive_Y_feedback(A)

    if A[3,2] == 0 || A[2,3] == 0 || A[2,3] == -1 || A[3,2] == -1 || A[4,3] == -1 || A[4,2] == -1
        return false
    end

    # SIGNED NETWORK RULE
    #The net impact of feedback cycle must be positive on Z

    return A[3,2] * A[4,3] > 0 || A[2,3] * A[4,2] > 0 

end





############################################################
# MODULE 4: Combined Sensitization Filter
############################################################

function is_sensitizing(A)

    return direct_positive_XZ(A) &&
           positive_Y_feedback(A) &&
           feedback_reachable_from_x(A)
           
end


############################################################
# MODULE 5: Motif Generation
############################################################

function generate_all_labeled_motifs(motif_vectors)

    perms = collect(permutations(1:4))
    all_motifs = Set{NTuple{16,Int}}()

    for v in eachrow(motif_vectors)

        A = vec_to_mat(v)

        for p in perms
            Ap = A[p, p]
            push!(all_motifs, mat_to_tuple(Ap))
        end
    end

    return all_motifs
end


############################################################
# MODULE 6: Redundancy Removal
############################################################

function remove_redundant_motifs(all_motifs)

    canonical_set = Set{NTuple{16,Int}}()

    for key in all_motifs
        A = vec_to_mat(key)
        push!(canonical_set, canonical_form(A))
    end

    return canonical_set
end


############################################################
# MODULE 7: Filtering Sensitizing Motifs
############################################################

function filter_sensitizing_motifs(motif_set)

    final_set = NTuple{16,Int}[]

    for key in motif_set
        A = vec_to_mat(key)

        if is_sensitizing(A)
            push!(final_set, key)
        end
    end

    return final_set
end


############################################################
# MODULE 8: Output Writer
############################################################

function write_output(filename, motifs)

    open(filename, "w") do io
        id = 1
        for k in sort(motifs)
            println(io, join((id, k...), ","))
            id += 1
        end
    end
end


############################################################
# MAIN PIPELINE
############################################################

function generate_motifs()

    println("Reading input...")

    data = readdlm("4N_2L_topo.txt", ',', Int)
    motif_vectors = data[:, 2:end]

    # Step 1: Generate labeled motifs
    all_motifs = generate_all_labeled_motifs(motif_vectors)
    println("Total labeled motifs = ", length(all_motifs))

    # Step 2: Remove redundancy
    canonical_set = remove_redundant_motifs(all_motifs)
    println("After removing redundancy = ", length(canonical_set))

    # Step 3: Apply sensitization filters
    final_set = filter_sensitizing_motifs(canonical_set)
    println("Final sensitizing motifs = ", length(final_set))

    # Step 4: Write output
    write_output("4n_2l_labelled_sensitizing_motifs.txt", final_set)

    println("Output written.")

end


############################################################
# RUN
############################################################

generate_motifs()