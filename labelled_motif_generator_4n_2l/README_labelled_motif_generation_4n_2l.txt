README: 4-Node Sensitizing Motif Enumeration (Labelled Motifs)

============================================================

OVERVIEW
------------------------------------------------------------
This pipeline generates and filters 4-node network motifs based on
biologically inspired sensitization criteria.

Final Output:
- 3276 labelled sensitizing motifs

Pipeline:
1. Read input motifs
2. Generate labeled variants
3. Remove redundancy (symmetry)
4. Apply structural filters
5. Write output

============================================================

EXECUTION LOG
------------------------------------------------------------
Reading input...
Total labeled motifs = 528384
After removing redundancy = 265248
Final sensitizing motifs = 3276
Output written.

Interpretation:
- 528384 → All possible labeled permutations generated
- 265248 → Unique motifs after removing symmetry (Y1 ↔ Y2)
- 3276   → Final biologically valid sensitizing motifs

============================================================

CONCEPTUAL MODEL
------------------------------------------------------------
Nodes:
X (1)  : Input
Y1 (2) : Intermediate
Y2 (3) : Intermediate
Z (4)  : Output

Representation:
- 4×4 adjacency matrix
- 1  → activation
- -1 → inhibition
- 0  → no interaction

Matrix Convention:
A[i, j] represents interaction from node j → node i
(rows = targets, columns = sources)

============================================================

INPUT FILE
------------------------------------------------------------
File: 4N_2L_topo.txt

Format:
- CSV
- First column = ID (ignored)
- Next 16 values = flattened 4×4 matrix

============================================================

CODE STRUCTURE
------------------------------------------------------------

MODULE 1: Matrix Utilities
- vec_to_mat(v): vector → matrix
- mat_to_tuple(A): matrix → tuple

------------------------------------------------------------

MODULE 2: Canonical Representation
- Removes redundancy due to Y1 ↔ Y2 symmetry
- Uses permutations:
  [1,2,3,4]
  [1,3,2,4]

- Keeps minimum lexicographic form

------------------------------------------------------------

MODULE 3: Structural Conditions

1. Direct Positive X → Z
   A[4,1] == 1

2. Feedback reachable from X
   - Direct: X → Y1 or Y2
   - Indirect: X → Z → Y1 or Y2

3. Positive Feedback (Y1 ↔ Y2)
   - Both directions must exist
   - No inhibitory edges in key paths
   - Net effect on Z must be positive

------------------------------------------------------------

MODULE 4: Combined Filter
- is_sensitizing(A)
- Returns true if all conditions are satisfied

------------------------------------------------------------

MODULE 5: Motif Generation
- Generates all 4! = 24 permutations
- Stores unique labeled motifs

------------------------------------------------------------

MODULE 6: Redundancy Removal
- Converts motifs to canonical form
- Removes duplicates

------------------------------------------------------------

MODULE 7: Filtering
- Applies sensitization rules
- Keeps valid motifs

------------------------------------------------------------

MODULE 8: Output Writer
- Writes to:
  4n_2l_labelled_sensitizing_motifs.txt

- Format:
  ID, 16 matrix entries

============================================================

PIPELINE EXECUTION
------------------------------------------------------------
Function: generate_motifs()

Steps:
1. Read input
2. Generate labeled motifs
3. Remove redundancy
4. Apply filters
5. Write output

============================================================

OUTPUT
------------------------------------------------------------
File:
4n_2l_labelled_sensitizing_motifs.txt

Format:
ID, a11, a12, ..., a44

Result:
3276 labelled sensitizing motifs

============================================================

INTUITION
------------------------------------------------------------
- X activates Z
- Y1 and Y2 form a reinforcing loop
- Feedback amplifies signal to Z

Result:
Small input → Strong output (sensitization)

============================================================

HOW TO RUN
------------------------------------------------------------
Requirements:
- Julia
- Packages:
  using DelimitedFiles
  using Combinatorics

Run:
generate_motifs()

============================================================

KEY DESIGN CHOICES
------------------------------------------------------------
- Canonicalization removes duplicates
- Biological constraints ensure meaningful motifs
- Modular design for easy modification
- Tuple storage enables efficient hashing

============================================================

EXTENSIONS
------------------------------------------------------------
- Modify feedback rules
- Add new biological conditions
- Track SI / PRR metrics
- Extend to larger networks

============================================================

NOTES
------------------------------------------------------------
- Node mapping fixed:
  X=1, Y1=2, Y2=3, Z=4
- Input must be correctly formatted
- Feedback ensures positive net effect on Z

============================================================

SUMMARY
------------------------------------------------------------
The code:
- Explores all motif structures
- Removes symmetry redundancy
- Applies biological filters
- Outputs valid sensitizing motifs

Final Result:
3276 motifs
