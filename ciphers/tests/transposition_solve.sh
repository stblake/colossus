#!/bin/bash
# End-to-end test: solve a keyed columnar transposition with the general
# (AZDecrypt-style) permutation solver. The ciphertext was produced by writing
# the plaintext row-wise into 6 columns and reading the columns out in the
# order [3,0,5,1,4,2]; the solver does not know the period or column order.
#
# The solver hill-climbs the full permutation key, seeding restarts from
# columnar layouts, detecting the period from the key, and reordering columns
# with period-targeted swaps. It should recover the plaintext in
# transposition_solve_solution.txt (score ~7.07). It is stochastic — increase
# -nrestarts / -nhillclimbs for a harder cipher.
SRC="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
"$SRC/colossus" -type transposition -cipher "$HERE/transposition_solve.txt" \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
  -nrestarts 800 -nhillclimbs 5000 2>/dev/null \
  | grep -E 'Result Score|^[A-Z]{30}' | head -3
