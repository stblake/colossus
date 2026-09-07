#!/bin/bash
# End-to-end test: solve a pure transmatrix (double grid rotation) transposition
# by optimization.
#
# The ciphertext decrypts (w1=9, w2=11, direction=cw) to the known English
# plaintext in transmatrix_solve_solution.txt. The solver should recover
# w1 = 9, w2 = 11, direction = cw and the exact plaintext. The (w1,w2,dir) space
# is larger than transperoffset's, so this uses more restarts.
SRC="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
"$SRC/colossus" -type transmatrix -cipher "$HERE/transmatrix_solve.txt" \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
  -nrestarts 200 -nhillclimbs 3000 2>/dev/null | grep -E 'transmatrix:|Result Score'
