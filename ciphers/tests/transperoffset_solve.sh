#!/bin/bash
# End-to-end test: solve a pure transperoffset transposition by optimization.
#
# The ciphertext decrypts (period d=7, offset n=11) to the known English plaintext
# in transperoffset_solve_solution.txt. The solver should recover period = 7; the
# offset is ambiguous up to a cyclic rotation without a crib (the n-gram score is
# rotation-invariant), so the plaintext comes out correct but possibly rotated.
# See transperoffset_crib_solve.sh for the crib-pinned offset variant.
SRC="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
"$SRC/colossus" -type transperoffset -cipher "$HERE/transperoffset_solve.txt" \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
  -nrestarts 40 -nhillclimbs 3000 2>/dev/null | grep -E 'transperiodoffset:|Result Score'
