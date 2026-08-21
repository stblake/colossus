#!/bin/bash
# End-to-end test: solve a single keyed columnar transposition with the dedicated
# columnar solver (-type transcol). The ciphertext transcol_single_tb.txt was
# produced by writing the plaintext row-wise into 9 columns and reading the
# columns out top-to-bottom in a keyed order (incomplete final row). The solver
# does not know the column count or order; it sweeps -mincols..-maxcols and
# hill-climbs the column-order permutation. It should recover transcol_solution.txt.
# Stochastic — raise -nrestarts / -nhillclimbs for a harder cipher.
SRC="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
out=$("$SRC/colossus" -type transcol -cipher "$HERE/transcol_single_tb.txt" \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
  -nrestarts 40 -nhillclimbs 8000 2>/dev/null)
echo "$out" | grep -E 'Result Score|columnar,'
pt=$(echo "$out" | grep -A3 '^Result Score' | sed -n '3p')
if [ "$pt" = "$(cat "$HERE/transcol_solution.txt")" ]; then echo "PASS"; else echo "FAIL: $pt"; exit 1; fi
