#!/bin/bash
# End-to-end test: solve a DOUBLE columnar transposition (-type transcol2).
# transcol_double.txt was produced by two stacked columnar stages: stage 1 wrote
# the plaintext into 7 columns (read top-to-bottom, keyed), then stage 2 wrote
# that into 9 columns (read top-to-bottom, keyed). The solver randomises (K1,K2)
# over -mincols..-maxcols each restart and anneals both column orders. Double
# columnar is inherently harder and stochastic — raise -nrestarts / -nhillclimbs
# if it does not converge. It should recover transcol_solution.txt.
SRC="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
out=$("$SRC/colossus" -type transcol2 -cipher "$HERE/transcol_double.txt" \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
  -mincols 4 -maxcols 12 -nrestarts 400 -nhillclimbs 12000 2>/dev/null)
echo "$out" | grep -E 'Result Score|columnar,'
pt=$(echo "$out" | grep -A3 '^Result Score' | sed -n '3p')
if [ "$pt" = "$(cat "$HERE/transcol_solution.txt")" ]; then echo "PASS"; else echo "FAIL: $pt"; exit 1; fi
