#!/bin/bash
# End-to-end test: a columnar transposition whose columns were read BOTTOM-TO-TOP
# (transcol_single_bt.txt, 11 columns). Read-direction variants are opt-in, so the
# solver must be told with -readdir bt (or -readdir both); the default top-to-bottom
# search will NOT recover it. Confirms the variant flag works and is required.
SRC="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
sol="$(cat "$HERE/transcol_solution.txt")"

out=$("$SRC/colossus" -type transcol -cipher "$HERE/transcol_single_bt.txt" \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
  -readdir bt -nrestarts 60 -nhillclimbs 8000 2>/dev/null)
echo "$out" | grep -E 'Result Score|columnar,'
pt=$(echo "$out" | grep -A3 '^Result Score' | sed -n '3p')
[ "$pt" = "$sol" ] && echo "PASS (-readdir bt recovered)" || { echo "FAIL: $pt"; exit 1; }

# Negative control: the default top-to-bottom search should not recover a bt cipher.
out2=$("$SRC/colossus" -type transcol -cipher "$HERE/transcol_single_bt.txt" \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
  -nrestarts 60 -nhillclimbs 8000 2>/dev/null)
pt2=$(echo "$out2" | grep -A3 '^Result Score' | sed -n '3p')
[ "$pt2" != "$sol" ] && echo "PASS (default tb did not crack bt, as expected)" \
  || { echo "FAIL: default tb unexpectedly recovered a bt cipher"; exit 1; }
