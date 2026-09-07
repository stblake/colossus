#!/bin/bash
# End-to-end test: transperoffset solve with a crib that pins the offset.
#
# transperoffset_crib.txt reveals the first 10 plaintext letters (THEUNITEDS) at
# positions 0-9. The blended crib score forces the rotation that places those
# letters at the start, so the solver should report period = 7, offset = 11 and
# the plaintext beginning THEUNITEDSTATESCONSTITUTION...
SRC="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
"$SRC/colossus" -type transperoffset -cipher "$HERE/transperoffset_solve.txt" \
  -crib "$HERE/transperoffset_crib.txt" \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
  -nrestarts 300 -nhillclimbs 8000 2>/dev/null | grep -E 'transperiodoffset:|Result Score'
