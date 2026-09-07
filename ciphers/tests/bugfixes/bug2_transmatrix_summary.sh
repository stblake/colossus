#!/bin/bash
# Bug 2: the -transmatrix '>>>' one-line summary, in the NO-dictionary branch,
# printed trans_period/trans_offset (the -transperoffset fields, here stale)
# instead of trans_w1/trans_w2/trans_clockwise. The dictionary branch was fine.
#
# The no-dictionary branch only fires when OxfordEnglishWords.txt is NOT auto-
# loaded, so this runs the binary from a scratch dir that lacks it. We pass
# -transmatrix 7 11 cw; the summary must now show "7, 11, 1" (w1, w2, cw),
# not the old stale "<period>, <offset>".
SRC="$(cd "$(dirname "$0")/../../.." && pwd)"
WORK="$(mktemp -d)"
cp "$SRC/ciphers/tests/bugfixes/bug2_cipher.txt" "$WORK/cipher.txt"
cd "$WORK" || exit 1
echo "Running -transmatrix 7 11 cw with no dictionary; expect '>>> <score>, <type>, 7, 11, 1, ...'"
"$SRC/colossus" -type vig -cipher cipher.txt \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
  -transmatrix 7 11 cw -nhillclimbs 20 -nrestarts 10 2>/dev/null | grep '>>>' | cut -c1-45
rm -rf "$WORK"
