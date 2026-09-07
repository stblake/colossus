#!/bin/bash
# Bug 4: a cipher/ngram/crib path longer than MAX_FILENAME_LEN overflowed the
# fixed char[] in ColossusConfig (main() strcpy's the CLI arg in unbounded),
# corrupting the config struct and crashing with SIGILL ("Illegal instruction").
# This bit any absolute path past the old 100-byte limit; it was first hit by the
# transposition test scripts, which pass absolute paths.
#
# Fix: MAX_FILENAME_LEN raised to 4096. This test passes a deliberately long
# (>100 char) absolute path and asserts the solver runs and emits its '>>>'
# summary instead of crashing.
SRC="$(cd "$(dirname "$0")/../../.." && pwd)"
WORK="$(mktemp -d)"
LONGDIR="$WORK/$(printf 'd%.0s' $(seq 1 130))"   # 130-char directory component
mkdir -p "$LONGDIR"
cp "$SRC/ciphers/tests/transperoffset_solve.txt" "$LONGDIR/c.txt"
CIPHER="$LONGDIR/c.txt"
echo "cipher path length = ${#CIPHER} (was crashing when > 100)"
N=$("$SRC/colossus" -type transperoffset -cipher "$CIPHER" \
      -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" \
      -nrestarts 10 -nhillclimbs 1000 2>/dev/null | grep -c '>>>')
if [ "$N" -ge 1 ]; then echo "PASS: solver ran without crashing (found $N summary line)"; else echo "FAIL: no summary line — likely crashed"; fi
rm -rf "$WORK"
