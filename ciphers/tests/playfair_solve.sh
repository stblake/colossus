#!/bin/bash
#
# End-to-end test: solve a Playfair cipher (-type playfair) and verify the recovered
# plaintext byte-for-byte against playfair_pride.solution.
#
# The 976-character ciphertext was enciphered with the keyword KRYPTOSABCDEF over the
# 25-letter alphabet (J merged into I). Playfair is genuinely near the limit of a
# quadgram attack on short texts, so this uses a long message, a fixed seed, the
# AZDecrypt-style log-probability scoring (-logprob), and the tuned anneal schedule.
# Playfair grids are unique only up to cyclic row/column rotation, but every such grid
# decrypts identically, so the recovered PLAINTEXT is what we check.
SRC="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"

out=$("$SRC/colossus" -type playfair -cipher "$HERE/playfair_pride.txt" \
  -ngramsize 4 -ngramfile "$SRC/ngram_data/english/english_quadgrams.txt" -logprob -seed 1 \
  -nrestarts 6 -nhillclimbs 400000 -inittemp 0.08 -backtrackprob 0.3 2>/dev/null)

echo "$out" | grep -E 'Result Score|^>>> ' | cut -c1-72
pt=$(echo "$out" | grep '^>>> ' | tail -1 | awk -F', ' '{print $NF}')
expected=$(tr -d '[:space:]' < "$HERE/playfair_pride.solution")

if [ "$pt" = "$expected" ]; then
    echo "PASS"
else
    echo "FAIL"
    echo "    expected: ${expected:0:80}..."
    echo "    got:      ${pt:0:80}..."
    exit 1
fi
