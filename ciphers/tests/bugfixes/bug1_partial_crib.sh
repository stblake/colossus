#!/bin/bash
# Bug 1: the trailing partial-crib match row used to index the positionally-packed
# crib_indices array (and test it against an uninitialized -1), printing garbage
# digits/'*' at EVERY position. The fix indexes the crib by cipher position via
# cribtext_str ('_' = no crib).
#
# bug1_cipher.txt is "ITISATRUTH...NEIGHBOURHOOD" enciphered with Vigenere key
# KRYPTOS (len 7). bug1_crib.txt reveals UNIVERSALLY, FORTUNE and NEIGHBOURHOOD.
# When solved, the partial-crib row should read all 0s (exact match) under each
# crib word and '_' everywhere else -- perfectly aligned with the crib line above.
cd "$(dirname "$0")/../../.." || exit 1
./colossus -type vig \
  -cipher ciphers/tests/bugfixes/bug1_cipher.txt \
  -crib   ciphers/tests/bugfixes/bug1_crib.txt \
  -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt \
  -nhillclimbs 300 -nrestarts 400
