#!/bin/bash
# Bit-identical regression oracle for the performance-optimization work.
# Runs a spread of cipher types/modes with a FIXED seed and prints the ">>>"
# summary lines. Output must be byte-for-byte identical before and after each
# optimization (that is the agreed correctness bar). Usage:
#
#   tests/perf_golden.sh > /tmp/after.txt && diff /tmp/golden.txt /tmp/after.txt
#
# Run from the repo root (the binary loads ngram/dict/cipher from cwd).
set -u
BIN=./colossus
NG="-ngramsize 4 -ngramfile english_quadgrams.txt"
# Quintgram + -logprob row (below) exercises the memory-bound big-table scoring path and
# the g_alpha==26 quintgram fast path (with its software prefetch) -- keep it here so any
# scorer optimisation is verified bit-identical on that path, not just the quadgram one.
QNG="-ngramsize 5 -ngramfile english_quintgrams.txt -logprob"
run() { $BIN "$@" -seed 42 2>/dev/null | grep '^>>>'; }

echo "Q3   :";    run -type quag3 -cipher ciphers/kryptos/K1.txt $NG -nhillclimbs 300 -nrestarts 60 -backtrackprob 0.15
echo "Q3+crib:";  run -type quag3 -cipher cipher.txt -crib crib.txt $NG -keywordlen 7 -cyclewordlen 7 -nhillclimbs 300 -nrestarts 60 -backtrackprob 0.15
echo "Vig  :";    run -type vig   -cipher ciphers/kryptos/K1.txt $NG -nhillclimbs 300 -nrestarts 60
echo "Beau :";    run -type beau  -cipher ciphers/kryptos/K1.txt $NG -nhillclimbs 300 -nrestarts 60
echo "Porta:";    run -type porta -cipher ciphers/kryptos/K1.txt $NG -nhillclimbs 300 -nrestarts 60
echo "Auto3:";    run -type auto3 -cipher ciphers/kryptos/K1.txt $NG -nhillclimbs 300 -nrestarts 30 -maxcyclewordlen 12
echo "Q3stoch:";  run -type quag3 -cipher ciphers/kryptos/K1.txt $NG -nhillclimbs 300 -nrestarts 60 -stochasticcycle
echo "Bifid5:";   run -type bifid -cipher ciphers/tests/bifid_pride.txt $QNG -period 7 -nrestarts 1 -nhillclimbs 120000 -inittemp 0.08 -backtrackprob 0.3 -nthreads 1
