#!/bin/bash

for i in {1..100}; do
../quagmire -type 4 -cipher q4_p130.txt -ngramsize 5 -ngramfile ../ngram_data/english/english_quintgrams.txt -nsigmathreshold 1. -nhillclimbs 10000 -nrestarts 10000 -backtrackprob 0.15 -slipprob 0.0005 -plaintextkeywordlen 4 -ciphertextkeywordlen 7 -cyclewordlen 5 -verbose
done
