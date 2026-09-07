#!/bin/bash

for i in {1..1000}; do
../quagmire -type 3 -cipher q3_p127.txt -ngramsize 5 -ngramfile ../ngram_data/english/english_quintgrams.txt -nsigmathreshold 1. -nhillclimbs 5000 -nrestarts 10000 -backtrackprob 0.15 -slipprob 0.0005 -plaintextkeywordlen 8 -ciphertextkeywordlen 8 -cyclewordlen 8 -verbose
done
