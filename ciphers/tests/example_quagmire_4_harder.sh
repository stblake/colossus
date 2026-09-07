#!/bin/bash

./quagmire -type 4 -cipher cipher_quagmire_4_harder.txt -crib crib.txt -ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt -nsigmathreshold 1. -nhillclimbs 3000 -nrestarts 30000 -backtrackprob 0.15 -slipprob 0.000 -verbose -plaintextkeywordlen 5 -ciphertextkeywordlen 3 -cyclewordlen 7

