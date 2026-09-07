#!/bin/bash

./quagmire -type 4 -cipher cipher_quagmire_4_longer.txt -ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt -nsigmathreshold 1. -nhillclimbs 5000 -nrestarts 15000 -backtrackprob 0.15 -slipprob 0.0005 -verbose -maxcyclewordlen 12 -plaintextkeywordlen 5 -ciphertextkeywordlen 6 -cyclewordlen 6
