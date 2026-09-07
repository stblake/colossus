#!/bin/bash

./quagmire -type 2 -cipher cipher_quagmire_2_longer.txt -ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt -nsigmathreshold 1. -nhillclimbs 2500 -nrestarts 15000 -backtrackprob 0.25 -slipprob 0.0005 -verbose -ciphertextkeywordlen 6 -cyclewordlen 7

