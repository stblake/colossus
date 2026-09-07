#!/bin/bash

./quagmire -type 1 -cipher cipher_quagmire_1_longer.txt -ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt -nsigmathreshold 1. -nhillclimbs 2500 -nrestarts 15000 -backtrackprob 0.25 -slipprob 0.0005 -verbose -plaintextkeywordlen 5 -cyclewordlen 7

