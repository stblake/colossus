#!/bin/bash

../quagmire -type 4 -cipher q4_p131.txt -ngramsize 5 -ngramfile ../ngram_data/english/english_quintgrams.txt -nsigmathreshold 1. -nhillclimbs 10000 -nrestarts 10000 -backtrackprob 0.15 -slipprob 0.0005 -plaintextkeywordlen 5 -ciphertextkeywordlen 5 -cyclewordlen 5 -verbose

