#!/bin/bash

../quagmire -type 2 -cipher q2_p125_2.txt -ngramsize 5 -ngramfile ../ngram_data/english/english_quintgrams.txt -nsigmathreshold 1. -nhillclimbs 2500 -nrestarts 10000 -backtrackprob 0.15 -slipprob 0.0005 -plaintextkeywordlen 6 -ciphertextkeywordlen 6 -verbose
