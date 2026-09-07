#!/bin/bash

../quagmire -type 2 -cipher cipher_p133.txt -ngramsize 5 -ngramfile ../ngram_data/english/english_quintgrams.txt -nsigmathreshold 1. -nhillclimbs 10000 -nrestarts 10000 -backtrackprob 0.15 -slipprob 0.0005 -plaintextkeywordlen 7 -ciphertextkeywordlen 7 -cyclewordlen 4 -verbose

