#!/bin/bash

./colossus -type 7 -cipher autokey_len97_wl21.txt -ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt -nhillclimbs 500 -nrestarts 1000 -backtrackprob 0.15 -slipprob 0.0005 -cyclewordlen 21 -verbose

