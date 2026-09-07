#!/bin/bash

./colossus -type 7 -cipher porta_aca.txt -ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt -nhillclimbs 500 -nrestarts 100 -backtrackprob 0.15 -slipprob 0.05 -cyclewordlen 11 -stochasticcycle -verbose
