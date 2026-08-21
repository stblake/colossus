#!/bin/bash

./colossus -type quag3 -cipher cipher.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -keywordlen 7 -cyclewordlen 7 -nhillclimbs 500 -nrestarts 1000 -backtrackprob 0.15 -verbose

