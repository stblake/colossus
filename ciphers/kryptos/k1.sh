#!/bin/bash

./colossus -type quag3 -cipher k1.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose


