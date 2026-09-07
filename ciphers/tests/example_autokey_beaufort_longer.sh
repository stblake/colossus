#!/bin/bash

./colossus -type autobeau -cipher example_autokey_beaufort_longer.txt -ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt -nhillclimbs 500 -nrestarts 1000 -backtrackprob 0.15 -slipprob 0.0005 -verbose -cyclewordlen 14

