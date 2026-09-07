#!/bin/bash

./colossus -type 0 -cipher kryptos_model_cipher.txt -ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt -nhillclimbs 500 -nrestarts 100 -backtrackprob 0.15 -slipprob 0.0005 -verbose

