#!/bin/bash

./colossus -type auto1 -cipher cipher_autokey_quag1_variant.txt -ngramsize 5 -ngramfile ngram_data/english/english_quintgrams.txt -nhillclimbs 500 -nrestarts 1000 -backtrackprob 0.15 -slipprob 0.0005 -verbose -cyclewordlen 7 -variant
-plaintextkeyword KRYPTOS

