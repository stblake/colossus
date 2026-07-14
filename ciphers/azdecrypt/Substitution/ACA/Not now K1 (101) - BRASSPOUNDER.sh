#!/bin/bash
home=../../../..
$home/colossus -type aristocrat -multiline \
    -cipher "Not now K1 (101) - BRASSPOUNDER.txt" -nthreads 8 \
    -ngramsize 5 -ngramfile $home/english_quintgrams.txt -logprob \
    -verbose \
    -spaces -spacesngramsize 5 -spacesngramfile $home/5-grams_english+spaces_jarlve_reddit_v1912.txt \
    -check-solution-file "Not now K1 (101) - BRASSPOUNDER.solution"
