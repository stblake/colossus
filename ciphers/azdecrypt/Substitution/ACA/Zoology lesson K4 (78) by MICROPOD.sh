#!/bin/bash
home=../../../..
$home/colossus -type aristocrat -multiline \
    -cipher "Zoology lesson K4 (78) by MICROPOD.txt" -nthreads 8 \
    -ngramsize 5 -ngramfile $home/english_quintgrams.txt -logprob \
    -verbose \
    -spaces -spacesngramsize 5 -spacesngramfile $home/5-grams_english+spaces_jarlve_reddit_v1912.txt \
    -check-solution-file "Zoology lesson K4 (78) by MICROPOD.solution"
