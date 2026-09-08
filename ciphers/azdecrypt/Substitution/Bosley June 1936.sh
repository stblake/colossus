#!/bin/bash
home=../../..
$home/colossus -type aristocrat -multiline \
    -cipher "Bosley June 1936.txt" -nthreads 8 \
    -ngramsize 5 -ngramfile $home/english_quintgrams.txt -logprob \
    -nrestarts 1200 -nhillclimbs 50000 \
    -verbose \
    -spaces -spacesngramsize 5 -spacesngramfile $home/5-grams_english+spaces_jarlve_reddit_v1912.txt \
    -check-solution-file "Bosley June 1936.solution"

