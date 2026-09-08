#!/bin/bash
home=../../..
$home/colossus -type homophonic -multiline \
    -cipher "_pi 1.txt" -nthreads 8 \
    -ngramsize 5 -ngramfile $home/english_quintgrams.txt -logprob  \
    -nrestarts 1500 -nhillclimbs 50000 \
    -inittemp 0.02 -weightmono 1.5 \
    -backtrackprob 0.15 -verbose \
    -check-solution-file "_pi 1.solution" \
    -spaces -spacesngramsize 5 -spacesngramfile $home/5-grams_english+spaces_jarlve_reddit_v1912.txt

