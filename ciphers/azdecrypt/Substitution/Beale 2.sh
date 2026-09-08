#!/bin/bash
home=../../..
$home/colossus -type homophonic \
    -cipher "Beale 2.txt" -nthreads 8 \
    -ngramsize 5 -ngramfile $home/english_quintgrams.txt -logprob \
    -nrestarts 1200 -nhillclimbs 1500000 \
    -weightentropy 1.5 -weightmono 0 -inittemp 0.8 -mintemp 0.005  \
    -backtrackprob 0.20 \
    -verbose \
    -delimiter , \
    -spaces -spacesngramsize 5 -spacesngramfile $home/5-grams_english+spaces_jarlve_reddit_v1912.txt \
    -check-solution-file "Beale 2.solution"
