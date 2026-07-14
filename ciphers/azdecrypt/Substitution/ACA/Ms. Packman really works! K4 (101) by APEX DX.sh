#!/bin/bash
home=../../../..
$home/colossus -type aristocrat -multiline \
    -cipher "Ms. Packman really works! K4 (101) by APEX DX.txt" -nthreads 8 \
    -ngramsize 5 -ngramfile $home/english_quintgrams.txt -logprob \
    -verbose \
    -spaces -spacesngramsize 5 -spacesngramfile $home/5-grams_english+spaces_jarlve_reddit_v1912.txt \
    -check-solution-file "Ms. Packman really works! K4 (101) by APEX DX.solution"
