#!/bin/bash

# Z408 is short (408 symbols) and carries the Zodiac's own misspellings, so the
# homophonic climb has sticky near-miss basins. Tuned schedule (vs the defaults):
#   -inittemp 0.02   a cooler anneal settles into the true basin instead of
#                    wandering at high temperature (default 0.10 stalls ~1/6 seeds);
#   -weightmono 1.5  a slightly stronger monogram (anti-collapse) penalty makes the
#                    real plaintext the clear global maximum over the collapse fixed
#                    point that high-frequency-letter folding would otherwise win.
# With the incremental scorer each iteration is cheap, so ~7.5M iterations land the
# full plaintext in a couple of seconds, reliably across seeds (was 8/8 in testing).
../../colossus -type homophonic -multiline \
    -cipher z408.txt \
    -ngramsize 5 -ngramfile ../../ngram_data/english/english_quintgrams.txt -logprob  \
    -nrestarts 150 -nhillclimbs 50000 \
    -inittemp 0.02 -weightmono 1.5 \
    -backtrackprob 0.15 -seed 1 -verbose

