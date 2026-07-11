#!/bin/bash

# Z340 with its transposition already undone -> a pure homophonic substitution
# (340 symbols, 63 distinct, only ~5.4 occurrences/symbol, so weaker per-symbol
# signal than Z408). Tuned schedule (vs the defaults):
#   -inittemp 0.015  a cool anneal; Z340 is harder than Z408 and needs to settle
#                    rather than wander (defaults 0.10 and even 0.02 stall lower);
#   -weightmono 1.4  monogram (anti-collapse) weight ~1.3-1.5 makes the true
#                    plaintext the stable attractor (>=1.7 converges ~2pts worse).
# With the incremental scorer ~30M iterations land a stable ~95% solve in ~25s,
# reliably across seeds. The recovered text is plainly the Zodiac message:
#   IHOPEYOUAREHA[V]INGLOTSOF[F]UN...TRYINGTOCATCHME...GASCHAMBER...PARADICE...
# The residual ~16 characters are genuinely ambiguous homophones (V<->T, D<->L,
# K<->T classes that quintgram statistics cannot separate) -- the fitness ceiling
# for a crib-free attack. Known cribs or 6-grams would be needed to reach 100%.
../../colossus -type homophonic -multiline \
    -cipher z340untransposed.txt \
    -ngramsize 5 -ngramfile ../../english_quintgrams.txt -logprob \
    -nrestarts 600 -nhillclimbs 50000 \
    -inittemp 0.015 -weightmono 1.4 \
    -backtrackprob 0.15 -verbose
