#!/bin/bash

# Beale cipher 2 is a book cipher on the U.S. Declaration of Independence (each number
# indexes a word; the plaintext letter is that word's first letter) -- which is exactly
# a HOMOPHONIC substitution: every number is one of the homophones for its plaintext
# letter. Fed to the solver as symbols via -delimiter ',' (762 positions, 182 distinct
# symbols, ~4.2 reps each -- far less over-determined than Zodiac Z408's 54 symbols).
#
# Tuned schedule (vs the defaults, and vs Z408):
#   -weightentropy 1.5  AZDecrypt's multiplicative fitness ngram*H^w (H = Shannon entropy
#                       of the decrypt). With so many symbols the true plaintext beats the
#                       homophonic "collapse" fixed point (fold many symbols onto E/T/A) by
#                       only ~0.06 on raw quintgrams -- the search just falls into collapse.
#                       Scaling the n-gram score by entropy^w suppresses the low-entropy
#                       collapse PROPORTIONALLY and makes the true basin dominant. This is
#                       the opposite regime to Z408, which instead needs the ADDITIVE
#                       monogram penalty (-weightmono 1.5); here -weightmono is turned OFF.
#   -inittemp 0.8       the entropy score scale is O(10-50) (not the ~-5 of mean-log-prob),
#                       so the anneal needs an O(1) starting temperature, not ~0.02.
#   long slow anneal    fewer, MUCH longer restarts (1200 x 1.5M climbs, cooling to 0.005):
#                       the true basin is a narrow needle at 182 symbols -- many short
#                       restarts freeze ~2.6 below it in a ~50%-correct near-miss, EVERY
#                       seed; a long slow cool escapes it and lands the true basin.
#
# ~8.5 min on 8 threads. Recovers the full treasure-contents message blind (no cribs):
# "I HAVE DEPOSITED IN THE COUNTY OF BEDFORD ABOUT FOUR MILES FROM BUFORD'S ...". A few
# letters stay wrong -- ambiguous rare-symbol homophones plus Beale's own documented
# transcription errors -- but the message reads end to end.
../../colossus -type homophonic -cipher Beale_2.txt -delimiter ',' \
    -ngramsize 5 -ngramfile ../../english_quintgrams.txt -logprob \
    -weightentropy 1.5 -weightmono 0 \
    -inittemp 0.8 -mintemp 0.005 \
    -nrestarts 1200 -nhillclimbs 1500000 \
    -backtrackprob 0.20 -nthreads 8 -seed 1 -verbose
