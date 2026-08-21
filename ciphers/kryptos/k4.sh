#!/bin/bash

./colossus -type vig -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type beaufort -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type porta -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type quag1 -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type quag2 -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type quag3 -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type quag4 -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

# variants

./colossus -type vig -variant -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type beaufort -variant -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type porta -variant -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type quag1 -variant -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type quag2 -variant -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type quag3 -variant -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose

./colossus -type quag4 -variant -cipher k4.txt -crib crib.txt -ngramsize 4 -ngramfile ngram_data/english/english_quadgrams.txt -nhillclimbs 500 -nrestarts 10000 -backtrackprob 0.15 -verbose


# autokey





