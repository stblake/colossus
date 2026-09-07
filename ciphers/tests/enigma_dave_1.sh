# Dave's Enigma test cipher.  Uses same key settings as Gillogly's example.

#!/bin/bash
home=../..
$home/colossus -type enigma \
    -cipher $home/ciphers/tests/enigma_dave_1.txt \
    -ngramsize 6 \
    -ngramfile $home/ngram_data/english/english_sixgrams.txt \
    -logprob \
    -nthreads 8 \
    -check-solution-file $home/ciphers/tests/enigma_dave_1.solution \
    -spaces -spacesngramsize 4 \
    -spacesngramfile $home/ngram_data/english/english_spaces_quadgrams.txt \
    -bombe -crib $home/ciphers/tests/enigma_dave_1.crib
