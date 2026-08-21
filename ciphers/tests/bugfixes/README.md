# Bug-fix demonstrations

End-to-end tests for the three "Known issues" that were listed in `CLAUDE.md`.
Each script builds/uses the in-tree `./colossus` and shows the behaviour
after the fix. Run them from anywhere; they resolve the repo root themselves.
Build first: `make` (or `gcc -Wall -O3 -o colossus *.c`).

## Bug 1 — partial-crib match row (`bug1_partial_crib.sh`)

The trailing `_`/digit/`*` row printed under the crib indexed the positionally
*packed* `crib_indices` array and tested it against an uninitialized `-1`, so it
emitted digits/`*` at every position regardless of where cribs were. Fixed to
index the crib by cipher position through `cribtext_str` (`'_'` = no crib).

`bug1_cipher.txt` is an English plaintext enciphered with Vigenere key `KRYPTOS`;
`bug1_crib.txt` reveals `UNIVERSALLY`, `FORTUNE`, `NEIGHBOURHOOD`. After a correct
solve the row reads all `0` (exact match) under each crib word and `_` elsewhere,
aligned with the crib line above it. Before the fix the row was garbage.

## Bug 2 — `-transmatrix` summary, no-dictionary branch (`bug2_transmatrix_summary.sh`)

The `>>>` CSV summary's no-dictionary branch printed `trans_period`/`trans_offset`
(the unrelated `-transperoffset` fields) instead of `w1`/`w2`/`clockwise`. The
dictionary branch was already correct. The script runs from a scratch dir (so the
Oxford dictionary is not auto-loaded) with `-transmatrix 7 11 cw`; the summary now
shows `7, 11, 1` instead of the old stale values.

## Bug 3 — `load_ngrams` `while(!feof)` (`bug3_ngram_load.sh`)

The loader looped on `!feof`, re-reading the final line; on a trailing/malformed
line `fscanf` left `freq` stale and assigned it to a spurious ngram index,
corrupting the table and the normalization total (which scales every score).
Fixed to loop on `fscanf(...) == 2`.

`bug3_ng_test.c` replicates the loader both ways. On `bug3_quad_malformed.txt`
(final line `ZZZZ` with no count) the old loader injects a bogus `ZZZZ` entry and
inflates the total; the new loader ignores it. On the well-formed file the two are
identical, confirming the fix is a no-op for valid input (e.g. the real
`ngram_data/english/english_quadgrams.txt`).
