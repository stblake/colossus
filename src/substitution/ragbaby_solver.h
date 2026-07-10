#ifndef RAGBABY_SOLVER_H
#define RAGBABY_SOLVER_H
#include "colossus.h"

// Ragbaby (type 77): a keyed 24-letter alphabet in which each plaintext LETTER is shifted
// forward by its word-position number (mod 24). The per-letter numbering is FIXED by the
// (spaced) ciphertext, so the only unknown is the keyed alphabet; the solver anneals it as an
// ACA keyed alphabet (keyword prefix + ordered tail), scored by plaintext n-grams. Parses the
// spaced ciphertext_str directly (like the digit-stream solvers) rather than the A..Z decode.
void solve_ragbaby(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Keyed-alphabet search internals, exposed for the solver tests. The 24-cell twin of
// fracmorse_canonicalize / fracmorse_move_seq: canonicalize re-sorts the tail seq[kw..23]
// ascending after the keyword prefix changes; move_seq applies one structure-preserving move
// (keeping seq a permutation of 0..23 with seq[kw..23] strictly ascending).
void ragbaby_canonicalize(int *seq, int kw);
void ragbaby_move_seq(int *seq, int *kw);
#endif
