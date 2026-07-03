#ifndef FRACMORSE_SOLVER_H
#define FRACMORSE_SOLVER_H
#include "colossus.h"

void solve_fracmorse(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Keyed-alphabet search internals, exposed for the solver tests. fracmorse_canonicalize
// re-sorts the tail seq[kw..25] (ascending) after the keyword prefix changes;
// fracmorse_move_seq applies one structure-preserving neighbour move to the 26-cell keyed
// alphabet (keyword prefix length *kw), updating *kw on a length move. Both preserve the
// invariant "seq is a permutation of 0..25 with seq[kw..25] strictly ascending".
void fracmorse_canonicalize(int *seq, int kw);
void fracmorse_move_seq(int *seq, int *kw);
#endif
