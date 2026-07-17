#ifndef SEQUENCE_TRANSPOSITION_SOLVER_H
#define SEQUENCE_TRANSPOSITION_SOLVER_H
#include "colossus.h"

// Sequence Transposition solver (TYPE sequence-transposition). Recovers the 10-bucket read-order
// permutation (a small transposition climb) with the chain-addition primer either supplied
// (cfg->seq_primer, the ACA convention transmits it) or recovered by a blind primer pre-pass.
void solve_sequence_transposition(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Blind primer pre-pass (exposed for the solver tests): enumerate the 10^primer_len primer
// space, rank each primer by a quick best-read-order mini-solve, and write the top-K primers
// (descending by score, each stride SEQ_TRANS_MAX_PRIMER ints) into out_primers. Returns the
// number kept (<= K). out_primers must hold K * SEQ_TRANS_MAX_PRIMER ints.
int seq_trans_rank_primers(const int cipher[], int n, int primer_len,
    const float *ngram, int ngram_size, int K, int out_primers[]);

#endif
