#ifndef POLLUX_SOLVER_H
#define POLLUX_SOLVER_H
#include "colossus.h"

// =====================================================================
//  Pollux solver (TYPE pollux) -- deterministic EXHAUSTIVE search
// =====================================================================
//
//  A Pollux key is a map digit 0..9 -> {DOT, DASH, X}, so the ENTIRE keyspace is
//  3^10 = 59,049 -- small enough to enumerate exhaustively in a fraction of a
//  second. And the landscape is a needle (re-assigning one digit re-parses the whole
//  Morse stream, so there is no local gradient), which is exactly the case where a
//  hill climber flails but enumeration is guaranteed. So, like the Period Column
//  solver, Pollux is a dedicated deterministic solver, NOT a hill-climbing
//  CipherModel: it decodes + scores every key and keeps the global best. There is no
//  engine schedule / -method to tune (the only knob is the validity weight below,
//  calibrated by the capability test).
//
//  Fitness is the shared n-gram score of the decoded plaintext (length-fair: it is a
//  MEAN log-prob) plus a Morse-VALIDITY reward (fraction of tokens that are legal
//  codewords), the analogue of Fractionated Morse's validity term -- it breaks near
//  ties and suppresses degenerate keys (e.g. mostly-X maps that decode to a few noisy
//  letters).

// Default weight of the Morse-validity reward folded into the score. Tuned against
// tests/test_pollux_solver.c (the reward-only quadgram table).
#define POLLUX_VALID_WEIGHT 3.0

// Letter emitted for an invalid Morse token ('X'), like Fractionated Morse.
#define POLLUX_FILLER 23

// Exhaustive search: decode digits[0..clen-1] under every 3^10 key, score each with
// state_score(ngram) + valid_weight * (n_valid/n_tokens), and return the best score,
// writing the winning key -> best_key[0..9] (each PX_DOT/PX_DASH/PX_X), the decoded
// plaintext -> best_pt[0..*best_n-1], and its token / valid counts. Exposed for the
// solver tests (which sweep valid_weight and length).
double pollux_search(const int *digits, int clen,
    float *ngram_data, int ngram_size, float weight_ngram, double valid_weight,
    int best_key[10], int best_pt[], int *best_n, int *best_nt, int *best_nv);

// Full solver entry, dispatched from solve_cipher().
void solve_pollux(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);
#endif
