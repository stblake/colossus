#ifndef PERIOD_COLUMN_SOLVER_H
#define PERIOD_COLUMN_SOLVER_H
#include "colossus.h"

// =====================================================================
//  Period column order transposition solver (PERIOD_COLUMN)
// =====================================================================
//
//  A dedicated, DETERMINISTIC solver for AZdecrypt's "Period column order"
//  transposition and short COMPOSITIONS of it. Because a single stage's key is
//  tiny (a complete-grid width dx | len, a period p in [2, dx-1], and a
//  transpose/untranspose flag), the whole single-stage keyspace is enumerable,
//  and a composition of two stages is enumerable as the product -- so instead of
//  hill-climbing a needle (a wrong stage yields no n-gram gradient), the solver
//  exhaustively tries every stage (depth 1) and every ordered pair of stages
//  (depth 2), scoring each decrypt with the shared n-gram (+ crib) fitness and
//  keeping the global best. This reproduces AZdecrypt's stacked solutions (the
//  motivating 168-char cipher decrypts as 4x42 TP P:3 then 56x3 UTP P:2) and,
//  unlike a stochastic search, is guaranteed to find the optimum inside the
//  searched (width, period, depth<=2) space.

// One recovered / applied stage.
typedef struct { int dx, period, utp; } PCStage;

// Exhaustive nested search. Applies up to `max_depth` (1 or 2) period-column
// stages to cipher[0..len-1], returns the best n-gram(+crib) score and writes the
// best decrypt to best_decrypted[0..len-1] and the winning stages to
// stages[0..*n_stages-1]. Widths are the complete-grid divisors dx | len in
// [3, len]. weight_ngram/weight_crib pick the state_score blend (crib blend only
// engages when n_cribs > 0). ngram_data/ngram_size are the loaded n-gram table.
double period_column_search(const int *cipher, int len, int max_depth,
    float *ngram_data, int ngram_size,
    int crib_indices[], int crib_positions[], int n_cribs,
    float weight_ngram, float weight_crib,
    int best_decrypted[], PCStage stages[], int *n_stages);

// Full solver entry, dispatched from solve_cipher().
void solve_period_column(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs);
#endif
