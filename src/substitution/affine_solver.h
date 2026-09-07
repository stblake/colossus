#ifndef AFFINE_SOLVER_H
#define AFFINE_SOLVER_H
#include "colossus.h"

// =====================================================================
//  Affine cipher solver (TYPE affine / af)
// =====================================================================
//
// A DETERMINISTIC-EXHAUSTIVE solve (a "needle", like Pollux / Morbit): the keyspace is only
// 12 coprime multipliers x 26 shifts = 312 keys, and one key change re-decrypts the whole
// stream (no local gradient), so we enumerate every key, decrypt, score by n-gram, and keep
// the best -- guaranteed the global optimum, no hill-climb needed. Word divisions are carried
// through cipher_indices as sentinels and restored in the report (Aristocrat-style).
void solve_affine(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);
#endif
