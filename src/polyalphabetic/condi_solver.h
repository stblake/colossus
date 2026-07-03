#ifndef CONDI_SOLVER_H
#define CONDI_SOLVER_H
#include "colossus.h"

// Condi solver (TYPE condi). Rides the cipher-agnostic engine. Condi is a plaintext-feedback
// cipher over a keyed alphabet sigma (a 26-permutation): the shift for each letter is the position
// of the preceding plaintext letter in sigma, the first letter using a starter offset. There is NO
// period and no per-column decoupling; the shipped attack is a free-permutation sigma anneal with
// the 26 starter values enumerated as engine configs. NOTE (see condi_solver.c / the solver test):
// the plaintext feedback makes the true sigma an ISOLATED NEEDLE (one cell swap cascades the whole
// downstream decrypt -> no basin), so local search does NOT recover Condi blind at any budget -- a
// documented structural limitation; the tractable attack is crib-anchored constraint solving.
// -startkey pins a starter. score_adjust stays 0. Cribs are passed through (positional).
void solve_condi(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

#endif
