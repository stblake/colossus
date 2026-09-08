// Public API for the Chaocipher solver (see chaocipher_solver.c).

#ifndef CHAOCIPHER_SOLVER_H
#define CHAOCIPHER_SOLVER_H

#include "colossus.h"
#include "chaocipher.h"

// Recover the two STARTING alphabets (the key) by annealing them under the shared engine.
// Two regimes, selected by whether known plaintext is supplied:
//   * BLIND ciphertext-only (no cribs): n-gram fitness over the decrypt (SHAPE_ANNEAL).
//   * KNOWN-PLAINTEXT / crib (cribs present): folds Lasry's AGGREGATE DISPLACEMENT ERROR
//     into score_adjust as a smooth partial-credit gradient (reliable from ~50 chars).
// Fixed-position -crib and global -cribdrag also blend through state_score.
void solve_chaocipher(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Perturbation move on a candidate state (a swap-dominant move on the two 26-letter starting
// alphabets held in ct_keyword[]/pt_keyword[]). Exposed for the solver move-invariant test.
void chaocipher_perturb_state(SolverState *st);

#endif // CHAOCIPHER_SOLVER_H
