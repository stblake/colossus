#ifndef CHECKERBOARD_SOLVER_H
#define CHECKERBOARD_SOLVER_H
#include "colossus.h"

// Checkerboard solver (checkerboard_solver.c). Recovers the keyed 5x5 square from a Checkerboard
// ciphertext (plaintext letter -> (row label, column label) digraph). The case is auto-detected
// per axis from the ciphertext: an axis showing more than `side` distinct labels is COMPLEX (two
// labels per line, homophonic), else SIMPLE (one label). Label ORDER is not identifiable
// ciphertext-only (absorbed by a row/column permutation of the square), so:
//   SIMPLE  -> the observed labels fix the merged code space; the search is a free 25<->25
//              bijection (an Aristocrat over the 25 merged codes) run on the homophonic
//              incremental fast path.
//   COMPLEX -> a square-independent per-axis chi-squared homogeneity pre-pass ranks the label
//              PAIRINGS (which 2 labels share a line); the top-K per axis are crossed into engine
//              configs, and each config's merged 25-code stream is searched as above (generic path).
// Complex recovery needs far more text than the ACA 60-90 range -- see the module header.
void solve_checkerboard(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Testing/calibration hook: the RANK (0 = best) of the TRUE label pairing among all enumerated
// pairings for one axis (0 = row, 1 = column), under the chi-squared homogeneity pre-pass over the
// letters-only digraph stream cipher[0..clen-1]. true_line[label letter] gives that label's true
// line (0..side-1). Rank 0 means the pre-pass ranks the truth first (the complex case is then
// solvable at that length). Deterministic; exposed for test_checkerboard_solver.c.
int checkerboard_pairing_rank(const int cipher[], int clen, int axis, const int true_line[]);

#endif
