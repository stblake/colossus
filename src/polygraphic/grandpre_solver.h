#ifndef GRANDPRE_SOLVER_H
#define GRANDPRE_SOLVER_H
#include "colossus.h"

// =====================================================================
//  Grandpre solver (TYPE grandpre) -- see grandpre_solver.c
// =====================================================================
//
// Grandpre decoding is code -> letter (unique), so cryptanalytically it is a HOMOPHONIC
// substitution over the <= N^2 distinct 2-digit codes -> 26 letters. The solver interns
// the observed codes into symbol ids and hill-climbs the many-to-one symbol -> letter map
// against the n-gram score on the HOMOPHONIC INCREMENTAL FAST PATH (a byte-for-byte twin
// of homophonic_solver.c: monogram-flattening seed, single-symbol reassignment moves,
// chi-squared anti-collapse penalty, incremental score_neighbor over the touched windows).
// Unlike the Checkerboard the row/col labels are the grid indices, so the recovered
// square is reported directly (not merely up to a label permutation).

void solve_grandpre(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

#endif
