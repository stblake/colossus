#ifndef SYLLABARY_SOLVER_H
#define SYLLABARY_SOLVER_H
#include "colossus.h"

// =====================================================================
//  Syllabary solver (TYPE syllabary) -- see syllabary_solver.c
// =====================================================================
//
// Syllabary decoding is a UNIQUE BIJECTION code -> token over the fixed 100-element syllabary
// alphabet, so cryptanalytically it is a SUBSTITUTION over 100 numeric codes -> 100 known
// multi-character tokens, with a LENGTH-CHANGING decode (tokens are 1-3 letters). The solver
// interns the observed codes and hill-climbs the code->token permutation (SHAPE_ANNEAL,
// two-code swap moves) against the n-gram score of the concatenated letters, TILED into a
// fixed scoring buffer so the mean n-gram is length-fair (the Fractionated Morse pattern). The
// row/column label order folds into the map, so this single composite-map search subsumes all
// three ACA "Known/Unknown Coordinates x Known/Unknown Keysquare" variants.
//
// CAPABILITY. A 100-token bijection over the ACA 110-154 code range is severely undersampled
// (~1-1.5 samples/code) and the isolog spellings flatten frequency, so blind recovery sits
// BELOW the ACA range -- a documented limitation (like the Checkerboard complex case /
// Tridigital), characterized (a length/accuracy curve) not asserted at 99%. Best with -logprob
// and quintgrams.

void solve_syllabary(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

#endif
