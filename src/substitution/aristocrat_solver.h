#ifndef ARISTOCRAT_SOLVER_H
#define ARISTOCRAT_SOLVER_H
#include "colossus.h"

// Aristocrat (type 79) / Patristocrat (type 80): a simple MONOALPHABETIC substitution -- one
// fixed 26-letter bijection sends each plaintext letter to a ciphertext letter. The Aristocrat
// preserves word divisions (the ciphertext keeps its spaces); the Patristocrat drops them
// (5-letter groups). Both share ONE solver core: the search is over a letters-only stream (the
// >=0 entries of the decoded ciphertext), the key is a 26-permutation climbed by n-gram score
// with the homophonic-style incremental fast path, and the two types differ ONLY in how the
// report lays out the recovered plaintext. cfg->cipher_type selects the presentation.
//
// Scoring is n-gram only (word divisions are used for the spaced report, the dictionary word
// count, and keyword recovery -- NOT as a solve gradient). Cribs are not used (the letters-only
// stream does not carry absolute cipher positions), matching the Ragbaby solver's choice.
void solve_aristocrat(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Substitution move, exposed for the solver test's move-invariant check: swap two ciphertext
// letters' plaintext images in key[0..alpha-1], keeping it a permutation of 0..alpha-1.
void aristocrat_move(int *key, int alpha);
#endif
