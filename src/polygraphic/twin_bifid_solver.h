#ifndef TWIN_BIFID_SOLVER_H
#define TWIN_BIFID_SOLVER_H
#include "colossus.h"

// Twin Bifid (ACA): two Bifid messages sharing one keyed Polybius square at different
// periods. The primary ciphertext (message 1) arrives as cipher_indices/cipher_len; the
// SECOND ciphertext is supplied via cfg->twincipher_str (-cipher2), decoded here. The
// solver anneals a single shared square, scoring the two decrypts jointly.
void solve_twin_bifid(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

#endif
