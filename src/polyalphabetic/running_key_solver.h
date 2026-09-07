#ifndef RUNNING_KEY_SOLVER_H
#define RUNNING_KEY_SOLVER_H
#include "colossus.h"

//
// Running Key solver (TYPE running-key / runningkey / rk).
//
// Three modes, one CipherModel:
//   * ACA self-keyed blind (default): the plaintext is one passage of length 2N, its
//     first half is the running key, its second half is the plaintext; recover the
//     whole passage. Both halves must read as English (that joint constraint is what
//     pins the otherwise-underdetermined solution) plus a seam reward at the K|P join.
//   * independent-key blind (-indepkey): the key is an unrelated English text; recover
//     both streams (no seam, report the plaintext half).
//   * known-key (-runningkeyfile <file>): deterministic decrypt against a supplied key
//     text (the "drag K1/K2/K3 as a running key" capability); family swept, best kept.
//
// The four letter-keyed families (Vigenere/Beaufort/Variant/Porta) are swept and the
// n-gram picks the winner, unless pinned by -variant / -beaufort.
//
void solve_running_key(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Deterministic beam warm start (exposed for the solver tests): recover the running-key
// stream for one family by a left-to-right beam that jointly maximises the n-gram
// fitness of BOTH the key stream and the derived plaintext. Writes N key letters into
// key_out and returns the joint beam score (raw window sum, un-normalised).
double rk_beam_warmstart(const int cipher[], int n, int family,
    const float *ngram, int ngram_size, int beam_width, int key_out[]);

#endif
