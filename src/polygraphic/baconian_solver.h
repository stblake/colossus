#ifndef BACONIAN_SOLVER_H
#define BACONIAN_SOLVER_H
#include "colossus.h"

// Baconian solver (TYPE baconian). The plaintext letters are biliteral-encoded (five
// a/b symbols each, fixed 24-letter table) and then CONCEALED in normal English cover
// text: a hidden per-symbol CLASSIFIER labels each cover unit a or b. The solver
// searches that classifier -- a 26-letter a/b labelling -- over two grouping modes
// (per-letter / per-word), decodes via the fixed table, and scores the recovered
// plaintext (n-gram) with a biliteral-VALIDITY reward folded in. See baconian_solver.c
// for the full design; blind recovery at the ACA <=25-letter maximum is gaming-limited
// (characterized in tests/test_baconian_solver.c), the reliable case being a canonical
// classifier (caught exactly by the sweep configs).
void solve_baconian(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Tuning constants (calibrated against tests/test_baconian_solver.c).
#define BAC_VALID_WEIGHT 2.0    // biliteral-validity reward weight (fraction of valid groups)
#define BAC_FILLER       23     // letter 'X' index emitted for an invalid (>=24) group

#endif
