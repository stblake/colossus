// Public API for the Enigma solver (see enigma_solver.c / enigma_bombe.c).

#ifndef ENIGMA_SOLVER_H
#define ENIGMA_SOLVER_H

#include "colossus.h"
#include "enigma.h"

// Ciphertext-only (Gillogly IoC + ring refinement + plugboard climb) OR, when cfg->enigma
// pins enough of the key, a deterministic known-key decrypt. Dispatches to the Bombe when
// cfg->enigma_bombe is set and a crib is present.
void solve_enigma(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Turing-Welchman menu/diagonal-board known-plaintext attack (enigma_bombe.c). Recovers
// rotor config + plugboard from a crib. Returns true if a solution was reported.
bool solve_enigma_bombe(ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

// Shared report helper (used by the solver, the known-key path, and the Bombe): prints the
// human block + the ">>>" CSV summary and fills *result (may be NULL). `decrypted` is the
// recovered plaintext for `key`.
void enigma_emit_report(ColossusConfig *cfg, SharedData *shared,
    int cipher[], int cipher_len, const EnigmaKey *key, int decrypted[], char *cribtext,
    double score, SolveResult *result);

// Enumerate the wheel-order templates the search will try (a pinned -rotors => 1; otherwise
// all 60 ordered triples of the default pool I..V). Each template carries reflector,
// n_wheels, rotor ids, and (M4) the Greek. Returns the count. Shared by the ciphertext-only
// phase 1 and the Bombe.
int enigma_enumerate_wheel_orders(const ColossusConfig *cfg, EnigmaKey *tmpl, int cap);

// Ring refinement + plugboard climb over a list of candidate base keys (wheel order + start
// positions, rings AAA). Shared by the ciphertext-only path and the Bombe. Honours cfg
// ring/plugboard pins, reports, and fills *result.
void enigma_attack_from_bases(ColossusConfig *cfg, SharedData *shared,
    int cipher[], int cipher_len, char *cribtext,
    int crib_indices[], int crib_positions[], int n_cribs,
    const EnigmaKey bases[], int n_bases, int maxplugs, SolveResult *result);

// Greedy plugboard recovery (Gillogly): from the fixed machine `key` (its plugboard is
// ignored/overwritten), greedily add the best-improving stecker pair by n-gram fitness up
// to `maxplugs`, writing the recovered involution into out_plug[0..25]. Returns the final
// n-gram score. Reused by the ciphertext-only climb's warm start and the Bombe's stop
// completion.
double enigma_greedy_plugboard(const EnigmaKey *key, int cipher[], int cipher_len,
    float *ngram_data, int ngram_size, int maxplugs, int out_plug[26]);

#endif // ENIGMA_SOLVER_H
