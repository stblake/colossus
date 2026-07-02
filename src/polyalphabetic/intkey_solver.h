#ifndef INTKEY_SOLVER_H
#define INTKEY_SOLVER_H
#include "colossus.h"

// Interrupted Key solver (TYPE intkey / intkey-var / intkey-beau). Rides the cipher-agnostic
// engine. The key is a P-letter keyword (per-column base shifts 0..25, in the cycleword lane)
// whose key index RESETS to 0 at break points. IoC period estimation is useless through the
// interruption (columns are not monoalphabetic), so the PERIOD is swept.
//
// How break points arise is the interruption STRATEGY, enumerated by the solver (n-gram score
// picks the true one):
//   - IK_STRAT_CT   (blind): reset after a chosen CIPHERTEXT letter -> break positions known from
//                    the ciphertext alone (keyword-independent), so per-column monogram fit over
//                    the known columns recovers the keyword like a plain Vigenere. Enumerated over
//                    the 26 letters x swept period, cheaply PRE-SCORED (monogram-derived keyword,
//                    n-gram of the decrypt), only the top-K annealed (Gromark-style).
//   - IK_STRAT_PT   (blind): reset after a chosen PLAINTEXT letter -> breaks causal (keyword-
//                    dependent), warm-started by an EM loop (fit keyword, causal-decrypt, refit).
//                    Also pre-scored + top-K, competing with CT on the same n-gram scale.
//   - IK_STRAT_BREAKS (reliable, -breaks): user-supplied group-start positions (random / word-
//                    division). Mask known -> decouples like Vigenere. One config.
//   - IK_STRAT_JOINT  (blind, -intscheme joint): keyword + an N-bit break-mask searched jointly
//                    in the engine (mask in the key lane). Characterized, not guaranteed.
//
// Default blind enumeration is CT + PT (fast, reliable). BREAKS is used iff -breaks is given;
// JOINT iff -intscheme joint. -interruptor / -intscheme / -period pin the respective axis. Cribs
// are supported (positional: decrypted[i] is plaintext[i]). The three base types (Vigenere /
// Variant / Beaufort) share the solver, branched on cfg->cipher_type. See intkey_solver.c.
void solve_intkey(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

#endif
