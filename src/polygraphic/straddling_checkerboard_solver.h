#ifndef STRADDLING_CHECKERBOARD_SOLVER_H
#define STRADDLING_CHECKERBOARD_SOLVER_H
#include "colossus.h"
#include "straddling_checkerboard.h"

// =====================================================================
//  Straddling Checkerboard solver (TYPE straddling-checkerboard / sc)
// =====================================================================
//
// A CipherModel anneal (SHAPE_ANNEAL), the digit-output cousin of Fractionated Morse:
// variable-length decode -> cyclic-tile to a fixed length -> a token-validity reward
// folded into score_adjust. Two structural additions handle a straddling board:
//
//   * A 45-config STRUCTURAL SWEEP over the row-indicator pair {r0,r1}. The digit-stream
//     tokenization depends ONLY on which two digit-values are the indicators (a digit equal
//     to an indicator starts a 2-digit token), so there are exactly C(10,2)=45 token
//     structures. No cheap statistic (IoC / monogram-greedy) ranks them -- both are
//     anti-discriminative -- so a real per-config MINI-SOLVE pre-pass (a short substitution
//     hill-climb) ranks the 45 and keeps the top-K (STRADDLING_KEEP, or -nprimers), each
//     warm-started from its mini-solve map. Wrong indicator pairs mis-segment and cannot
//     reach an English score.
//
//   * The annealed state is the FREE code->cell substitution `map[28]`: a bijection from
//     the 28 token codes (8 single-digit + 20 double-digit) to the 28 board cells (26
//     letters + a FIG figure-shift marker + a NULL). Keyed column labels fold into this map
//     (they are not separately identifiable), and the figure-shift digit is derivable from
//     the token's own digit -- so this ONE substitution captures the whole general cipher
//     (keyed labels + figure-shift) with no redundant lane. Cell-swap anneal (à la the
//     simple-substitution / square solvers).
//
// The ciphertext DIGIT stream is parsed from the raw string (not the A-Z decode). The
// plaintext alphabet is the 36-symbol ADFGVX set (A-Z + 0-9) so figure-shifted digits are
// recovered and scored (at negligible weight). Effectively needs -logprob (it stops the
// reward-only table preferring a common-quadgram near-miss over the true rare-word text).
// Cribs are unused (the length change breaks the positional mapping).

// Solver-internal board-cell symbols carried in the map (distinct from the 0..35 plaintext
// alphabet): 0..25 letters, SC_MAP_FIG the figure-shift marker, SC_MAP_NULL the unused cell.
#define SC_MAP_CELLS  28
#define SC_MAP_FIG    26
#define SC_MAP_NULL   27

#define STRADDLING_FILLER       23           // letter emitted for an invalid token ('X')
#define STRADDLING_VALID_WEIGHT 3.0          // token-validity reward weight (calibrated)
#define STRADDLING_N_CONFIGS    45           // C(10,2) indicator pairs
#define STRADDLING_KEEP         12           // top-K configs annealed (override with -nprimers)
#define STRADDLING_PREPASS_ITERS 36000       // mini-solve hill-climb length per config (ranking)

void solve_straddling_checkerboard(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

#endif
