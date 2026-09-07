#ifndef LAYERED_SOLVER_H
#define LAYERED_SOLVER_H
#include "colossus.h"

// =====================================================================
//  Layered cipher: outer Quagmire III  o  inner columnar transposition
//  (QUAG_TRANS, aliases quagtrans/qtrans) -- the Paradigm PK3/PK4/PK6 forms.
// =====================================================================
//
//  Cryptanalytic model (confirmed empirically on the Paradigm challenge,
//  see memory/paradigm-pk-series.md and the module .c header):
//
//    CT = Quag_P( Trans( PT ) )
//
//  where the OUTERMOST layer is a single Quagmire III over a keyed alphabet
//  (KRYPTOS) with period P, and the layer underneath is a columnar
//  transposition of the plaintext (depth 0 = none / PK3, depth 1 = single
//  columnar / PK4, depth 2 = double columnar / PK6). A "double Quagmire III"
//  over one alphabet collapses to a single Quagmire of period lcm (Quagmire
//  III is a Vigenere in keyed-index space, so stacked layers sum), so the
//  outer stage is always one Quagmire regardless of how many Q layers the
//  puzzle notation lists -- pass P = lcm of their periods via -cyclewordlen.
//
//  The attack exploits that a transposition PRESERVES MONOGRAMS: the outer
//  cycleword is recovered by the monogram statistic (derive_optimal_cycleword)
//  even though the intermediate is scrambled English; the inner transposition
//  is then solved by n-gram; finally the cycleword is polished by n-gram
//  coordinate-ascent THROUGH the recovered transposition (the monogram warm
//  start under-resolves the ~P thin columns). depth 2 reuses the divide-and-
//  conquer double-transposition core (dct_solve_core).
//
//  Params: -cyclewordlen P (outer Quag period; else estimated by periodic IoC),
//  -plaintextkeyword KRYPTOS (keyed alphabet; required), -depth {0,1,2} (inner
//  columnar stage count), -mincols/-maxcols (inner column-count sweep),
//  -readdir tb|bt|both. -logprob is effectively required. Cribs (-crib) are on
//  final-plaintext positions and steer the inner solve + cycleword refine.

void solve_quag_trans(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs);

// Strip a plain (no-transposition) period-P Quagmire III of English from `cipher`:
// derive the cycleword by the monogram statistic, then polish it by n-gram
// coordinate ascent. Fills cw_out[0..P-1] (letter indices) and pt_out[0..len-1],
// and returns the n-gram fitness of the recovered plaintext. Shared with the
// Hill-Quag solver, which calls it on the Hill-stripped intermediate.
double quag_strip_and_refine(ColossusConfig *cfg, const int *cipher, int len,
    int pt_kw[], int ct_kw[], int P, float *ngram,
    int *crib_indices, int *crib_positions, int n_cribs,
    int cw_out[], int pt_out[]);

#endif
