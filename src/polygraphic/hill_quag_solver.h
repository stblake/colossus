#ifndef HILL_QUAG_SOLVER_H
#define HILL_QUAG_SOLVER_H
#include "colossus.h"

// =====================================================================
//  Layered cipher: outer Hill(k x k)  o  inner Quagmire III  (HILL_QUAG,
//  aliases hillquag/hq) -- the Paradigm PK7 form  Q(6)H(3x3).
// =====================================================================
//
//    CT = Hill_k( Quag_P( PT ) )
//
//  The outer Hill flattens all periodicity, so unlike the quagtrans forms the
//  raw ciphertext shows no period-P peak. The attack exploits a clean
//  DECOUPLING: with the Hill DECRYPTION matrix D, D.CT = Quag_P(PT), whose
//  columns (positions == c mod P) are each a keyed-shift of an English column,
//  hence peaked -- and this "columns look like English after their best shift"
//  statistic is INDEPENDENT of the Quagmire cycleword.
//
//  When k divides P (e.g. k=3, P=6), each mod-P column of D.CT is produced by
//  exactly ONE row of D, so the rows are recovered INDEPENDENTLY: for each row,
//  exhaustively score all 26^k candidate rows by the monogram-correlation of the
//  columns it generates, keep the top few, then combine (final n-gram breaks the
//  multiplicative-unit ambiguity). Otherwise a full-matrix anneal over the mean
//  period-P columnar IoC is used. With D fixed, strip the Hill and solve the
//  plain period-P Quagmire (quag_strip_and_refine). -logprob is effectively
//  required; blind Hill recovery is length-limited (see the module .c header).
//
//  Params: -period k (Hill block, default 3), -cyclewordlen P (inner Quag
//  period; required -- Hill hides it), -plaintextkeyword KRYPTOS. Cribs steer
//  the final Quagmire solve (positions are preserved PT<->intermediate).

void solve_hill_quag(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs);

#endif
