#ifndef PERIOD_COLUMN_SPACE_SOLVER_H
#define PERIOD_COLUMN_SPACE_SOLVER_H
#include "colossus.h"
#include "period_column_solver.h"   // PCStage, period_column_search

// =====================================================================
//  Period column order (space-robust) transposition solver
//  (PERIOD_COLUMN_SPACE)
// =====================================================================
//
//  A generalisation of the deterministic period-column solver for POORLY
//  ENCIPHERED ciphertext -- one where a few characters were accidentally
//  dropped (or, symmetrically, where the grid alignment is a little off). The
//  observed ciphertext is kept EXACTLY, spaces and all: the spaces are
//  structurally significant (they ride the transposition as real grid cells and
//  a single misplaced space re-shuffles every downstream letter), so they are
//  never stripped or moved. Instead the solver INSERTS additional blank "gap"
//  cells (synthesised spaces) at SEARCHED positions:
//
//    * a gap models a dropped character -- it restores the grid coordinates of
//      every following cell, re-aligning the columns; and
//    * the extra cells let the padded length N+g reach a richer set of
//      complete-grid factorisations than the observed length N alone.
//
//  The corruption may go BOTH ways, so the search is symmetric:
//    * INSERT up to max_ins blank gap cells (repair DROPPED characters -- the
//      observed cipher is too short); and
//    * DELETE up to max_dels observed cells (repair spuriously ADDED characters
//      -- the observed cipher is too long).
//  A gap models a lost cell; a deletion removes an extraneous one. Both fix the
//  downstream grid alignment and let the edited length reach a complete-grid
//  factorisation. For each (n_ins, n_dels) pair in [0,max_ins] x [0,max_dels]
//  (0,0 is exactly the plain period-column search, so this solver never scores
//  below it) the edit POSITIONS are searched FROM SCRATCH -- random seeds over
//  the whole cipher, not nudges of the observed spacing -- by a restart+anneal
//  hill climb whose fitness is the fast exhaustive period-column search
//  (inner_depth, typically 1) on the edited stream. The single best edited
//  cipher then gets a full final_depth (<= 2) exhaustive finish. Honours
//  -nrestarts / -nhillclimbs / -inittemp / -mintemp like the other stochastic
//  solvers; -maxgaps sets the insertion budget, -maxdels the deletion budget,
//  -depth the final composition depth.
//
//  Cribs are NOT supported here (inserting gaps shifts plaintext positions, so a
//  positional crib no longer lines up); use the plain period-column solver when
//  the cipher is known-clean and cribs are available.

// Run the edit search. cipher[0..len-1] is the observed ciphertext (letters as
// 0..25, spaces/punctuation as negative sentinels). Writes the best edited cipher
// (length *out_len = len + inserted gaps - deleted cells) to out_edited, its
// decrypt to out_decrypt, the winning period-column stages to stages[0..*n_stages-1],
// the recovered gap insertion slots (0..len) to gap_positions[0..*n_gaps-1], and
// the recovered deleted observed indices (0..len-1, ascending) to
// del_positions[0..*n_dels-1]. Returns the best n-gram(+word) score. out_edited /
// out_decrypt must hold at least len + max_ins ints; gap_positions at least max_ins,
// del_positions at least max_dels.
double period_column_space_search(const int *cipher, int len,
    int max_ins, int max_dels, int inner_depth, int final_depth,
    int n_restarts, int n_hillclimbs, double init_temp, double min_temp,
    float *ngram_data, int ngram_size,
    float weight_ngram, float weight_word,
    SharedData *shared, bool use_words,
    int out_edited[], int *out_len,
    int out_decrypt[], PCStage stages[], int *n_stages,
    int gap_positions[], int *n_gaps,
    int del_positions[], int *n_dels, bool verbose);

// Full solver entry, dispatched from solve_cipher().
void solve_period_column_space(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs);
#endif
