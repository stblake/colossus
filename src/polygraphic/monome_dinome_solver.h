#ifndef MONOME_DINOME_SOLVER_H
#define MONOME_DINOME_SOLVER_H
#include "colossus.h"
#include "monome_dinome.h"

// =====================================================================
//  Monome-Dinome solver (TYPE monome-dinome / md)
// =====================================================================
//
// A CipherModel anneal (SHAPE_ANNEAL), the letter-only cousin of the Straddling
// Checkerboard: variable-length digit-stream decode -> cyclic-tile to a fixed length ->
// a token-validity reward folded into score_adjust. Two structural pieces:
//
//   * A 45-config STRUCTURAL SWEEP over the row-indicator pair {r0,r1}, pre-filtered by
//     TOKEN VALIDITY. The digit-stream tokenization depends ONLY on which two of the ten
//     digits are the row indicators (a digit equal to an indicator starts a 2-digit
//     "dinome"; any other digit is a 1-digit "monome"), so there are exactly C(10,2)=45
//     token structures. Crucially, validity is MAP-INDEPENDENT and purely structural: a
//     dinome whose second digit is itself an indicator is impossible on a real board, so
//     the TRUE indicator pair is always EXACTLY 100% valid while a wrong pair drops below
//     as soon as a dinome's second digit lands on an indicator. So the pre-pass ranks all
//     45 by structural validity and keeps only the fully-valid set (1-6 configs typically,
//     ALWAYS including the true pair -- the "decoupling reward" pattern), padded to
//     MD_MIN_CANDIDATES for a slightly indel-corrupted ciphertext. Each kept config is then
//     mini-solved (a short substitution anneal) to warm-start its map. This is far cheaper
//     and safer than an n-gram ranking of all 45 (which n-gram GAMING can defeat at short
//     lengths -- a free substitution over a mis-segmented stream out-scores the constrained
//     truth on quadgrams).
//
//   * The annealed state is the FREE code->letter substitution `map[24]`: a bijection from
//     the 24 token codes (8 monome + 8+8 dinome) to the 24 board letters. Keyed column
//     labels fold into this map (they are not separately identifiable). Unlike the
//     Straddling Checkerboard there is NO figure-shift and NO null cell -- the box is a
//     clean 24-letter fill -- so the map is a plain permutation of 0..23 (no sentinels),
//     the smallest, cleanest member of the digit-stream family. Cell-swap anneal.
//
//   * CROSS-CONFIG SELECTION BY DICTIONARY COVERAGE. n-gram score alone does NOT identify the
//     true config: a free substitution over a mis-segmented (wrong-indicator) stream reaches
//     a HIGHER n-gram score than the true segmentation's true map (measured: a gamed config
//     out-scores real English by ~0.2 quintgram at 160 letters). So run_solver's score-based
//     global-best picks a gamed config. The fix: DRIVE THE ENGINE ONCE PER kept config and
//     select the winner by DICTIONARY COVERAGE of its best decode -- the fraction of chars
//     covered by dictionary words (greedy longest match, min length 4). Real plaintext covers
//     ~0.6-0.8; a gamed mis-segmentation ~0.2-0.3, a gap n-gram gaming cannot fake. Coverage
//     is evaluated only ~1-6 times (once per kept config), so it is cheap. With no dictionary
//     loaded the selector falls back to the n-gram score (and gaming can then mislead).
//
// The ciphertext DIGIT stream is parsed from the raw string (not an A-Z decode). The
// plaintext alphabet is the 24-letter Monome-Dinome set (A..Z with J->I, Z->Y), forced by
// -type md. Effectively NEEDS -logprob with QUINTGRAMS (-ngramsize 5): the free
// substitution over a mis-segmented stream games a quadgram score, so quintgrams (which
// resist gaming) are what let the true config's true map win. Recovers reliably from ~120
// letters; the ACA low end (60-100) is below the blind floor (a documented limitation).
// Cribs are unused (the length change breaks the positional mapping).

#define MD_MAP_CELLS      24       // the code->letter map is a permutation of 0..23
#define MD_FILLER         22       // letter emitted for an invalid token ('X' in the 24-set)
#define MD_VALID_WEIGHT   3.0      // token-validity reward weight (calibrated)
#define MD_N_CONFIGS      45       // C(10,2) indicator pairs
#define MD_KEEP           12       // cap on kept configs (override with -nprimers)
#define MD_MIN_CANDIDATES 4        // keep at least this many by validity (indel-corruption insurance)
#define MD_VALID_EPS      1e-6     // "fully valid" tolerance around the max validity
#define MD_PREPASS_ITERS  36000    // mini-solve anneal length per kept config (warm-start)

void solve_monome_dinome(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

#endif
