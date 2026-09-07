#ifndef TRIDIGITAL_SOLVER_H
#define TRIDIGITAL_SOLVER_H
#include "colossus.h"
#include "tridigital.h"

// =====================================================================
//  Tridigital solver (TYPE tridigital / td)
// =====================================================================
//
// A CipherModel anneal (SHAPE_ANNEAL) over the digit-stream Tridigital cipher. Tridigital is
// the ONLY member of the digit-stream family whose decode is genuinely AMBIGUOUS: a letter's
// row does not affect its digit, so each of the 10 digits stands for up to 3 plaintext letters
// (9 "letter" digits carry 26 letters as 8 columns of 3 + 1 of 2; the 10th digit is the word
// separator). There is therefore no unique inverse -- recovering plaintext is a language-model
// problem. Three structural pieces:
//
//   * The SEARCHED KEY is a partition of the 26 letters into 9 groups (`key[L]` = group 0..8),
//     one group per non-separator digit. Because 26 letters in 9 groups capped at 3 forces the
//     sizes to 8x3 + 1x2, the perturbations are size-preserving (swap two letters' groups, or
//     relocate a letter into the current 2-group so the 2-group can wander). The separator
//     digit is a per-config parameter, not part of the key.
//
//   * The DECODE is an inner DP / beam Viterbi. Given the key + separator, each ciphertext
//     digit yields <=3 candidate letters (the separator yields a word break). We choose one
//     letter per position to maximize the n-gram score of the resulting letter stream. Because
//     ngram_score treats spaces as TRANSPARENT (letters-only projection), the DP carries the
//     previous n-1 letters ACROSS word boundaries. A small beam (TRI_BEAM) keeps it cheap
//     during the anneal; the final winner is re-decoded with a wider beam (TRI_FINAL_BEAM).
//
//   * CONFIG DISCRIMINATION over the 10 candidate separator digits, then CROSS-CONFIG SELECTION
//     BY DICTIONARY COVERAGE (not raw n-gram). A map-independent pre-filter ranks the separator
//     digits by how well the induced word-length histogram (gaps between separators) fits
//     english_word_length_frequencies[]; the top TRI_KEEP are mini-solved to warm-start a key.
//     The engine is then driven once per kept config and the winner chosen by dictionary
//     coverage of its decode -- a signal n-gram GAMING cannot fake. Tridigital is a DENSE
//     polyphonic cipher (9 symbols for 26 letters) with a free max-over-disambiguations decode,
//     so gaming is severe; word boundaries are known (from the separator), which makes the
//     dictionary signal strong. Effectively NEEDS -logprob with QUINTGRAMS and a -dictionary.
//     The ACA low end (75-100 letters) is below the blind floor (a documented limitation).

#define TRI_MAP_CELLS      26     // the key is a group id (0..8) per letter
#define TRI_NGROUPS        9      // 9 letter-groups (one per non-separator digit)
#define TRI_GROUP_CAP      3      // max letters per group (forces sizes 8x3 + 1x2 over 26)
#define TRI_N_CONFIGS      10     // 10 candidate separator digits
#define TRI_KEEP           4      // separator configs kept after the word-length pre-filter (-nprimers)
#define TRI_MIN_CANDIDATES 3      // keep at least this many (insurance for a noisy word-length fit)
#define TRI_PREPASS_ITERS  18000  // mini-solve anneal length per kept config (warm-start)
#define TRI_WORD_WEIGHT    3.0    // exact-dictionary-word reward in score_adjust (fights n-gram gaming)
#define TRI_MIN_WORD       4      // min word length for the wordset (shorter words hit by chance)
#define TRI_BEAM           4      // decode beam width during the search (DP; balances cost/gradient)
#define TRI_FINAL_BEAM     32     // decode beam width for the reported/stashed winner
#define TRI_BEAM_MAX       64     // hard cap on beam scratch (>= both beams)

void solve_tridigital(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);

#endif
