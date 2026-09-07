#ifndef KEYPHRASE_SOLVER_H
#define KEYPHRASE_SOLVER_H
#include "colossus.h"
#include "keyphrase.h"

// =====================================================================
//  Key Phrase solver (TYPE keyphrase / kp, ACA)
// =====================================================================
//
// Cryptanalytically the TWIN of Tridigital: an ambiguous many-to-one decode resolved by an
// inner beam-Viterbi, but SINGLE-CONFIG. Because the ACA Key Phrase retains word divisions in
// the ciphertext, there is no separator to discover (unlike Tridigital's separator-digit
// search), so one CipherModel anneal solves it directly.
//
//   * SEARCHED KEY -- a partition of the 26 plaintext letters among the K DISTINCT ciphertext
//     letters that occur: grp[p] in {0..K-1} = which observed-ct-letter group plaintext letter
//     p maps to. (The phrase induces this partition; the group labels are the observed ct
//     letters, read straight off the ciphertext.) Moves reassign / swap letters between groups,
//     keeping every group non-empty (each occurring ct letter needs >= 1 plaintext letter).
//
//   * DECODE -- an inner beam Viterbi over the letters-only ciphertext: at each position the ct
//     letter's group yields its candidate plaintext letters, and we choose per position to
//     maximize the n-gram score. ngram_score is letters-only-transparent, so the search runs on
//     the spaces-stripped stream (word divisions are restored only in the report).
//
// Blind recovery is markedly better than Tridigital (K ~ 15-19 groups for a 75-100 letter cipher
// vs Tridigital's 26-in-9), but residual errors persist where a ct letter covers several
// plaintext letters with weak context. -logprob (+ quintgrams) recommended.

#define KP_BEAM        8    // decode beam width during the anneal (cost / gradient balance)
#define KP_FINAL_BEAM  32   // decode beam width for the reported winner
#define KP_BEAM_MAX    64   // hard cap on beam scratch (>= both beams)

// Exposed for the solver test's move-invariant check: one non-empty-preserving move on grp[].
void keyphrase_move(int grp[KEYPHRASE_LEN], int K);

void solve_keyphrase(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result);
#endif
