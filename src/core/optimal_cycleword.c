//
// Optimal cycleword derivation
//

#include "colossus.h"

/*
   derive_optimal_cycleword
   ========================

   Determines the statistically most likely cycleword (key) for a given set of
   plaintext and ciphertext alphabets. Instead of perturbing the cycleword
   stochastically, this routine deterministically solves for the optimal key
   character for each column of the period.

   ## Mathematical Model

   Let L be the period (cycleword length) and C the ciphertext of length N.
   Partition C into L columns, the k-th column being all C_i with i = k (mod L).
   For each column we pick the key character that maximizes the correlation
   between the decrypted column's letter frequencies and English monograms E.

   For a candidate shift s, the column score is the dot product of the decrypted
   frequencies with E:

       S_s = sum_i f^(s)_i * E_i

   ## Implementation (column histogram * shift-weight correlation)

   The decryption D(c, s) of a ciphertext char c under shift s does not depend on
   which column c sits in, so the per-shift, per-char monogram weight

       weight[s][c] = E[ D(c, s) ]

   is built ONCE per call (a 26x26 table). Each column is then scored by
   histogramming its raw ciphertext characters once and forming

       S_s = sum_c hist[c] * weight[s][c].

   This replaces the previous "re-decrypt the whole column for each of the 26
   shifts" inner loop: per-column work drops from O(26 * column_length) to
   O(column_length + 26*26). derive_optimal_cycleword is the dominant cost of the
   default -optimalcycle hill climb, so this is the hot path.

   ## Decryption functions D(c, s)

   * Vigenere : D = (c - s) mod 26          (variant: (s - c) mod 26)
   * Beaufort : D = (s - c) mod 26
   * Porta    : S = floor(s/2); c in [0,12] -> (c + S) mod 13 + 13,
                                 c in [13,25] -> (c - 13 - S) mod 13
   * Quagmire : p = position of c in the CT keyed alphabet;
                D = A_pt[(p - s) mod 26]     (variant: A_pt[(p + s) mod 26])
*/

void derive_optimal_cycleword(
    ColossusConfig *cfg,
    int cipher_indices[], int cipher_len,
    int plaintext_keyword_indices[], int ciphertext_keyword_indices[],
    int cycleword_state[], int cycleword_len, int *hist_by_col) {

    int i, c, s, col, row;

    // Position of each ciphertext character within the CT keyed alphabet
    // (ciphertext_keyword_indices is always a permutation of 0..25).
    int ct_key_lookup[ALPHABET_SIZE];
    for (i = 0; i < g_alpha; i++) ct_key_lookup[ciphertext_keyword_indices[i]] = i;

    // Quagmire I-IV are handled by a factored fast path (below); the cached 26x26
    // weight table is only built for the keyword-free ciphers. Gronsfeld has straight
    // alphabets like Vigenere, so it takes the keyword-free path too.
    int is_quag = !(cfg->cipher_type == PORTA || cfg->cipher_type == BEAUFORT ||
                    cfg->cipher_type == VIGENERE || cfg->cipher_type == GRONSFELD);

    // Gronsfeld shifts are key digits 0..9, so its per-column shift search is bounded
    // to that domain (vs the full 0..25); every other type sweeps the whole alphabet.
    int smax = (cfg->cipher_type == GRONSFELD) ? GRONSFELD_DIGITS : g_alpha;

    // Quagmire factoring.
    // ===================
    // For Quagmire the per-shift, per-char monogram weight is
    //
    //     weight[s][c] = monogram[ pt_kw[ (pos[c] -/+ s) mod 26 ] ]
    //
    // where pos[c] = ct_key_lookup[c] is c's position in the CT keyed alphabet and
    // -/+ selects standard/variant. The right-hand monogram depends on the PT
    // keyword ONLY through the 26-vector
    //
    //     M[q] = monogram[ pt_kw[q] ],
    //
    // so weight[s][c] = M[(pos[c] -/+ s) mod 26]. Building M (26 entries) per call
    // and indexing it directly replaces the previous 676-entry weight-table
    // rebuild, which never cached for Quagmire because the keyword changes on
    // essentially every hill-climb iteration -- the dominant cost of the default
    // -optimalcycle climb. To drop the per-access `mod 26`, M is laid out twice
    // back-to-back in Mext[0..51] (Mext[i] = M[i % 26]); then for nonzero column
    // entries at position p = pos[c],
    //
    //     standard: weight = Mext[p + 26 - s]   (p+26-s in 1..51)
    //     variant : weight = Mext[p + s]        (p+s    in 0..50).
    //
    // The looked-up double is bit-for-bit the old weight[s][c] (same monogram
    // array element), and the per-shift dot product keeps the same operands and
    // summation order, so every column's argmax -- and the derived cycleword -- is
    // unchanged.
    double Mext[2 * ALPHABET_SIZE];
    if (is_quag) {
        for (i = 0; i < g_alpha; i++) {
            double m = g_monograms[plaintext_keyword_indices[i]];
            Mext[i] = m;
            Mext[i + g_alpha] = m;
        }
    }

    // The 26x26 monogram-weight table for the keyword-free ciphers (Vigenere/
    // Beaufort/Porta) depends only on (cipher_type, variant); their keywords never
    // change, so it is built once for the whole run and cached. (We still compare
    // the cached keywords so a stray batch run with a different type rebuilds.)
    static _Thread_local double weight[ALPHABET_SIZE][ALPHABET_SIZE];
    static _Thread_local int w_valid = 0, w_type = -1, w_variant = -1;
    static _Thread_local int w_pt[ALPHABET_SIZE], w_ct[ALPHABET_SIZE];

    if (!is_quag) {
        int rebuild = !w_valid || w_type != cfg->cipher_type || w_variant != (int)cfg->variant;
        if (!rebuild) {
            for (i = 0; i < g_alpha; i++) {
                if (w_pt[i] != plaintext_keyword_indices[i] ||
                    w_ct[i] != ciphertext_keyword_indices[i]) { rebuild = 1; break; }
            }
        }
        if (rebuild) {
            for (s = 0; s < g_alpha; s++) {
                for (c = 0; c < g_alpha; c++) {
                    int pt_char;
                    if (cfg->cipher_type == PORTA) {
                        int porta_shift = s / 2;
                        if (c < 13) pt_char = (c + porta_shift) % 13 + 13;
                        else        pt_char = (c - 13 - porta_shift + ALPHABET_SIZE) % 13;
                    } else if (cfg->cipher_type == BEAUFORT) {
                        pt_char = (s - c + g_alpha) % g_alpha;
                    } else {
                        // Vigenere. Must match vigenere_decrypt: standard
                        // P = (C - K), variant P = (C + K) mod 26.
                        if (cfg->variant) pt_char = (c + s) % g_alpha;
                        else              pt_char = (c - s + g_alpha) % g_alpha;
                    }
                    weight[s][c] = g_monograms[pt_char];
                }
            }
            w_valid = 1; w_type = cfg->cipher_type; w_variant = (int)cfg->variant;
            for (i = 0; i < g_alpha; i++) {
                w_pt[i] = plaintext_keyword_indices[i];
                w_ct[i] = ciphertext_keyword_indices[i];
            }
        }
    }

    // Solve each column independently for the shift best matching English. The
    // dropped per-column normalization (1/column_length) is constant across
    // shifts, so it does not affect the argmax; ties go to the lowest shift, as
    // before.
    int local_hist[ALPHABET_SIZE];
    for (col = 0; col < cycleword_len; col++) {
        int *hist;
        if (hist_by_col) {
            hist = hist_by_col + col * ALPHABET_SIZE;
        } else {
            // Standalone path: histogram this column's ciphertext chars locally.
            for (c = 0; c < g_alpha; c++) local_hist[c] = 0;
            for (row = 0; row * cycleword_len + col < cipher_len; row++)
                local_hist[cipher_indices[row * cycleword_len + col]]++;
            hist = local_hist;
        }

        // Compact the column's nonzero histogram entries once (in increasing c
        // order). Columns hold only ~cipher_len/cycleword_len characters, so most
        // of the 26 counts are zero. Skipping the zero terms is exact -- adding
        // hist[c]*weight == 0.0 is a no-op in IEEE arithmetic -- so the per-shift
        // scores, and hence the argmax, are bit-for-bit unchanged. This turns the
        // inner dot product from 26 multiply-adds into (#distinct chars) of them,
        // and derive_optimal_cycleword is the dominant cost of the optimal climb.
        int nz_c[ALPHABET_SIZE], nz_n[ALPHABET_SIZE], n_nz = 0;
        for (c = 0; c < g_alpha; c++) {
            if (hist[c]) { nz_c[n_nz] = c; nz_n[n_nz] = hist[c]; n_nz++; }
        }

        double best_score = -1.0;
        int best_shift = 0;
        if (is_quag) {
            // Index Mext directly (see the factoring note above). pp[k] folds in
            // the per-entry CT-keyed-alphabet position so the shift loop is a flat
            // gather; the resulting dot products equal the old weight[s][.] ones.
            int pp[ALPHABET_SIZE];
            if (cfg->variant) {
                for (int k = 0; k < n_nz; k++) pp[k] = ct_key_lookup[nz_c[k]];
                for (s = 0; s < g_alpha; s++) {
                    double score = 0.0;
                    for (int k = 0; k < n_nz; k++) score += nz_n[k] * Mext[pp[k] + s];
                    if (score > best_score) { best_score = score; best_shift = s; }
                }
            } else {
                for (int k = 0; k < n_nz; k++) pp[k] = ct_key_lookup[nz_c[k]] + g_alpha;
                for (s = 0; s < g_alpha; s++) {
                    double score = 0.0;
                    for (int k = 0; k < n_nz; k++) score += nz_n[k] * Mext[pp[k] - s];
                    if (score > best_score) { best_score = score; best_shift = s; }
                }
            }
        } else {
            for (s = 0; s < smax; s++) {
                double score = 0.0;
                const double *ws = weight[s];
                for (int k = 0; k < n_nz; k++) score += nz_n[k] * ws[nz_c[k]];
                if (score > best_score) { best_score = score; best_shift = s; }
            }
        }

        // State stores the cycleword CHARACTER (from the CT keyed alphabet).
        cycleword_state[col] = ciphertext_keyword_indices[best_shift];
    }
}
