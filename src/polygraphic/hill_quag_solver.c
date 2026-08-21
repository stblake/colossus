#include "hill_quag_solver.h"
#include "layered_solver.h"
#include "scoring.h"
#include <string.h>
#include <math.h>

// =====================================================================
//  HILL_QUAG: outer Hill(k) o inner Quagmire III. See the header for the model.
//
//  Blind Hill recovery is intrinsically length-limited (the standard blind
//  Hill-3x3 floor is ~1200 letters), but the per-row monogram-correlation
//  decoupling below is far denser than an n-gram climb: each row's ~2 columns
//  carry ~2N/P samples, so the true row's columns visibly out-correlate English
//  monograms. When k | P the rows are recovered by an EXACT 26^k-per-row sweep;
//  otherwise a full-matrix anneal over the period-P columnar IoC is the (weaker,
//  experimental) fallback. Reported as best-effort for the shortest ACA lengths.
// =====================================================================

#define TOPT    40                    // top candidate rows kept per row position (anchors)

// Position of each ciphertext letter within the CT keyed alphabet.
static void build_ct_lookup(const int ct_kw[], int lookup[]) {
    for (int i = 0; i < g_alpha; i++) lookup[ct_kw[i]] = i;
}

// Best monogram correlation of a column (values in vals[] at start, start+stride, ...)
// over the 26 Quagmire de-shifts: max_s sum_c hist[c]*monogram[pt_kw[(pos_ct(c)-s) mod 26]].
// Mp[q] = monogram[pt_kw[q]] is precomputed once by the caller.
static double best_col_corr(const int *vals, int start, int stride, int total,
                            const double *Mp, const int *ct_lookup) {
    int hist[ALPHABET_SIZE]; memset(hist, 0, sizeof(hist));
    for (int i = start; i < total; i += stride) hist[vals[i]]++;
    int nz_c[ALPHABET_SIZE], nz_n[ALPHABET_SIZE], nnz = 0;
    double norm2 = 0.0;
    for (int c = 0; c < g_alpha; c++) if (hist[c]) {
        nz_c[nnz] = ct_lookup[c]; nz_n[nnz] = hist[c]; nnz++;
        norm2 += (double) hist[c] * hist[c];
    }
    double best = -1e18;
    for (int s = 0; s < g_alpha; s++) {
        double sc = 0.0;
        for (int j = 0; j < nnz; j++) {
            int q = nz_c[j] - s; if (q < 0) q += g_alpha;
            sc += nz_n[j] * Mp[q];
        }
        if (sc > best) best = sc;
    }
    // Normalize by the column histogram's L2 norm (cosine-style). Without this a row
    // that collapses a column onto one high-frequency letter games the raw dot product;
    // the norm makes an English-SHAPED distribution win over a spike.
    return (norm2 > 0.0) ? best / sqrt(norm2) : best;
}

// ---- exact per-row recovery (k | P) --------------------------------------
// For row position r, evaluate every 26^k candidate row: apply it to each block's
// within-block position r, bucket outputs into the P/k columns {r, r+k, ...} they
// feed, and score by the summed best_col_corr. Keep the TOPT best rows in cand[].
typedef struct { int row[HILL_MAX_K]; double score; } RowCand;

static void recover_row(const int *cipher, int len, int k, int P, int r,
                        const int *pt_kw, const int *ct_lookup, const int *to_letter,
                        RowCand cand[TOPT]) {
    int nblocks = len / k;
    int ncols = P / k;                        // number of mod-P columns this row feeds
    for (int t = 0; t < TOPT; t++) cand[t].score = -1e18;

    double Mp[ALPHABET_SIZE];
    for (int q = 0; q < g_alpha; q++) Mp[q] = g_monograms[pt_kw[q]];

    int row[HILL_MAX_K]; memset(row, 0, sizeof(row));
    int outs[MAX_CIPHER_LENGTH];              // block outputs (mapped to letters) for this row
    long total = 1; for (int j = 0; j < k; j++) total *= g_alpha;   // 26^k

    for (long idx = 0; idx < total; idx++) {
        for (int b = 0; b < nblocks; b++) {
            const int *v = cipher + b * k;
            int acc = 0;
            for (int j = 0; j < k; j++) acc += row[j] * v[j];
            outs[b] = to_letter[((acc % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE];
        }
        double sc = 0.0;
        for (int c = 0; c < ncols; c++)       // column c = blocks b == c (mod ncols)
            sc += best_col_corr(outs, c, ncols, nblocks, Mp, ct_lookup);
        if (sc > cand[TOPT - 1].score) {      // insert into the sorted TOPT list
            int p = TOPT - 1;
            while (p > 0 && cand[p - 1].score < sc) { cand[p] = cand[p - 1]; p--; }
            memcpy(cand[p].row, row, k * sizeof(int));
            cand[p].score = sc;
        }
        for (int j = k - 1; j >= 0; j--) { if (++row[j] < g_alpha) break; row[j] = 0; }
    }
}

// Quick fitness of a Hill decryption matrix D (in the working space): apply it, map
// back to letters, monogram-strip the inner Quagmire (no refine), and return the
// n-gram score. Singular matrices score -inf. `scr` is caller scratch (>= len ints).
static double score_D(ColossusConfig *cfg, const int *D, int k, const int *wct, int len,
                      const int *to_letter, int P, int pt_kw[], int ct_kw[],
                      float *ngram, int *scr) {
    if (hill_mod_inverse(hill_det_mod(D, k), ALPHABET_SIZE) == 0) return -1e18;
    int Iw[MAX_CIPHER_LENGTH], I[MAX_CIPHER_LENGTH], scw[MAX_CYCLEWORD_LEN];
    hill_mat_mul_blocks(D, k, wct, len, Iw);
    for (int i = 0; i < len; i++) I[i] = to_letter[Iw[i]];
    derive_optimal_cycleword(cfg, I, len, pt_kw, ct_kw, scw, P, NULL);
    quagmire_decrypt(scr, I, len, pt_kw, ct_kw, scw, P, cfg->variant);
    return ngram_score(scr, len, ngram, cfg->ngram_size);
}

// Robust matrix finder for one convention. The per-row cosine ranking usually puts
// each true row in the top-TOPT, but at ~2N/P samples per column noise can push one
// row deep (rank >> TOPT). So for EACH choice of a "free" row position, anchor the
// other k-1 rows from their top-TOPT candidates (beam by the n-gram of the partial
// decrypt), then search the free row EXHAUSTIVELY over 26^k by the same n-gram --
// with two correct anchors the third row is pinned regardless of its cosine rank.
// Fills bestD (best over all free-position choices) and returns its score.
#define ANCHOR_B 10                   // anchored (k-1)-row combos kept per free position
static double find_matrix_conv(ColossusConfig *cfg, const int *wct, int len,
        const int *to_letter, int k, int P, int pt_kw[], int ct_kw[], const int *ct_lookup,
        float *ngram, int *bestD) {
    RowCand cand[HILL_MAX_K][TOPT];
    for (int r = 0; r < k; r++)
        recover_row(wct, len, k, P, r, pt_kw, ct_lookup, to_letter, cand[r]);

    int scr[MAX_CIPHER_LENGTH];
    double best = -1e18;

    for (int f = 0; f < k; f++) {                    // f = the exhaustively-searched row
        // Anchor the other k-1 rows: beam over their top-TOPT candidates by the n-gram
        // of the full decrypt (free row f held at its best-cosine guess).
        double topb_s[ANCHOR_B]; int topb_sel[ANCHOR_B][HILL_MAX_K]; int nb = 0;
        for (int i = 0; i < ANCHOR_B; i++) topb_s[i] = -1e18;

        int sel[HILL_MAX_K]; memset(sel, 0, sizeof(sel));   // index into cand[pos] for each pos!=f
        long anchors = 1; for (int r = 0; r < k; r++) if (r != f) anchors *= TOPT;
        for (long a = 0; a < anchors; a++) {
            int D[HILL_MAX_KEY];
            for (int r = 0; r < k; r++)
                if (r != f) memcpy(D + r * k, cand[r][sel[r]].row, k * sizeof(int));
            // Score the anchored (k-1) rows with a few free-row guesses; a good pair may
            // be singular with one guess but not another, so keep the best.
            double s = -1e18;
            for (int g = 0; g < 3 && g < TOPT; g++) {
                memcpy(D + f * k, cand[f][g].row, k * sizeof(int));
                double sg = score_D(cfg, D, k, wct, len, to_letter, P, pt_kw, ct_kw, ngram, scr);
                if (sg > s) s = sg;
            }
            if (s > topb_s[ANCHOR_B - 1]) {
                int p = ANCHOR_B - 1;
                while (p > 0 && topb_s[p-1] < s) { topb_s[p]=topb_s[p-1]; memcpy(topb_sel[p],topb_sel[p-1],sizeof(topb_sel[0])); p--; }
                topb_s[p] = s; memcpy(topb_sel[p], sel, sizeof(sel));
                if (nb < ANCHOR_B) nb++;
            }
            for (int r = k - 1; r >= 0; r--) { if (r == f) continue; if (++sel[r] < TOPT) break; sel[r] = 0; }
        }

        // For each strong anchored (k-1)-row set, search the free row over 26^k.
        for (int b = 0; b < nb; b++) {
            int D[HILL_MAX_KEY];
            for (int r = 0; r < k; r++)
                if (r != f) memcpy(D + r * k, cand[r][topb_sel[b][r]].row, k * sizeof(int));
            int row[HILL_MAX_K]; memset(row, 0, sizeof(row));
            long tot = 1; for (int j = 0; j < k; j++) tot *= g_alpha;
            for (long idx = 0; idx < tot; idx++) {
                memcpy(D + f * k, row, k * sizeof(int));
                double s = score_D(cfg, D, k, wct, len, to_letter, P, pt_kw, ct_kw, ngram, scr);
                if (s > best) { best = s; memcpy(bestD, D, k * k * sizeof(int)); }
                for (int j = k - 1; j >= 0; j--) { if (++row[j] < g_alpha) break; row[j] = 0; }
            }
        }
    }
    return best;
}

// ---- full-matrix anneal fallback (k does not divide P) --------------------
static double matrix_period_ioc(const int *D, int k, const int *cipher, int len, int P) {
    int I[MAX_CIPHER_LENGTH];
    hill_mat_mul_blocks(D, k, cipher, len, I);
    double sum = 0.0; int nc = 0;
    for (int c = 0; c < P; c++) {
        int hist[ALPHABET_SIZE]; memset(hist, 0, sizeof(hist)); int n = 0;
        for (int i = c; i < len; i += P) { hist[I[i]]++; n++; }
        if (n < 2) continue;
        double s = 0.0; for (int q = 0; q < g_alpha; q++) s += (double) hist[q] * (hist[q] - 1);
        sum += s / ((double) n * (n - 1)); nc++;
    }
    double ic = nc ? sum / nc : 0.0;
    if (hill_mod_inverse(hill_det_mod(D, k), ALPHABET_SIZE) == 0) ic -= 1.0;   // singular penalty
    return ic;
}

static void anneal_matrix(ColossusConfig *cfg, const int *cipher, int len, int k, int P, int D[]) {
    int restarts = cfg->n_restarts > 0 ? cfg->n_restarts : 400;
    int climbs   = cfg->n_hill_climbs > 0 ? cfg->n_hill_climbs : 20000;
    int n = k * k, best[HILL_MAX_KEY], cur[HILL_MAX_KEY];
    double best_score = -1e18;
    for (int r = 0; r < restarts; r++) {
        for (int i = 0; i < n; i++) cur[i] = rand_int(0, g_alpha);
        double cs = matrix_period_ioc(cur, k, cipher, len, P);
        for (int h = 0; h < climbs; h++) {
            int i = rand_int(0, n), old = cur[i];
            cur[i] = rand_int(0, g_alpha);
            double s = matrix_period_ioc(cur, k, cipher, len, P);
            if (s >= cs) cs = s; else cur[i] = old;
        }
        if (cs > best_score) { best_score = cs; memcpy(best, cur, n * sizeof(int)); }
    }
    memcpy(D, best, n * sizeof(int));
}

// ---- reporting ------------------------------------------------------------
static void report_hill_quag(ColossusConfig *cfg, SharedData *shared,
        const int *cipher, int len, char *cribtext_str, int k, int P,
        const int *D, int pt_kw[], const int *cw, const int *pt, double score, int n_cribs) {
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(pt[i]);
    plaintext_string[len] = '\0';
    int n_words = 0;
    if (cfg->dictionary_present && shared->dict != NULL)
        n_words = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);

    printf("\nhillquag: Hill(%dx%d) o Quagmire III (period %d)\n", k, k, P);
    printf("\nResult Score: %.2f | Words: %d\n", score, n_words);
    print_text((int *) cipher, len); printf("\n");
    print_text((int *) pt, len);     printf("\n");
    print_spaces_line(g_spaces_table, (int *) pt, len);
    if (n_cribs > 0) printf("%s\n", cribtext_str);

    // Recovered Hill ENCRYPTION matrix (invert the decryption matrix D).
    int E[HILL_MAX_KEY]; int ok = hill_mat_inverse(D, k, E);
    printf("hill decryption matrix:");
    for (int i = 0; i < k * k; i++) printf(" %d", D[i]);
    printf("\nhill encryption matrix:");
    if (ok) for (int i = 0; i < k * k; i++) printf(" %d", E[i]);
    else printf(" (singular)");
    printf("\nkeyed alphabet: "); print_text(pt_kw, g_alpha);
    printf("\ncycleword (len %d):", P);
    for (int i = 0; i < P; i++) printf(" %c", index_to_char(cw[i]));
    printf("\n\n");

    if (cfg->dictionary_present) printf(">>> %.2f, %d, %d, ", score, n_words, cfg->cipher_type);
    else                         printf(">>> %.2f, %d, ", score, cfg->cipher_type);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text((int *) cipher, len); printf(", ");
    print_text((int *) pt, len);     printf("\n");
    fflush(stdout);
}

// ---- entry ----------------------------------------------------------------
void solve_hill_quag(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: hillquag needs the 26-letter alphabet (Hill is mod 26).\n\n");
        return;
    }
    if (cipher_len < 16) { printf("\n\nERROR: ciphertext too short.\n\n"); return; }

    int pt_kw[ALPHABET_SIZE], ct_kw[ALPHABET_SIZE];
    if (cfg->user_plaintext_keyword_present)      make_keyed_alphabet(cfg->user_plaintext_keyword, pt_kw);
    else if (cfg->user_ciphertext_keyword_present) make_keyed_alphabet(cfg->user_ciphertext_keyword, pt_kw);
    else { printf("\n\nERROR: hillquag needs the keyed alphabet (-plaintextkeyword KRYPTOS).\n\n"); return; }
    memcpy(ct_kw, pt_kw, ALPHABET_SIZE * sizeof(int));

    int k = cfg->period_present ? cfg->period : 3;
    if (k < 2 || k > HILL_MAX_K) { printf("\n\nERROR: Hill block size k must be 2..%d.\n\n", HILL_MAX_K); return; }
    if (!cfg->cycleword_len_present) {
        printf("\n\nERROR: hillquag needs -cyclewordlen P (the inner Quagmire period; "
               "the Hill layer hides it).\n\n");
        return;
    }
    int P = cfg->cycleword_len;

    printf("hillquag: Hill k=%d, inner Quagmire period P=%d, keyed alphabet pinned\n", k, P);

    int ct_lookup[ALPHABET_SIZE];
    build_ct_lookup(ct_kw, ct_lookup);
    float *ngram = shared->ngram_data;

    int bestD[HILL_MAX_KEY];
    int cw[MAX_CYCLEWORD_LEN], pt[MAX_CIPHER_LENGTH], I[MAX_CIPHER_LENGTH];
    double best_score = -1e18;
    bool found = false;
    int best_conv = -1;

    if (k <= 4 && P % k == 0) {
        // EXACT-ish per-row recovery, tried in BOTH Hill spaces: the matrix may multiply
        // the ciphertext in the plain A..Z index (conv 0) or the KRYPTOS keyed index (conv
        // 1, "Alphabet: KRYPTOS"). find_matrix_conv recovers the Hill decryption matrix in
        // each space; keep the better full Quagmire solve.
        int identity[ALPHABET_SIZE];
        for (int i = 0; i < g_alpha; i++) identity[i] = i;
        int wctA[MAX_CIPHER_LENGTH], wctB[MAX_CIPHER_LENGTH];
        for (int i = 0; i < cipher_len; i++) { wctA[i] = cipher_indices[i]; wctB[i] = ct_lookup[cipher_indices[i]]; }

        for (int conv = 0; conv < 2; conv++) {
            int *wct       = (conv == 0) ? wctA : wctB;
            int *to_letter = (conv == 0) ? identity : ct_kw;   // matmul result -> letter
            int D[HILL_MAX_KEY];
            find_matrix_conv(cfg, wct, cipher_len, to_letter, k, P, pt_kw, ct_kw, ct_lookup, ngram, D);

            // Full Quagmire solve (strip + component refine) on the recovered matrix.
            int Iw[MAX_CIPHER_LENGTH];
            hill_mat_mul_blocks(D, k, wct, cipher_len, Iw);
            for (int j = 0; j < cipher_len; j++) I[j] = to_letter[Iw[j]];
            int tcw[MAX_CYCLEWORD_LEN], tpt[MAX_CIPHER_LENGTH];
            double s = quag_strip_and_refine(cfg, I, cipher_len, pt_kw, ct_kw, P, ngram,
                                             crib_indices, crib_positions, n_cribs, tcw, tpt);
            if (s > best_score) {
                best_score = s; found = true; best_conv = conv;
                memcpy(bestD, D, k * k * sizeof(int));
                memcpy(cw, tcw, P * sizeof(int));
                memcpy(pt, tpt, cipher_len * sizeof(int));
            }
        }
        printf("hillquag: recovered via Hill space = %s\n",
               best_conv == 0 ? "plain A..Z" : best_conv == 1 ? "KRYPTOS keyed" : "(fallback)");
    }

    if (!found) {
        // Fallback: full-matrix anneal over the period-P columnar IoC.
        anneal_matrix(cfg, cipher_indices, cipher_len, k, P, bestD);
        hill_mat_mul_blocks(bestD, k, cipher_indices, cipher_len, I);
        best_score = quag_strip_and_refine(cfg, I, cipher_len, pt_kw, ct_kw, P, ngram,
                                           crib_indices, crib_positions, n_cribs, cw, pt);
    }

    report_hill_quag(cfg, shared, cipher_indices, cipher_len, cribtext_str,
                     k, P, bestD, pt_kw, cw, pt, best_score, n_cribs);
}
