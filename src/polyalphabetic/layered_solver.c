#include "layered_solver.h"
#include "scoring.h"
#include "double_transposition_solver.h"
#include <string.h>
#include <math.h>

// =====================================================================
//  QUAG_TRANS: outer Quagmire III (keyed alphabet) o inner columnar
//  transposition (depth 0/1/2). See layered_solver.h for the model.
//
//  Pipeline:
//    1. build the keyed alphabet (KRYPTOS, pinned); determine period P.
//    2. STRIP the outer Quagmire by the monogram statistic
//       (derive_optimal_cycleword) -- transposition-invariant, so it needs
//       no knowledge of the inner transposition. Gives the intermediate I,
//       a columnar-transposed English text (with ~1/P of the columns still
//       mis-shifted, since each holds only ~N/P samples).
//    3. SOLVE the inner transposition on I by n-gram:
//         depth 0 -> I is already the plaintext (PK3, a plain Quagmire).
//         depth 1 -> a single columnar (exhaustive for small K, else anneal).
//         depth 2 -> a double columnar via dct_solve_core (divide & conquer).
//    4. REFINE the cycleword by n-gram coordinate-ascent THROUGH the fixed
//       inner transposition (the un-transposed text is now English-ordered,
//       so n-grams -- not just monograms -- resolve every column). One or two
//       outer rounds re-solve the transposition on the cleaned intermediate.
// =====================================================================

// Description of the recovered inner transposition, enough to re-apply it.
typedef struct {
    int depth;                              // 0 (none) / 1 (single) / 2 (double)
    int K,  order[MAX_COLS],  dir;          // depth 1
    int L1, order1[MAX_COLS], dir1;         // depth 2, first-applied (inner)
    int L2, order2[MAX_COLS], dir2;         // depth 2, second-applied (outer)
} InnerTrans;

// Un-apply the inner transposition: intermediate `in` (= Trans(PT)) -> `out` (PT).
// depth 0 copies through. Convention matches columnar_solver / dct: encryption
// applies stage 1 (inner) then stage 2 (outer), so decryption undoes 2 then 1.
static void inner_untranspose(const InnerTrans *t, const int *in, int len, int *out) {
    if (t->depth == 0) {
        memcpy(out, in, len * sizeof(int));
    } else if (t->depth == 1) {
        decrypt_columnar((int *) in, len, t->K, (int *) t->order, t->dir, out);
    } else {
        int tmp[MAX_CIPHER_LENGTH];
        decrypt_columnar((int *) in, len, t->L2, (int *) t->order2, t->dir2, tmp);
        decrypt_columnar(tmp,        len, t->L1, (int *) t->order1, t->dir1, out);
    }
}

// Composite cycleword: G component periods summed in keyed-index (shift) space.
// The effective shift at position i is sum_g sh[g][i mod len[g]] (mod 26); a single
// free period is G=1. This parameterization is what makes a long effective period
// (e.g. 45 = lcm(5,9)) recoverable -- each component shift is constrained by ~N/len
// samples instead of ~N/lcm, so the n-gram coordinate ascent has a real gradient.
typedef struct { int G; int len[8]; int P; int sh[8][MAX_CYCLEWORD_LEN]; } CycleKey;

// Synthesize the effective P-length cycleword (letter indices) from the components.
static void ck_synth(const CycleKey *ck, int ct_kw[], int cw_eff[]) {
    for (int i = 0; i < ck->P; i++) {
        int s = 0;
        for (int g = 0; g < ck->G; g++) s += ck->sh[g][i % ck->len[g]];
        cw_eff[i] = ct_kw[((s % g_alpha) + g_alpha) % g_alpha];
    }
}

// Build the component CycleKey from cfg. A composed multi-Quagmire (-cyclewordlens
// 5,9) yields those components warm-started to zero shift (the coordinate ascent has
// a strong per-component gradient); a single period is one component warm-started
// from the monogram strip mono_cw (P letter indices), which reduces to the plain
// free-period behaviour.
static void ck_from_cfg(ColossusConfig *cfg, int P, int ct_kw[], const int *mono_cw, CycleKey *ck) {
    int ct_lookup[ALPHABET_SIZE];
    for (int i = 0; i < g_alpha; i++) ct_lookup[ct_kw[i]] = i;
    memset(ck, 0, sizeof(*ck));
    ck->P = P;
    if (cfg->n_cycleword_lens > 1) {
        ck->G = cfg->n_cycleword_lens;
        for (int g = 0; g < ck->G; g++) ck->len[g] = cfg->cycleword_lens[g];
    } else {
        ck->G = 1; ck->len[0] = P;
        if (mono_cw) for (int j = 0; j < P; j++) ck->sh[0][j] = ct_lookup[mono_cw[j]];
    }
}

// Build the keyed PT/CT alphabets (Quagmire III => identical). KRYPTOS is pinned
// via -plaintextkeyword or -ciphertextkeyword. Returns false if neither is given.
static bool build_alphabet(ColossusConfig *cfg, int pt_kw[], int ct_kw[]) {
    if (cfg->user_plaintext_keyword_present)
        make_keyed_alphabet(cfg->user_plaintext_keyword, pt_kw);
    else if (cfg->user_ciphertext_keyword_present)
        make_keyed_alphabet(cfg->user_ciphertext_keyword, pt_kw);
    else
        return false;
    memcpy(ct_kw, pt_kw, ALPHABET_SIZE * sizeof(int));   // Quagmire III: PT-kw == CT-kw
    return true;
}

// Estimate the outer period by the periodic index of coincidence: pick the P in
// [2, Pmax] whose columns (positions == c mod P) have the highest mean IoC. Used
// only when -cyclewordlen is absent; the challenge supplies P directly.
static int estimate_period_ioc(const int *cipher, int len, int Pmax) {
    int best_P = 1; double best_ic = -1.0;
    for (int P = 2; P <= Pmax && P < len / 2; P++) {
        double sum_ic = 0.0; int ncols = 0;
        for (int c = 0; c < P; c++) {
            int hist[ALPHABET_SIZE]; memset(hist, 0, sizeof(hist));
            int n = 0;
            for (int i = c; i < len; i += P) { hist[cipher[i]]++; n++; }
            if (n < 2) continue;
            double s = 0.0;
            for (int k = 0; k < g_alpha; k++) s += (double) hist[k] * (hist[k] - 1);
            sum_ic += s / ((double) n * (n - 1)); ncols++;
        }
        if (ncols == 0) continue;
        double ic = sum_ic / ncols;
        // Prefer the SMALLEST period at the top IoC plateau (harmonics inflate IoC
        // from fewer samples), so require a clear margin to switch.
        if (ic > best_ic + 1e-6) { best_ic = ic; best_P = P; }
    }
    return best_P;
}

// Score a full decrypt (PT) by the configured fitness (n-gram + optional crib).
static double score_pt(ColossusConfig *cfg, const int *pt, int len,
                       float *ngram, int *crib_indices, int *crib_positions, int n_cribs) {
    float w_ngram = (cfg->weight_ngram != 0.f) ? cfg->weight_ngram : 1.f;
    return state_score((int *) pt, len, crib_indices, crib_positions, n_cribs,
                       ngram, cfg->ngram_size, w_ngram, cfg->weight_crib, 0.f, 0.f);
}

// ---- inner single columnar (depth 1) --------------------------------------
// Solve the columnar transposition of the (possibly slightly corrupted) English
// intermediate I. Exhaustive over K! for small K (guarantees the true order even
// when a few of I's letters are wrong); a short anneal for large K.

// One scored candidate column order (K columns, direction dir).
typedef struct { double score; int K; int dir; int order[MAX_COLS]; } ColCand;
static int colcand_desc(const void *a, const void *b) {
    double d = ((const ColCand *) b)->score - ((const ColCand *) a)->score;
    return (d < 0) ? -1 : (d > 0) ? 1 : 0;
}

// Forward declaration (defined below).
static double refine_components(ColossusConfig *cfg, const int *cipher, int len,
    int pt_kw[], int ct_kw[], CycleKey *ck, const InnerTrans *trans,
    float *ngram, int *ci, int *cp, int nc, int *pt_out, int max_sweeps);

// Append every K! column order for one (K, dir), each scored on the intermediate I
// (Heap's algorithm). Used to build the beam.
static void collect_orders(ColossusConfig *cfg, const int *I, int len, int K, int kk,
                           int *order, int dir, float *ngram, int *ci, int *cp, int nc,
                           ColCand *cands, int *ncand, int cap) {
    if (kk == 1) {
        if (*ncand < cap) {
            int pt[MAX_CIPHER_LENGTH];
            decrypt_columnar((int *) I, len, K, order, dir, pt);
            ColCand *c = &cands[(*ncand)++];
            c->score = score_pt(cfg, pt, len, ngram, ci, cp, nc);
            c->K = K; c->dir = dir; memcpy(c->order, order, K * sizeof(int));
        }
        return;
    }
    for (int i = 0; i < kk; i++) {
        collect_orders(cfg, I, len, K, kk - 1, order, dir, ngram, ci, cp, nc, cands, ncand, cap);
        int j = (kk % 2 == 0) ? i : 0;
        int t = order[j]; order[j] = order[kk - 1]; order[kk - 1] = t;
    }
}

// Short simulated-anneal over the column order (large K, where K! is infeasible).
static void anneal_order(ColossusConfig *cfg, const int *I, int len, int K, int dir,
                         float *ngram, int *ci, int *cp, int nc,
                         double *best, InnerTrans *out) {
    int restarts = cfg->n_restarts > 0 ? cfg->n_restarts : 200;
    int climbs   = cfg->n_hill_climbs > 0 ? cfg->n_hill_climbs : 2000;
    if (restarts > 400) restarts = 400;         // K is small enough that this is plenty
    int order[MAX_COLS], cur[MAX_COLS], pt[MAX_CIPHER_LENGTH];
    for (int r = 0; r < restarts; r++) {
        for (int c = 0; c < K; c++) cur[c] = c;
        shuffle(cur, K);
        decrypt_columnar((int *) I, len, K, cur, dir, pt);
        double cs = score_pt(cfg, pt, len, ngram, ci, cp, nc);
        for (int h = 0; h < climbs; h++) {
            memcpy(order, cur, K * sizeof(int));
            int a = rand_int(0, K), b = rand_int(0, K);
            int t = order[a]; order[a] = order[b]; order[b] = t;
            decrypt_columnar((int *) I, len, K, order, dir, pt);
            double s = score_pt(cfg, pt, len, ngram, ci, cp, nc);
            if (s >= cs) { cs = s; memcpy(cur, order, K * sizeof(int)); }
        }
        if (cs > *best) {
            *best = cs; out->depth = 1; out->K = K; out->dir = dir;
            memcpy(out->order, cur, K * sizeof(int));
        }
    }
}

// Beam budgets: rank all column orders on the (thin-strip) intermediate, then run the
// full component cycleword refine on the top BEAM_M0 -- only the true order refines to
// clean English, so it wins. (M0 must exceed the true order's rank on the raw strip,
// empirically the top few percent, so it is generous.)
#define BEAM_CAP   120000   // max column orders enumerated (K<=8: 8!*2 = 80640)
#define BEAM_M0      2500   // orders full-refined after the raw n-gram rank
#define BEAM_FULL      14   // refine sweeps (coordinate ascent early-exits at convergence)

// depth-1 solve. Ranks every column order on the monogram-stripped intermediate
// (cw_in = free-period monogram cycleword), then uses the COMPONENT cycleword refine
// -- which only reaches clean English through the true order -- to surface it. The
// refine warm-starts from ck_tmpl. Fills out_trans / out_ck / pt_out; returns score.
static double solve_depth1(ColossusConfig *cfg, const int *cipher, int len,
                           int pt_kw[], int ct_kw[], const int *cw_in, int P,
                           const CycleKey *ck_tmpl, float *ngram, int *ci, int *cp, int nc,
                           InnerTrans *out_trans, CycleKey *out_ck, int *pt_out) {
    int I[MAX_CIPHER_LENGTH];
    quagmire_decrypt(I, (int *) cipher, len, pt_kw, ct_kw, (int *) cw_in, P, cfg->variant);

    int Klo = cfg->min_cols > 2 ? cfg->min_cols : 2;
    int Khi = cfg->max_cols;
    int cap = len / 2;
    if (cap > MAX_COLS) cap = MAX_COLS;
    if (Khi > cap) Khi = cap;
    if (Klo > Khi) Klo = Khi;
    int dirs[2], ndir;
    if (cfg->read_direction == COL_READ_BOTH) { dirs[0] = COL_READ_TB; dirs[1] = COL_READ_BT; ndir = 2; }
    else { dirs[0] = cfg->read_direction; ndir = 1; }

    ColCand *cands = malloc(sizeof(ColCand) * BEAM_CAP);
    int ncand = 0, order[MAX_COLS];
    for (int K = Klo; K <= Khi && K <= 8; K++)
        for (int d = 0; d < ndir; d++) {
            for (int c = 0; c < K; c++) order[c] = c;
            collect_orders(cfg, I, len, K, K, order, dirs[d], ngram, ci, cp, nc, cands, &ncand, BEAM_CAP);
        }
    for (int K = 9; K <= Khi; K++)                     // K>8: annealed samples only
        for (int d = 0; d < ndir && ncand < BEAM_CAP; d++) {
            InnerTrans t; memset(&t, 0, sizeof(t)); double bs = -1e18;
            anneal_order(cfg, I, len, K, dirs[d], ngram, ci, cp, nc, &bs, &t);
            ColCand *c = &cands[ncand++];
            c->score = bs; c->K = K; c->dir = t.dir; memcpy(c->order, t.order, K * sizeof(int));
        }

    if (ncand == 0) { free(cands); out_trans->depth = 0; return -1e18; }
    qsort(cands, ncand, sizeof(ColCand), colcand_desc);
    int m0 = ncand < BEAM_M0 ? ncand : BEAM_M0;

    int tpt[MAX_CIPHER_LENGTH];
    double best = -1e18;
    for (int i = 0; i < m0; i++) {                     // full component refine each beam order
        InnerTrans t; memset(&t, 0, sizeof(t));
        t.depth = 1; t.K = cands[i].K; t.dir = cands[i].dir;
        memcpy(t.order, cands[i].order, t.K * sizeof(int));
        CycleKey ck = *ck_tmpl;
        double s = refine_components(cfg, cipher, len, pt_kw, ct_kw, &ck, &t,
                                     ngram, ci, cp, nc, tpt, BEAM_FULL);
        if (s > best) {
            best = s; *out_trans = t; *out_ck = ck;
            memcpy(pt_out, tpt, len * sizeof(int));
        }
    }
    free(cands);
    return best;
}

// ---- cycleword refine (n-gram coordinate ascent through a fixed transposition) --
// Coordinate-ascend every COMPONENT shift (sum_g len[g] of them, each 0..25): for
// each, try all 26 values and keep the one giving the best final-plaintext n-gram
// after un-transposing. Iterate to convergence. Because a component shift touches
// ~N/len[g] positions (not ~N/P), the gradient is strong even when the effective
// period P is long. Returns the best score; ck and *pt_out are updated in place.
static double refine_components(ColossusConfig *cfg, const int *cipher, int len,
                               int pt_kw[], int ct_kw[], CycleKey *ck,
                               const InnerTrans *trans, float *ngram,
                               int *ci, int *cp, int nc, int *pt_out, int max_sweeps) {
    int cw[MAX_CYCLEWORD_LEN], I[MAX_CIPHER_LENGTH], pt[MAX_CIPHER_LENGTH];
    ck_synth(ck, ct_kw, cw);
    quagmire_decrypt(I, (int *) cipher, len, pt_kw, ct_kw, cw, ck->P, cfg->variant);
    inner_untranspose(trans, I, len, pt);
    double best = score_pt(cfg, pt, len, ngram, ci, cp, nc);
    memcpy(pt_out, pt, len * sizeof(int));

    for (int sweep = 0; sweep < max_sweeps; sweep++) {
        bool improved = false;
        for (int g = 0; g < ck->G; g++) {
            for (int j = 0; j < ck->len[g]; j++) {
                double col_best = best; int col_keep = ck->sh[g][j];
                for (int s = 0; s < g_alpha; s++) {
                    ck->sh[g][j] = s;
                    ck_synth(ck, ct_kw, cw);
                    quagmire_decrypt(I, (int *) cipher, len, pt_kw, ct_kw, cw, ck->P, cfg->variant);
                    inner_untranspose(trans, I, len, pt);
                    double sc = score_pt(cfg, pt, len, ngram, ci, cp, nc);
                    if (sc > col_best) { col_best = sc; col_keep = s; }
                }
                ck->sh[g][j] = col_keep;
                if (col_best > best + 1e-9) { best = col_best; improved = true; }
            }
        }
        if (!improved) break;
    }
    ck_synth(ck, ct_kw, cw);
    quagmire_decrypt(I, (int *) cipher, len, pt_kw, ct_kw, cw, ck->P, cfg->variant);
    inner_untranspose(trans, I, len, pt_out);
    return best;
}

// ---- reporting ------------------------------------------------------------
static void report_quag_trans(ColossusConfig *cfg, SharedData *shared,
        const int *cipher, int len, char *cribtext_str,
        int pt_kw[], int ct_kw[], const CycleKey *ck, const InnerTrans *trans,
        const int *pt, double score, int n_cribs) {
    int P = ck->P;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(pt[i]);
    plaintext_string[len] = '\0';

    int n_words = 0;
    if (cfg->dictionary_present && shared->dict != NULL)
        n_words = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);

    if (trans->depth == 0)
        printf("\nquagtrans: Quagmire III (period %d), no transposition\n", P);
    else if (trans->depth == 1)
        printf("\nquagtrans: Quagmire III (period %d) o single columnar (K=%d, dir=%s)\n",
            P, trans->K, trans->dir == COL_READ_BT ? "bt" : "tb");
    else
        printf("\nquagtrans: Quagmire III (period %d) o double columnar (K1=%d, K2=%d)\n",
            P, trans->L1, trans->L2);

    printf("\nResult Score: %.2f | Words: %d\n", score, n_words);
    print_text((int *) cipher, len); printf("\n");
    print_text((int *) pt, len);     printf("\n");
    print_spaces_line(g_spaces_table, (int *) pt, len);
    if (n_cribs > 0) printf("%s\n", cribtext_str);

    int cw[MAX_CYCLEWORD_LEN];
    ck_synth(ck, ct_kw, cw);
    printf("keyed alphabet: "); print_text(pt_kw, g_alpha); printf("\n");
    if (ck->G > 1) {
        for (int g = 0; g < ck->G; g++) {
            printf("component cycleword %d (len %d):", g + 1, ck->len[g]);
            for (int j = 0; j < ck->len[g]; j++) printf(" %c", index_to_char(ct_kw[ck->sh[g][j]]));
            printf("\n");
        }
    }
    printf("effective cycleword (len %d):", P);
    for (int i = 0; i < P; i++) printf(" %c", index_to_char(cw[i]));
    printf("\n");
    if (trans->depth == 1) {
        printf("columnar (K=%d, dir=%s) order:", trans->K, trans->dir == COL_READ_BT ? "bt" : "tb");
        for (int c = 0; c < trans->K; c++) printf(" %d", trans->order[c]);
        printf("\n");
    } else if (trans->depth == 2) {
        printf("stage1 (K=%d, dir=%s) order:", trans->L1, trans->dir1 == COL_READ_BT ? "bt" : "tb");
        for (int c = 0; c < trans->L1; c++) printf(" %d", trans->order1[c]);
        printf("\nstage2 (K=%d, dir=%s) order:", trans->L2, trans->dir2 == COL_READ_BT ? "bt" : "tb");
        for (int c = 0; c < trans->L2; c++) printf(" %d", trans->order2[c]);
        printf("\n");
    }
    printf("\n");

    if (cfg->dictionary_present) printf(">>> %.2f, %d, %d, ", score, n_words, cfg->cipher_type);
    else                         printf(">>> %.2f, %d, ", score, cfg->cipher_type);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text((int *) cipher, len); printf(", ");
    print_text((int *) pt, len);     printf("\n");
    fflush(stdout);
}

// ---- entry ----------------------------------------------------------------
void solve_quag_trans(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str;

    if (cipher_len < 16) {
        printf("\n\nERROR: ciphertext too short for a layered quag/transposition solve.\n\n");
        return;
    }

    int pt_kw[ALPHABET_SIZE], ct_kw[ALPHABET_SIZE];
    if (!build_alphabet(cfg, pt_kw, ct_kw)) {
        printf("\n\nERROR: quagtrans needs the keyed alphabet pinned "
               "(-plaintextkeyword KRYPTOS).\n\n");
        return;
    }

    int P = cfg->cycleword_len_present ? cfg->cycleword_len
                                       : estimate_period_ioc(cipher_indices, cipher_len,
                                                             cfg->max_cycleword_len);
    if (P < 1) P = 1;
    printf("quagtrans: keyed alphabet pinned, outer Quagmire period P=%d%s, depth=%d\n",
           P, cfg->cycleword_len_present ? "" : " (estimated)", cfg->trans_depth);

    int depth = cfg->trans_depth;
    if (depth < 0) depth = 0;
    if (depth > 2) depth = 2;

    // ---- Phase 1: monogram strip of the outer Quagmire (transposition-invariant).
    //      Used to rank the inner transposition; the cycleword itself is refined over
    //      its components afterwards. ----
    int cw[MAX_CYCLEWORD_LEN];
    derive_optimal_cycleword(cfg, cipher_indices, cipher_len, pt_kw, ct_kw, cw, P, NULL);
    CycleKey ck;
    ck_from_cfg(cfg, P, ct_kw, cw, &ck);

    float *ngram = shared->ngram_data;
    int I[MAX_CIPHER_LENGTH], pt[MAX_CIPHER_LENGTH];
    InnerTrans trans; memset(&trans, 0, sizeof(trans)); trans.depth = depth;
    double best_score = -1e18;

    if (depth == 0) {
        // ---- PK3: plain Quagmire, no transposition. Component refine. ----
        best_score = refine_components(cfg, cipher_indices, cipher_len, pt_kw, ct_kw, &ck,
                                       &trans, ngram, crib_indices, crib_positions, n_cribs,
                                       pt, BEAM_FULL);
    } else if (depth == 1) {
        // ---- PK4: single columnar under the Quagmire. Beam over the column order. ----
        best_score = solve_depth1(cfg, cipher_indices, cipher_len, pt_kw, ct_kw, cw, P, &ck,
                                  ngram, crib_indices, crib_positions, n_cribs, &trans, &ck, pt);
    } else {
        // ---- PK6: double columnar. The small period (e.g. 6) strips reliably, so the
        //      divide-and-conquer solve runs on the stripped intermediate; then refine. ----
        double *bg = dct_load_bigrams(cfg->ngram_file, cfg->ngram_size);
        if (!bg) { printf("\n\nERROR: could not build the bigram table.\n\n"); return; }
        // The small outer period (e.g. 6) strips reliably by monogram, so one DCT solve
        // on the stripped intermediate suffices; then the component refine polishes cw.
        int cwr[MAX_CYCLEWORD_LEN];
        ck_synth(&ck, ct_kw, cwr);
        quagmire_decrypt(I, cipher_indices, cipher_len, pt_kw, ct_kw, cwr, P, cfg->variant);
        DCTResult res;
        dct_solve_core(I, cipher_len, cfg, ngram, bg,
                       crib_indices, crib_positions, n_cribs, &res);
        trans.depth = 2;
        trans.L1 = res.L1; trans.dir1 = res.dir1; memcpy(trans.order1, res.order1, res.L1 * sizeof(int));
        trans.L2 = res.L2; trans.dir2 = res.dir2; memcpy(trans.order2, res.order2, res.L2 * sizeof(int));
        best_score = refine_components(cfg, cipher_indices, cipher_len, pt_kw, ct_kw,
                                       &ck, &trans, ngram, crib_indices, crib_positions,
                                       n_cribs, pt, BEAM_FULL);
        free(bg);
    }

    report_quag_trans(cfg, shared, cipher_indices, cipher_len, cribtext_str,
                      pt_kw, ct_kw, &ck, &trans, pt, best_score, n_cribs);
}

// Public helper (see layered_solver.h): monogram-strip + n-gram refine a plain
// period-P Quagmire (no transposition). Used by the Hill-Quag solver. Fills cw_out
// (the P effective cycleword letters) and pt_out; returns the n-gram fitness.
double quag_strip_and_refine(ColossusConfig *cfg, const int *cipher, int len,
    int pt_kw[], int ct_kw[], int P, float *ngram,
    int *crib_indices, int *crib_positions, int n_cribs,
    int cw_out[], int pt_out[]) {
    derive_optimal_cycleword(cfg, (int *) cipher, len, pt_kw, ct_kw, cw_out, P, NULL);
    CycleKey ck;                                     // single component of length P
    ck_from_cfg(cfg, P, ct_kw, cw_out, &ck);
    if (ck.G != 1) { ck.G = 1; ck.len[0] = P;        // Hill-Quag inner Quag is one component
        for (int j = 0; j < P; j++) { int s=0; for (int q=0;q<g_alpha;q++) if (ct_kw[q]==cw_out[j]) s=q; ck.sh[0][j]=s; } }
    InnerTrans t; memset(&t, 0, sizeof(t)); t.depth = 0;
    double sc = refine_components(cfg, cipher, len, pt_kw, ct_kw, &ck, &t, ngram,
                                  crib_indices, crib_positions, n_cribs, pt_out, BEAM_FULL);
    ck_synth(&ck, ct_kw, cw_out);
    return sc;
}
