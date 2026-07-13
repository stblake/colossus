#include "double_transposition_solver.h"
#include "scoring.h"

// =====================================================================
//  Double columnar transposition -- divide and conquer (IDP)
//
//  See double_transposition_solver.h for the method. The pipeline is:
//    1. build a 26x26 bigram log-freq table (dct_load_bigrams);
//    2. SCREEN every (L1, L2, dir2) length/direction hypothesis with a
//       short IDP hill-climb, ranking them by best IDP;
//    3. REFINE the top few with a longer IDP hill-climb on K2;
//    4. FINISH each: undo K2, then break the residual single columnar
//       transposition (K1) with a quadgram-scored hill-climb, warm-started
//       from the IDP's own greedy column reconstruction;
//    5. keep the plaintext with the best n-gram score across all configs.
//
//  Everything is single-threaded and deterministic under a fixed -seed.
// =====================================================================

#define DCT_NEG 1.0e30          // a very negative sentinel for max searches
#define DCT_MIN_ROWS 2          // need >= this many rows to score column bigrams

// ---- bigram table --------------------------------------------------------

double *dct_load_bigrams(const char *ngram_file, int ngram_size) {
    if (ngram_size < 2) return NULL;
    FILE *fp = fopen(ngram_file, "r");
    if (!fp) return NULL;

    // Accumulate raw bigram counts by marginalizing every internal bigram
    // position of each listed n-gram (a quadgram ABCD contributes AB, BC, CD).
    double *cnt = malloc(26 * 26 * sizeof(double));
    for (int i = 0; i < 26 * 26; i++) cnt[i] = 0.0;

    char ngram[MAX_NGRAM_SIZE];
    int freq;
    while (fscanf(fp, "%s\t%d", ngram, &freq) == 2) {
        if ((int) strlen(ngram) < ngram_size) continue;
        for (int p = 0; p + 1 < ngram_size; p++) {
            int a = toupper((unsigned char) ngram[p])     - 'A';
            int b = toupper((unsigned char) ngram[p + 1]) - 'A';
            if (a < 0 || a >= 26 || b < 0 || b >= 26) continue;
            cnt[a * 26 + b] += freq;
        }
    }
    fclose(fp);

    double total = 0.0;
    for (int i = 0; i < 26 * 26; i++) total += cnt[i];
    if (total <= 0.0) { free(cnt); return NULL; }

    // log10 probabilities with a rare-but-not-impossible floor for unseen pairs
    // (the same shape as load_ngrams' -logprob floor).
    double floor = log10(0.01 / total);
    double *bg = cnt;                       // reuse the buffer in place
    for (int i = 0; i < 26 * 26; i++)
        bg[i] = (cnt[i] > 0.0) ? log10(cnt[i] / total) : floor;
    return bg;
}

// ---- greedy max-weight Hamiltonian path over the neighbour matrix --------
//
// B[i*L + j] = score of column j being the immediate RIGHT neighbour of column i.
// Repeatedly take the highest-scoring edge whose tail has no right neighbour yet
// and whose head has no left neighbour yet, skipping edges that would close a
// cycle, until L-1 edges form one Hamiltonian path. Returns the summed weight
// (the IDP) and, if chain_out != NULL, the head-to-tail node order.
static double dct_greedy_path(const double *B, int L, int *chain_out) {
    static int outdeg[MAX_COLS], indeg[MAX_COLS];
    static int nxt[MAX_COLS], prv[MAX_COLS];
    static int seg_start[MAX_COLS], seg_end[MAX_COLS];
    for (int i = 0; i < L; i++) {
        outdeg[i] = indeg[i] = 0;
        nxt[i] = prv[i] = -1;
        seg_start[i] = seg_end[i] = i;
    }

    double total = 0.0;
    for (int e = 0; e < L - 1; e++) {
        double best = -DCT_NEG;
        int bi = -1, bj = -1;
        for (int i = 0; i < L; i++) {
            if (outdeg[i]) continue;                 // i already has a right neighbour
            for (int j = 0; j < L; j++) {
                if (j == i || indeg[j]) continue;    // j already has a left neighbour
                if (seg_start[i] == seg_start[j]) continue;  // same segment -> cycle
                if (B[i * L + j] > best) { best = B[i * L + j]; bi = i; bj = j; }
            }
        }
        if (bi < 0) break;                           // no legal edge (shouldn't happen)
        nxt[bi] = bj; prv[bj] = bi; outdeg[bi] = 1; indeg[bj] = 1;
        int a = seg_start[bi], b = seg_end[bj];      // join segment(bi) .. segment(bj)
        seg_end[a] = b; seg_start[b] = a;
        total += best;
    }

    if (chain_out) {
        int head = 0;
        for (int i = 0; i < L; i++) if (prv[i] == -1) { head = i; break; }
        int c = 0;
        for (int u = head; u != -1 && c < L; u = nxt[u]) chain_out[c++] = u;
        while (c < L) chain_out[c++] = 0;            // paranoia: pad a broken chain
    }
    return total;
}

// ---- the Index of Digraphic Potential ------------------------------------

double dct_idp(const int *cipher, int len, int L1, int L2, int dir2,
               const int *order2, const double *bg, int *work, int *chain_out) {
    if (L1 < 2 || L2 < 2 || L1 >= len) return -DCT_NEG;

    // 1. undo the second transposition -> intermediate text I (in `work`).
    int *I = work;
    decrypt_columnar((int *) cipher, len, L2, (int *) order2, dir2, I);

    // 2. geometry of the hypothesized single (K1) transposition of I.
    int q = len / L1;          // short-column height (== number of scored rows)
    int s = len % L1;          // number of tall columns (leftmost s in grid order)
    if (q < DCT_MIN_ROWS) return -DCT_NEG;

    // Floor score for a column edge with no scorable letter pair (both cells non-letter
    // in every paired row). Non-letters (spaces / punctuation carried as negative
    // sentinels when the transposition grid includes them, per index_to_char) are not
    // digraphs, so such an edge is neutrally LOW rather than the artificial 0 a raw
    // log-frequency sum would give. Table minimum; scanned once per IDP eval (676 ops,
    // negligible beside the O(L1^2 * q) column-pair scoring below).
    double bg_floor = 0.0;
    for (int i = 0; i < 26 * 26; i++) if (bg[i] < bg_floor) bg_floor = bg[i];

    // Per read-position start windows. Read-position i (the i-th column segment as
    // it appears in I) starts at i*q + t, where t = number of tall columns among
    // the i preceding it, bounded by how many tall/short columns remain.
    static int starts[MAX_COLS][MAX_COLS];
    static int nstart[MAX_COLS];
    for (int i = 0; i < L1; i++) {
        int tlo = i - (L1 - s); if (tlo < 0) tlo = 0;
        int thi = (i < s) ? i : s;
        int n = 0;
        for (int t = tlo; t <= thi; t++) starts[i][n++] = i * q + t;
        nstart[i] = n;
    }

    // 3. B[i][j] = best MEAN bigram log-frequency of column j placed to the RIGHT of
    // column i (pair row k of each, k = 0..q-1), maximized over both placements. Any
    // paired row touching a non-letter cell (a space / punctuation sentinel) is not a
    // digraph and is skipped; the mean is taken over the LETTER pairs only. For an
    // all-letter cipher every row is valid (nv == q), so this is bit-identical to the
    // old per-row sum/q.
    static double B[MAX_COLS * MAX_COLS];
    for (int i = 0; i < L1; i++) {
        for (int j = 0; j < L1; j++) {
            if (i == j) { B[i * L1 + j] = -DCT_NEG; continue; }
            double best = -DCT_NEG;
            for (int a = 0; a < nstart[i]; a++) {
                const int *ci = I + starts[i][a];
                for (int b = 0; b < nstart[j]; b++) {
                    const int *cj = I + starts[j][b];
                    double sum = 0.0; int nv = 0;
                    for (int k = 0; k < q; k++) {
                        int u = ci[k], v = cj[k];
                        if (u < 0 || v < 0) continue;   // space/punct: not a digraph
                        sum += bg[u * 26 + v]; nv++;
                    }
                    double val = (nv > 0) ? sum / nv : bg_floor;
                    if (val > best) best = val;
                }
            }
            B[i * L1 + j] = best;
        }
    }

    // 4. greedy digraph-adjacency reconstruction. Normalize by the edge count so the
    // score is a MEAN per-neighbour digraph potential -- comparable across different
    // L1 hypotheses (a raw sum has one term per column, so it is trivially larger for
    // smaller L1 and would bias config ranking toward too-few columns).
    return dct_greedy_path(B, L1, chain_out) / (L1 - 1);
}

// ---- crib-aware K2 fitness -----------------------------------------------
//
// The IDP is structural (no plaintext), so a crib cannot guide it directly. But once
// K2 is undone the residual is a single columnar under K1, which the IDP's greedy
// chain reconstructs UP TO a rotation -- so we can CHEAPLY probe how well a putative
// K2 leads to the crib: rotation-resolve K1 from the chain, decrypt, and take the best
// crib match over the L1 rotations (+ reversal). Folding that into the K2 objective
// makes a known-plaintext fragment pull the K2 climb toward the true key: far from the
// solution every K2's best crib fraction sits at the ~1/26 floor (no gradient, so the
// IDP drives globally), but as K2 approaches the truth the crib fraction climbs steeply
// and sharpens onto the exact permutation. crib_w scales the bonus (0 => pure IDP).
static double dct_k2_fitness(const int *cipher, int len, int L1, int L2, int dir2,
                             const int *order2, const double *bg, int *work,
                             int *crib_idx, int *crib_pos, int n_cribs, double crib_w) {
    if (n_cribs == 0)
        return dct_idp(cipher, len, L1, L2, dir2, order2, bg, work, NULL);

    // dct_idp leaves the intermediate text I in `work` and fills `chain`.
    int chain[MAX_COLS];
    double idp = dct_idp(cipher, len, L1, L2, dir2, order2, bg, work, chain);
    const int *I = work;
    static int dec[MAX_CIPHER_LENGTH];
    int order[MAX_COLS], rchain[MAX_COLS];
    for (int g = 0; g < L1; g++) rchain[g] = chain[L1 - 1 - g];
    double best_crib = 0.0;
    for (int pass = 0; pass < 2; pass++) {
        const int *ch = pass ? rchain : chain;
        for (int rot = 0; rot < L1; rot++) {
            for (int g = 0; g < L1; g++) order[ch[(g + rot) % L1]] = g;
            decrypt_columnar((int *) I, len, L1, order, dir2, dec);
            double cf = crib_score(dec, len, crib_idx, crib_pos, n_cribs);
            if (cf > best_crib) best_crib = cf;
        }
    }
    return idp + crib_w * best_crib;
}

// ---- permutation move + simulated-annealing helpers ----------------------

static void dct_perm_seed(int *order, int K) {
    for (int c = 0; c < K; c++) order[c] = c;
    shuffle(order, K);
}

// One permutation-preserving move. Uses the paper's SEGMENT-WISE transformations
// (segment slide, segment swap, 3-partite rotation) -- which Lasry et al. found far
// more effective on transposition keys than single-element swaps -- plus a single
// swap for fine adjustment. Segments are UNCAPPED (up to K/2), unlike the old size-8
// block moves, so the climb can cross the large basins a swap cannot.
static void dct_perm_move(int *order, int K) {
    if (K < 2) return;
    double r = frand();
    if (r < 0.35) {
        // Single swap (fine move).
        int a = rand_int(0, K), b = rand_int(0, K);
        int t = order[a]; order[a] = order[b]; order[b] = t;
    } else if (r < 0.65) {
        // Segment slide: cut a contiguous segment and re-insert it elsewhere.
        int len = rand_int(1, K / 2 + 1);
        int s = rand_int(0, K - len + 1);
        int d = rand_int(0, K - len + 1);
        if (d == s) return;
        int tmp[MAX_COLS];
        for (int a = 0; a < len; a++) tmp[a] = order[s + a];
        if (d < s) { for (int a = s - 1; a >= d; a--) order[a + len] = order[a]; }
        else       { for (int a = s + len; a < d + len; a++) order[a - len] = order[a]; }
        for (int a = 0; a < len; a++) order[d + a] = tmp[a];
    } else if (r < 0.85) {
        // Segment swap: exchange two equal-length non-overlapping segments.
        int len = rand_int(1, K / 2 + 1);
        if (2 * len > K) return;
        int a = rand_int(0, K - 2 * len + 1);
        int b = rand_int(a + len, K - len + 1);
        for (int t = 0; t < len; t++) { int u = order[a + t]; order[a + t] = order[b + t]; order[b + t] = u; }
    } else {
        // 3-partite rotation: 3-cycle three distinct positions (a<-c<-b<-a).
        int a = rand_int(0, K), b = rand_int(0, K), c = rand_int(0, K);
        if (a == b || b == c || a == c) return;
        int t = order[a]; order[a] = order[c]; order[c] = order[b]; order[b] = t;
    }
}

// Left-to-right best-improvement heuristic (Lasry et al. sec 4.6, step 2): for each
// position i, swap it with the position j>i that most improves the IDP; repeat full
// passes until no swap helps. A strong deterministic local optimizer for the K2 order,
// applied both to each restart's seed and as a final polish -- it drives a promising
// key to its swap-local optimum far faster than annealing alone. O(L^2) IDP evals per
// pass (a few hundred at L~23), negligible beside the anneal.
static double dct_lr_improve(const int *cipher, int len, int L1, int L2, int dir2,
                             const double *bg, int *order, int *work, double cur,
                             int *crib_idx, int *crib_pos, int n_cribs, double crib_w) {
    int improved = 1;
    while (improved) {
        improved = 0;
        for (int i = 0; i < L2 - 1; i++) {
            int best_j = -1; double best_gain = 1e-9;
            for (int j = i + 1; j < L2; j++) {
                int t = order[i]; order[i] = order[j]; order[j] = t;
                double s = dct_k2_fitness(cipher, len, L1, L2, dir2, order, bg, work,
                                          crib_idx, crib_pos, n_cribs, crib_w);
                t = order[i]; order[i] = order[j]; order[j] = t;    // undo
                if (s - cur > best_gain) { best_gain = s - cur; best_j = j; }
            }
            if (best_j >= 0) {
                int t = order[i]; order[i] = order[best_j]; order[best_j] = t;
                cur += best_gain; improved = 1;
            }
        }
    }
    return cur;
}

// Geometric temperature at step `it` of `iters` from hi -> lo.
static double dct_temp(double hi, double lo, int it, int iters) {
    if (iters <= 1) return lo;
    double f = (double) it / (double) (iters - 1);
    return hi * pow(lo / hi, f);
}

// ---- phase A: hill-climb K2 by the IDP -----------------------------------
//
// Anneals order2 (a permutation of 0..L2-1) to MAXIMIZE dct_idp. Returns the best
// IDP found and writes the best order into best_order. If seed_order != NULL the
// first restart starts from it (warm restart during refinement). When verbose, each
// new GLOBAL best IDP is streamed live (tag/t0 label it) -- the K2-climb analogue of
// the engine's report_verbose best-improvement dialog.
static double dct_idp_climb(const int *cipher, int len, int L1, int L2, int dir2,
                            const double *bg, int restarts, int iters,
                            const int *seed_order, int *best_order, int *work,
                            int *crib_idx, int *crib_pos, int n_cribs, double crib_w,
                            int verbose, const char *tag, clock_t t0) {
    int cur[MAX_COLS], cand[MAX_COLS];
    double global_best = -DCT_NEG;

    for (int rs = 0; rs < restarts; rs++) {
        if (rs == 0 && seed_order) { for (int c = 0; c < L2; c++) cur[c] = seed_order[c]; }
        else dct_perm_seed(cur, L2);

        double cur_score = dct_k2_fitness(cipher, len, L1, L2, dir2, cur, bg, work,
                                          crib_idx, crib_pos, n_cribs, crib_w);
        // Left-to-right local optimization of the seed (paper step 2) before annealing.
        cur_score = dct_lr_improve(cipher, len, L1, L2, dir2, bg, cur, work, cur_score,
                                   crib_idx, crib_pos, n_cribs, crib_w);
        double local_best = cur_score;
        int local_best_order[MAX_COLS];
        for (int c = 0; c < L2; c++) local_best_order[c] = cur[c];

        // Temperature matched to the PER-EDGE-normalized IDP scale (values ~ -2.3,
        // single-swap deltas ~ 0.02-0.2), so early iterations still climb rather than
        // random-walk. (A schedule tuned to the old un-normalized sum would be ~L1x
        // too hot here and waste the whole budget wandering.)
        double hi = 0.10, lo = 0.004;

        for (int it = 0; it < iters; it++) {
            for (int c = 0; c < L2; c++) cand[c] = cur[c];
            dct_perm_move(cand, L2);
            double cand_score = dct_k2_fitness(cipher, len, L1, L2, dir2, cand, bg, work,
                                               crib_idx, crib_pos, n_cribs, crib_w);
            double d = cand_score - cur_score;
            double temp = dct_temp(hi, lo, it, iters);
            if (d > 0 || frand() < exp(d / temp)) {
                for (int c = 0; c < L2; c++) cur[c] = cand[c];
                cur_score = cand_score;
                if (cur_score > local_best) {
                    local_best = cur_score;
                    for (int c = 0; c < L2; c++) local_best_order[c] = cur[c];
                }
            }
        }
        // Polish the restart's best with the left-to-right heuristic (paper step 3/4).
        double polished = dct_lr_improve(cipher, len, L1, L2, dir2, bg,
                                         local_best_order, work, local_best,
                                         crib_idx, crib_pos, n_cribs, crib_w);
        if (polished > local_best) local_best = polished;

        if (local_best > global_best) {
            global_best = local_best;
            for (int c = 0; c < L2; c++) best_order[c] = local_best_order[c];
            if (verbose)
                printf("    %s IDP=%.4f  (restart %d/%d)  %.0fs\n", tag, global_best,
                    rs + 1, restarts, (double)(clock() - t0) / CLOCKS_PER_SEC),
                fflush(stdout);
        }
    }
    return global_best;
}

// ---- phase B: break the residual single columnar (K1) by n-grams ---------
//
// Anneals order1 (a permutation of 0..L1-1) to MAXIMIZE the quadgram score of
// decrypt_columnar(I, L1, order1, dir1). Restart 0 warm-starts from seed_order
// (the IDP's greedy reconstruction), the rest are random. Returns the best score
// and writes the best plaintext + order.
static double dct_columnar_finish(const int *I, int len, int L1, int dir1,
                                  const int *seed_order, int restarts, int iters,
                                  float *ngram, int ngram_size,
                                  int *crib_idx, int *crib_pos, int n_cribs,
                                  float w_ngram, float w_crib,
                                  int *best_pt, int *best_order) {
    int cur[MAX_COLS], cand[MAX_COLS];
    int dec[MAX_CIPHER_LENGTH];
    double global_best = -DCT_NEG;

    for (int rs = 0; rs < restarts; rs++) {
        if (rs == 0 && seed_order) { for (int c = 0; c < L1; c++) cur[c] = seed_order[c]; }
        else dct_perm_seed(cur, L1);

        decrypt_columnar((int *) I, len, L1, cur, dir1, dec);
        // state_score == ngram_score when n_cribs==0; else blends the crib match (over
        // plaintext positions, which decrypt_columnar's output indexes directly).
        double cur_score = state_score(dec, len, crib_idx, crib_pos, n_cribs,
                                       ngram, ngram_size, w_ngram, w_crib, 0.f, 0.f);
        double hi = 0.10, lo = 0.002;

        for (int it = 0; it < iters; it++) {
            for (int c = 0; c < L1; c++) cand[c] = cur[c];
            dct_perm_move(cand, L1);
            decrypt_columnar((int *) I, len, L1, cand, dir1, dec);
            double cand_score = state_score(dec, len, crib_idx, crib_pos, n_cribs,
                                            ngram, ngram_size, w_ngram, w_crib, 0.f, 0.f);
            double d = cand_score - cur_score;
            double temp = dct_temp(hi, lo, it, iters);
            if (d > 0 || frand() < exp(d / temp)) {
                for (int c = 0; c < L1; c++) cur[c] = cand[c];
                cur_score = cand_score;
                if (cur_score > global_best) {
                    global_best = cur_score;
                    for (int c = 0; c < L1; c++) best_order[c] = cur[c];
                    decrypt_columnar((int *) I, len, L1, cur, dir1, best_pt);
                }
            }
        }
    }
    return global_best;
}

// The IDP greedy chain lists the read-positions in left-to-right column order, but a
// columnar transposition is recoverable only UP TO a cyclic rotation (a Hamiltonian
// path over the neighbour matrix has free endpoints, and the last column's spurious
// "best right neighbour" closes the chain into a near-cycle). So try all L1 rotations
// of the chain -- and its reversal, for the opposite read direction -- as columnar
// orders, decrypt each, and keep the best-scoring order1. The correct rotation cuts
// the near-cyclic chain at the true first/last-column boundary, recovering K1 almost
// exactly when K2 is correct (verified: an exact recovery on the true K2).
static double dct_rotations_best(const int *chain, int L1, const int *I, int len,
                                 int dir, float *ngram, int ngram_size,
                                 int *crib_idx, int *crib_pos, int n_cribs,
                                 float w_ngram, float w_crib, int *best_order) {
    static int dec[MAX_CIPHER_LENGTH];
    int order[MAX_COLS], rchain[MAX_COLS];
    for (int g = 0; g < L1; g++) rchain[g] = chain[L1 - 1 - g];
    double best = -DCT_NEG;
    for (int pass = 0; pass < 2; pass++) {
        const int *ch = pass ? rchain : chain;
        for (int rot = 0; rot < L1; rot++) {
            for (int g = 0; g < L1; g++) order[ch[(g + rot) % L1]] = g;
            decrypt_columnar((int *) I, len, L1, order, dir, dec);
            // state_score blends n-grams with the crib match; with n_cribs==0 it is
            // exactly ngram_score, so a crib-free solve is unchanged.
            double sc = state_score(dec, len, crib_idx, crib_pos, n_cribs,
                                    ngram, ngram_size, w_ngram, w_crib, 0.f, 0.f);
            if (sc > best) { best = sc; for (int c = 0; c < L1; c++) best_order[c] = order[c]; }
        }
    }
    return best;
}

// ---- config ranking scaffold ---------------------------------------------

typedef struct { int L1, L2, dir2; double idp; int order2[MAX_COLS]; } DCConfig;

static int dct_cmp_config(const void *a, const void *b) {
    double da = ((const DCConfig *) a)->idp, db = ((const DCConfig *) b)->idp;
    return (da < db) - (da > db);          // descending IDP
}

// ---- divide-and-conquer search core (no I/O) -----------------------------

void dct_solve_core(const int *cipher, int len, ColossusConfig *cfg,
                    float *ngram_data, const double *bg,
                    int *crib_indices, int *crib_positions, int n_cribs,
                    DCTResult *res) {
    res->len = len;
    res->score = -DCT_NEG;
    float w_ngram = cfg->weight_ngram, w_crib = cfg->weight_crib;
    // Crib bonus weight for the K2 climb: a full crib match (fraction 1) adds ~1.0,
    // well above the per-edge IDP spread (~0.3), so a crib sharpens the K2 climb onto
    // the exact key near the solution while the IDP still drives the search globally.
    double crib_w = 1.0;

    // Column-count range (shared with the columnar solvers), clamped to [2, len/2].
    int lo = cfg->min_cols, hi = cfg->max_cols;
    int cap = len / 2; if (cap > MAX_COLS - 1) cap = MAX_COLS - 1;
    if (lo < 2) lo = 2; if (hi > cap) hi = cap; if (lo > hi) lo = hi;

    // Read directions to try (both stages share cfg->read_direction).
    int dirs[2], ndir;
    if (cfg->read_direction == COL_READ_BOTH) { dirs[0] = COL_READ_TB; dirs[1] = COL_READ_BT; ndir = 2; }
    else { dirs[0] = cfg->read_direction; ndir = 1; }

    // Budgets. Screening is broad but shallow; refinement is deep on the top configs.
    // -nrestarts / -nhillclimbs scale the refinement stage; screening depth is tied to
    // it (a fraction), so a bigger run also ranks configs more reliably -- important at
    // long key lengths where the per-edge IDP margin is thin and the true (L1,L2) must
    // survive screening to reach the deep refinement.
    int refine_restarts = cfg->n_restarts, refine_iters = cfg->n_hill_climbs;
    if (refine_restarts < 8)    refine_restarts = 30;      // DC default if user left it tiny
    if (refine_iters   < 2000)  refine_iters   = 8000;
    int SCREEN_RESTARTS = refine_restarts / 4; if (SCREEN_RESTARTS < 6)    SCREEN_RESTARTS = 6;
    int SCREEN_ITERS    = refine_iters / 5;    if (SCREEN_ITERS   < 1500)  SCREEN_ITERS   = 1500;
    // The finish is warm-started from the rotation-resolved IDP chain (near-exact for
    // a correct K2), so it needs only a light polish anneal.
    const int FINISH_RESTARTS = 4, FINISH_ITERS = 8000;
    // Keep several top configs: the per-edge IDP margin between the true (L1,L2) and
    // spurious neighbours is thin, so the finish's real n-gram score -- not the
    // screening IDP alone -- makes the final call across the kept set.
    int keep = 8;

    static int work[MAX_CIPHER_LENGTH];
    int nconf = (hi - lo + 1) * (hi - lo + 1) * ndir;
    if (nconf < 1) nconf = 1;
    DCConfig *configs = malloc((size_t) nconf * sizeof(DCConfig));
    int nc = 0;

    clock_t t_start = clock();
    if (cfg->verbose) {
        printf("\ntranscol2-dc: screening %d (L1,L2,dir) configs, L in [%d,%d], "
               "screen %dx%d, refine %dx%d, keep %d\n",
               nconf, lo, hi, SCREEN_RESTARTS, SCREEN_ITERS,
               refine_restarts, refine_iters, keep);
        fflush(stdout);
    }

    // --- SCREEN: rank every (L1, L2, dir2) by a short IDP climb ---
    double screen_best = -DCT_NEG;
    for (int d = 0; d < ndir; d++)
        for (int L2 = lo; L2 <= hi; L2++)
            for (int L1 = lo; L1 <= hi; L1++) {
                DCConfig *cc = &configs[nc];
                cc->L1 = L1; cc->L2 = L2; cc->dir2 = dirs[d];
                cc->idp = dct_idp_climb(cipher, len, L1, L2, dirs[d],
                                        bg, SCREEN_RESTARTS, SCREEN_ITERS,
                                        NULL, cc->order2, work,
                                        crib_indices, crib_positions, n_cribs, crib_w,
                                        0, NULL, t_start);
                if (cc->idp > screen_best) screen_best = cc->idp;
                if (cfg->verbose) {
                    printf("  screen [%2d/%d] L1=%2d L2=%2d %s  IDP=%.4f  (best %.4f)  %.0fs\n",
                        nc + 1, nconf, L1, L2, dirs[d] == COL_READ_BT ? "bt" : "tb",
                        cc->idp, screen_best,
                        (double)(clock() - t_start) / CLOCKS_PER_SEC);
                    fflush(stdout);
                }
                nc++;
            }
    qsort(configs, nc, sizeof(DCConfig), dct_cmp_config);
    if (keep > nc) keep = nc;

    if (cfg->verbose) {
        printf("top %d screened configs (by IDP):\n", keep);
        for (int i = 0; i < keep; i++)
            printf("  #%d  L1=%d L2=%d dir2=%s  IDP=%.4f\n", i + 1,
                configs[i].L1, configs[i].L2,
                configs[i].dir2 == COL_READ_BT ? "bt" : "tb", configs[i].idp);
        fflush(stdout);
    }

    // --- REFINE + FINISH the top configs ---
    for (int i = 0; i < keep; i++) {
        int L1 = configs[i].L1, L2 = configs[i].L2, dir2 = configs[i].dir2;
        if (cfg->verbose) {
            printf("refine #%d/%d L1=%d L2=%d dir2=%s ...  %.0fs\n", i + 1, keep, L1, L2,
                dir2 == COL_READ_BT ? "bt" : "tb",
                (double)(clock() - t_start) / CLOCKS_PER_SEC);
            fflush(stdout);
        }
        int k2[MAX_COLS], chain[MAX_COLS];
        char tag[64];
        snprintf(tag, sizeof tag, "[refine #%d L1=%d L2=%d]", i + 1, L1, L2);
        dct_idp_climb(cipher, len, L1, L2, dir2,
                      bg, refine_restarts, refine_iters, configs[i].order2, k2, work,
                      crib_indices, crib_positions, n_cribs, crib_w,
                      cfg->verbose, tag, t_start);

        // Intermediate under the refined K2, plus the greedy column reconstruction.
        static int I[MAX_CIPHER_LENGTH];
        decrypt_columnar((int *) cipher, len, L2, k2, dir2, I);
        dct_idp(cipher, len, L1, L2, dir2, k2, bg, work, chain);

        for (int d = 0; d < ndir; d++) {
            int dir1 = dirs[d];
            static int pt[MAX_CIPHER_LENGTH];
            int o1[MAX_COLS], seed_order1[MAX_COLS];
            // Rotation-resolve the chain into a strong single-columnar warm start,
            // then polish with a short n-gram anneal.
            dct_rotations_best(chain, L1, I, len, dir1, ngram_data, cfg->ngram_size,
                               crib_indices, crib_positions, n_cribs, w_ngram, w_crib, seed_order1);
            double sc = dct_columnar_finish(I, len, L1, dir1, seed_order1,
                                            FINISH_RESTARTS, FINISH_ITERS,
                                            ngram_data, cfg->ngram_size,
                                            crib_indices, crib_positions, n_cribs, w_ngram, w_crib,
                                            pt, o1);
            int is_best = (sc > res->score);
            if (is_best) {
                res->score = sc;
                for (int t = 0; t < len; t++) res->plaintext[t] = pt[t];
                for (int c = 0; c < L1; c++) res->order1[c] = o1[c];
                for (int c = 0; c < L2; c++) res->order2[c] = k2[c];
                res->L1 = L1; res->L2 = L2; res->dir1 = dir1; res->dir2 = dir2;
            }
            // Live dialog (verbose only): on each new GLOBAL best, print a metrics block
            // in the same value-TAB-[label] convention as the polyalphabetic solver's
            // report_verbose, ending with the full decrypted plaintext.
            if (cfg->verbose && is_best) {
                double elapsed = (double)(clock() - t_start) / CLOCKS_PER_SEC;
                printf("\n%.2f\t[sec]\n", elapsed);
                printf("%d\t[K1 cols]\n", L1);
                printf("%d\t[K2 cols]\n", L2);
                printf("%s/%s\t[read dir K1/K2]\n",
                    dir1 == COL_READ_BT ? "bt" : "tb", dir2 == COL_READ_BT ? "bt" : "tb");
                printf("%.4f\t[IOC]\n", index_of_coincidence(pt, len));
                printf("%.4f\t[entropy]\n", entropy(pt, len));
                printf("%.2f\t[chi-squared]\n", chi_squared(pt, len));
                printf("%.2f\t[score]\n", sc);
                print_text(pt, len); printf("\n");
                fflush(stdout);
            }
        }
    }
    if (cfg->verbose) {
        printf("done in %.0fs. best score %.3f (K1=%d K2=%d)\n",
            (double)(clock() - t_start) / CLOCKS_PER_SEC, res->score, res->L1, res->L2);
        fflush(stdout);
    }
    free(configs);
}

// ---- solver entry --------------------------------------------------------

void solve_double_transposition(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str;

    if (cipher_len < 16) {
        printf("\n\nERROR: ciphertext too short for a double-transposition solve.\n\n");
        return;
    }

    double *bg = dct_load_bigrams(cfg->ngram_file, cfg->ngram_size);
    if (!bg) {
        printf("\n\nERROR: could not build the bigram table from '%s'.\n\n", cfg->ngram_file);
        return;
    }

    DCTResult res;
    dct_solve_core(cipher_indices, cipher_len, cfg, shared->ngram_data, bg,
                   crib_indices, crib_positions, n_cribs, &res);

    // --- report ---
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < cipher_len; i++) plaintext_string[i] = index_to_char(res.plaintext[i]);
    plaintext_string[cipher_len] = '\0';

    int n_words_found = 0;
    if (cfg->dictionary_present && shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);

    printf("\ntranscol2-dc: double columnar (divide & conquer), K1=%d x K2=%d cols\n",
        res.L1, res.L2);
    printf("\nResult Score: %.2f | Words: %d\n", res.score, n_words_found);

    print_text(cipher_indices, cipher_len); printf("\n");
    print_text(res.plaintext, cipher_len); printf("\n");
    print_spaces_line(g_spaces_table, res.plaintext, cipher_len);

    // Crib line + per-position match row (matches the columnar solver's report):
    // '_' where no crib, else |decrypted - crib| (a digit, or '*' if >= 10).
    if (n_cribs > 0) {
        printf("%s\n", cribtext_str);
        for (int i = 0; i < cipher_len; i++) {
            if (cribtext_str[i] == '_') { printf("_"); continue; }
            int diff = abs(res.plaintext[i] - g_char_to_idx[toupper((unsigned char)cribtext_str[i]) & 127]);
            if (diff < 10) printf("%d", diff); else printf("*");
        }
        printf("\n");
    }

    printf("K1 (L=%d, dir=%s) order:", res.L1, res.dir1 == COL_READ_BT ? "bt" : "tb");
    for (int c = 0; c < res.L1; c++) printf(" %d", res.order1[c]);
    printf("\nK2 (L=%d, dir=%s) order:", res.L2, res.dir2 == COL_READ_BT ? "bt" : "tb");
    for (int c = 0; c < res.L2; c++) printf(" %d", res.order2[c]);
    printf("\n\n");

    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, ", res.score, n_words_found, cfg->cipher_type);
    else
        printf(">>> %.2f, %d, ", res.score, cfg->cipher_type);
    printf("%d, %d, ", res.L1, res.L2);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(cipher_indices, cipher_len); printf(", ");
    print_text(res.plaintext, cipher_len); printf("\n");
    fflush(stdout);

    free(bg);
}
