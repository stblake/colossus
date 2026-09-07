#include "tridigital_solver.h"
#include "engine.h"
#include "scoring.h"
#include <math.h>

// =====================================================================
//  Tridigital solver -- see tridigital_solver.h
// =====================================================================

// Decode scratch (single-threaded per worker), off the stack.
static _Thread_local int g_tri_dec[MAX_CIPHER_LENGTH];               // a decode (letters + TRI_SPACE)
static _Thread_local unsigned char g_bp_sym[MAX_CIPHER_LENGTH][TRI_BEAM_MAX];   // 0..25 letter, 255 = space
static _Thread_local unsigned char g_bp_back[MAX_CIPHER_LENGTH][TRI_BEAM_MAX];  // parent beam index
// Parsed ciphertext digit stream, filled by solve_tridigital (read-only during the search).
static int g_tri_digits[MAX_CIPHER_LENGTH];

typedef struct {
    int   *digits;                            // ciphertext digit stream (length clen)
    int    clen;
    int    ngram_size;
    float *ngram;
    int    beam;                              // decode beam width during the search

    int    n_configs;                         // kept separator-digit configs
    int    cfg_sep[TRI_N_CONFIGS];            // separator digit per kept config (word-length ranked)
    int    warm[TRI_N_CONFIGS][TRI_MAP_CELLS];// mini-solve warm-start key per kept config

    // Per-config engine driving + gaming-resistant cross-config selection (cf. Monome-Dinome):
    // the engine is run once per kept config; tri_report stashes that config's best state here
    // (no printing while `quiet`), and solve_tridigital selects across configs by DICTIONARY
    // COVERAGE of the best decode -- a signal n-gram gaming cannot fake.
    int    active_config;
    bool   quiet;
    int    stash_sep;
    int    stash_key[TRI_MAP_CELLS];
    double stash_score;
    int    stash_decode[MAX_CIPHER_LENGTH];   // full decode (letters + TRI_SPACE), length clen
    int    stash_declen;                      // == clen

    // Dictionary word-set for the word-hit gradient + coverage selection (built once). NULL => none.
    char **wordset;
    int    wordset_cap;
    int    max_word_len;
} TridigitalScratch;

// ===================================================================
//  Dictionary word-set (FNV-1a open addressing; borrowed pointers into shared->dict)
// ===================================================================

static unsigned long tri_str_hash(const char *s) {          // FNV-1a
    unsigned long h = 1469598103934665603UL;
    for (; *s; s++) { h ^= (unsigned char) *s; h *= 1099511628211UL; }
    return h;
}

static char **tri_build_wordset(char **dict, int n_dict, int *cap_out) {
    if (!dict || n_dict <= 0) { *cap_out = 0; return NULL; }
    int cap = 1; while (cap < 4 * n_dict) cap <<= 1;         // <= 50% load
    char **tab = calloc((size_t) cap, sizeof(char *));
    if (!tab) { *cap_out = 0; return NULL; }
    for (int i = 0; i < n_dict; i++) {
        if ((int) strlen(dict[i]) < TRI_MIN_WORD) continue;
        unsigned long h = tri_str_hash(dict[i]) & (unsigned long) (cap - 1);
        while (tab[h]) h = (h + 1) & (unsigned long) (cap - 1);
        tab[h] = dict[i];
    }
    *cap_out = cap;
    return tab;
}

static bool tri_wordset_has(char **tab, int cap, const char *frag) {
    unsigned long h = tri_str_hash(frag) & (unsigned long) (cap - 1);
    while (tab[h]) { if (strcmp(tab[h], frag) == 0) return true; h = (h + 1) & (unsigned long) (cap - 1); }
    return false;
}

// ===================================================================
//  Inner decode: beam Viterbi over the ambiguous digit stream
// ===================================================================
//
// Given the key (group per letter) and the separator digit, each ciphertext digit yields <=3
// candidate letters (the separator yields a word break -> TRI_SPACE). We choose letters to
// maximize the n-gram score of the letter stream. Because ngram_score treats spaces as
// transparent (letters-only projection), the DP state -- the packed last (n-1) letters -- is
// carried ACROSS word breaks. A window value is added once, when the letter completing it is
// placed (i.e. from the n-th letter on), matching ngram_score's window walk exactly. A small
// beam keeps only the top-`beam` states per position (approximate for beam < 3^(n-1), exact
// otherwise). Fills out[0..clen-1] (letters + TRI_SPACE) and returns the letter count.
static int tri_decode(const TridigitalScratch *a, const int key[TRI_MAP_CELLS], int sep,
                      int beam, int out[]) {
    const int n = a->ngram_size;
    const int clen = a->clen;
    const float *ngram = a->ngram;
    if (beam > TRI_BEAM_MAX) beam = TRI_BEAM_MAX;

    long POW = 1; for (int i = 0; i < n - 1; i++) POW *= 26;   // 26^(n-1)

    // digit -> candidate letters. The 9 non-separator digits map (ascending) to groups 0..8.
    int digit_of_group[TRI_NGROUPS];
    { int r = 0; for (int d = 0; d < 10; d++) if (d != sep) digit_of_group[r++] = d; }
    int cand[10][TRI_GROUP_CAP], ncand[10];
    for (int d = 0; d < 10; d++) ncand[d] = 0;
    for (int L = 0; L < TRI_MAP_CELLS; L++) {
        int gg = key[L];
        if (gg < 0 || gg >= TRI_NGROUPS) continue;            // defensive
        int d = digit_of_group[gg];
        if (ncand[d] < TRI_GROUP_CAP) cand[d][ncand[d]++] = L;
    }

    long   cur_key[TRI_BEAM_MAX];
    double cur_score[TRI_BEAM_MAX];
    int    cur_cnt = 1;
    cur_key[0] = 0; cur_score[0] = 0.0;
    int have = 0;                                             // letters emitted so far

    for (int t = 0; t < clen; t++) {
        int d = a->digits[t];
        if (d < 0 || d > 9 || d == sep) {                     // word break: carry the beam unchanged
            for (int b = 0; b < cur_cnt; b++) { g_bp_sym[t][b] = 255; g_bp_back[t][b] = (unsigned char) b; }
            continue;
        }
        int nc = ncand[d];
        if (nc == 0) {                                        // no candidate (degenerate key): filler
            for (int b = 0; b < cur_cnt; b++) { g_bp_sym[t][b] = 0; g_bp_back[t][b] = (unsigned char) b; }
            have++;
            continue;
        }
        long   nk_top[TRI_BEAM_MAX];
        double ns_top[TRI_BEAM_MAX];
        int    from_top[TRI_BEAM_MAX], sym_top[TRI_BEAM_MAX], top = 0;
        int addwin = (have >= n - 1);
        for (int b = 0; b < cur_cnt; b++) {
            for (int c = 0; c < nc; c++) {
                int x = cand[d][c];
                long widx = cur_key[b] * 26 + x;
                double ns = cur_score[b] + (addwin ? ngram[widx] : 0.0);
                long nk = widx % POW;
                // Insert into the descending-by-score top list (beam is tiny, so O(beam) insert).
                if (top < beam) {
                    int p = top++;
                    while (p > 0 && ns_top[p - 1] < ns) {
                        ns_top[p] = ns_top[p - 1]; nk_top[p] = nk_top[p - 1];
                        from_top[p] = from_top[p - 1]; sym_top[p] = sym_top[p - 1]; p--;
                    }
                    ns_top[p] = ns; nk_top[p] = nk; from_top[p] = b; sym_top[p] = x;
                } else if (ns > ns_top[beam - 1]) {
                    int p = beam - 1;
                    while (p > 0 && ns_top[p - 1] < ns) {
                        ns_top[p] = ns_top[p - 1]; nk_top[p] = nk_top[p - 1];
                        from_top[p] = from_top[p - 1]; sym_top[p] = sym_top[p - 1]; p--;
                    }
                    ns_top[p] = ns; nk_top[p] = nk; from_top[p] = b; sym_top[p] = x;
                }
            }
        }
        for (int b = 0; b < top; b++) {
            cur_key[b] = nk_top[b]; cur_score[b] = ns_top[b];
            g_bp_sym[t][b]  = (unsigned char) sym_top[b];
            g_bp_back[t][b] = (unsigned char) from_top[b];
        }
        cur_cnt = top;
        have++;
    }

    // Best final beam entry, then backtrack to fill out[].
    int best = 0; for (int b = 1; b < cur_cnt; b++) if (cur_score[b] > cur_score[best]) best = b;
    int b = best, nletters = 0;
    for (int t = clen - 1; t >= 0; t--) {
        int sym = g_bp_sym[t][b];
        if (sym == 255) out[t] = TRI_SPACE;
        else { out[t] = sym; nletters++; }
        b = g_bp_back[t][b];
    }
    return nletters;
}

// Fraction of dictionary-word-eligible words (length >= TRI_MIN_WORD) that ARE dictionary
// words -- a cheap, sharp within-config gradient (word boundaries are known from the separator).
static double tri_word_hit_frac(const TridigitalScratch *a, const int *out, int clen) {
    if (!a->wordset) return 0.0;
    static _Thread_local char buf[MAX_CIPHER_LENGTH + 1];
    int den = 0, hits = 0, i = 0;
    while (i < clen) {
        if (out[i] < 0) { i++; continue; }
        int j = i; while (j < clen && out[j] >= 0) j++;
        int wlen = j - i;
        if (wlen >= TRI_MIN_WORD) {
            for (int k = 0; k < wlen; k++) buf[k] = index_to_char(out[i + k]);
            buf[wlen] = '\0';
            den++;
            if (tri_wordset_has(a->wordset, a->wordset_cap, buf)) hits++;
        }
        i = j;
    }
    return (den > 0) ? (double) hits / den : 0.0;
}

// Graded dictionary coverage (greedy longest non-overlapping match within each known word run,
// min length TRI_MIN_WORD) -- the gaming-resistant CROSS-CONFIG selection metric.
static double tri_coverage(const TridigitalScratch *a, const int *out, int clen) {
    if (!a->wordset) return 0.0;
    static _Thread_local char buf[MAX_CIPHER_LENGTH + 1];
    int maxL = a->max_word_len; if (maxL > 24) maxL = 24;
    int total = 0, covered = 0, i = 0;
    while (i < clen) {
        if (out[i] < 0) { i++; continue; }
        int j = i; while (j < clen && out[j] >= 0) j++;
        int wlen = j - i;
        for (int k = 0; k < wlen; k++) buf[k] = index_to_char(out[i + k]);
        total += wlen;
        int p = 0;
        while (p < wlen) {
            int hit = 0, Lmax = (wlen - p < maxL) ? wlen - p : maxL;
            for (int L = Lmax; L >= TRI_MIN_WORD; L--) {
                char c = buf[p + L]; buf[p + L] = '\0';
                bool in = tri_wordset_has(a->wordset, a->wordset_cap, buf + p);
                buf[p + L] = c;
                if (in) { covered += L; p += L; hit = 1; break; }
            }
            if (!hit) p++;
        }
        i = j;
    }
    return (total > 0) ? (double) covered / total : 0.0;
}

// Score a key under a separator: n-gram of its decode + the word-hit reward (the mini-solve
// objective, matching what the engine ranks by via state_score + score_adjust).
static double tri_score_key(const TridigitalScratch *a, const int key[TRI_MAP_CELLS], int sep) {
    tri_decode(a, key, sep, a->beam, g_tri_dec);
    double s = ngram_score(g_tri_dec, a->clen, a->ngram, a->ngram_size);
    s += TRI_WORD_WEIGHT * tri_word_hit_frac(a, g_tri_dec, a->clen);
    return s;
}

// ===================================================================
//  Key moves (size-preserving: the partition stays 8x3 + 1x2)
// ===================================================================

static void tri_group_counts(const int key[TRI_MAP_CELLS], int cnt[TRI_NGROUPS]) {
    for (int g = 0; g < TRI_NGROUPS; g++) cnt[g] = 0;
    for (int L = 0; L < TRI_MAP_CELLS; L++) cnt[key[L]]++;
}

static void tri_random_key(int key[TRI_MAP_CELLS]) {
    int perm[TRI_MAP_CELLS];
    for (int i = 0; i < TRI_MAP_CELLS; i++) perm[i] = i;
    shuffle(perm, TRI_MAP_CELLS);
    int idx = 0;
    for (int g = 0; g < TRI_NGROUPS; g++) {
        int sz = (g < TRI_NGROUPS - 1) ? 3 : 2;              // 8 groups of 3 + 1 group of 2 = 26
        for (int s = 0; s < sz; s++) key[perm[idx++]] = g;
    }
}

static void tri_perturb_key(int key[TRI_MAP_CELLS]) {
    if (frand() < 0.70) {                                    // swap two letters in different groups
        int i, j, tries = 0;
        do { i = rand_int(0, TRI_MAP_CELLS); j = rand_int(0, TRI_MAP_CELLS); tries++; }
        while (key[i] == key[j] && tries < 20);
        int t = key[i]; key[i] = key[j]; key[j] = t;
    } else {                                                 // relocate into the 2-group (it wanders)
        int cnt[TRI_NGROUPS]; tri_group_counts(key, cnt);
        int g2 = -1; for (int g = 0; g < TRI_NGROUPS; g++) if (cnt[g] == 2) { g2 = g; break; }
        if (g2 < 0) { int i = rand_int(0, TRI_MAP_CELLS), j = rand_int(0, TRI_MAP_CELLS);
                      int t = key[i]; key[i] = key[j]; key[j] = t; return; }
        int L, tries = 0;
        do { L = rand_int(0, TRI_MAP_CELLS); tries++; } while (cnt[key[L]] != 3 && tries < 40);
        if (cnt[key[L]] == 3) key[L] = g2;
    }
}

// ===================================================================
//  Config pre-pass: separator-digit discrimination + mini-solve warm start
// ===================================================================

// Map-independent word-length fit of a candidate separator digit: split the stream on the
// digit, histogram the run lengths, and score the fit to english_word_length_frequencies[]
// (negative L1 distance -- higher is better). The true separator's runs are English words, so
// its length distribution matches; a wrong digit (too few/many, wrong lengths) fits poorly.
static double tri_sep_fit(const int *digits, int clen, int sep) {
    double hist[25]; for (int b = 0; b < 25; b++) hist[b] = 0.0;
    int nruns = 0, run = 0;
    for (int i = 0; i < clen; i++) {
        if (digits[i] == sep) { if (run > 0) { int b = (run > 25) ? 25 : run; hist[b - 1] += 1.0; nruns++; run = 0; } }
        else run++;
    }
    if (run > 0) { int b = (run > 25) ? 25 : run; hist[b - 1] += 1.0; nruns++; }
    if (nruns == 0) return -1e9;                             // separator never occurs -> useless
    double l1 = 0.0;
    for (int b = 0; b < 25; b++) { hist[b] /= nruns; l1 += fabs(hist[b] - english_word_length_frequencies[b]); }
    return -l1;
}

static double tri_mini_solve(const TridigitalScratch *a, int sep, int iters, int best_out[TRI_MAP_CELLS]) {
    int nrst = 3, per = iters / nrst;
    double best = -1e18;
    int key[TRI_MAP_CELLS], save[TRI_MAP_CELLS], mp[TRI_MAP_CELLS];
    for (int r = 0; r < nrst; r++) {
        tri_random_key(key);
        double cur = tri_score_key(a, key, sep);
        double rbest = cur; for (int i = 0; i < TRI_MAP_CELLS; i++) mp[i] = key[i];
        double ratio = (per > 1) ? pow(0.02 / 0.30, 1.0 / (per - 1)) : 1.0, temp = 0.30;
        for (int it = 0; it < per; it++, temp *= ratio) {
            for (int i = 0; i < TRI_MAP_CELLS; i++) save[i] = key[i];
            tri_perturb_key(key);
            double s = tri_score_key(a, key, sep);
            double dscore = s - cur;
            if (dscore >= 0 || frand() < exp(dscore / temp)) {
                cur = s;
                if (s > rbest) { rbest = s; for (int i = 0; i < TRI_MAP_CELLS; i++) mp[i] = key[i]; }
            } else {
                for (int i = 0; i < TRI_MAP_CELLS; i++) key[i] = save[i];
            }
        }
        if (rbest > best) { best = rbest; for (int i = 0; i < TRI_MAP_CELLS; i++) best_out[i] = mp[i]; }
    }
    return best;
}

// Rank the 10 candidate separator digits by word-length fit, keep the top `keep` (padded to
// TRI_MIN_CANDIDATES), then mini-solve each for a warm-start key.
static void tri_rank_configs(TridigitalScratch *a, int keep, bool verbose) {
    int    sep[TRI_N_CONFIGS];
    double fit[TRI_N_CONFIGS];
    for (int d = 0; d < 10; d++) { sep[d] = d; fit[d] = tri_sep_fit(a->digits, a->clen, d); }
    // Order all 10 by fit DESC (selection sort).
    for (int i = 0; i < TRI_N_CONFIGS; i++) {
        int bi = i;
        for (int j = i + 1; j < TRI_N_CONFIGS; j++) if (fit[j] > fit[bi]) bi = j;
        double tf = fit[i]; fit[i] = fit[bi]; fit[bi] = tf;
        int ts = sep[i]; sep[i] = sep[bi]; sep[bi] = ts;
    }
    int nkeep = keep;
    if (nkeep < TRI_MIN_CANDIDATES) nkeep = TRI_MIN_CANDIDATES;
    if (nkeep > TRI_N_CONFIGS) nkeep = TRI_N_CONFIGS;

    a->n_configs = nkeep;
    for (int i = 0; i < nkeep; i++) {
        a->cfg_sep[i] = sep[i];
        tri_mini_solve(a, sep[i], TRI_PREPASS_ITERS, a->warm[i]);
    }
    if (verbose) {
        printf("\ntridigital pre-pass: separator digits ranked by word-length fit, keep %d:\n", nkeep);
        for (int i = 0; i < nkeep; i++) printf("  sep=%d  wordlen-fit=%.4f\n", sep[i], fit[i]);
    }
}

// ===================================================================
//  CipherModel hooks
// ===================================================================

static int tri_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const TridigitalScratch *a = (const TridigitalScratch *) ctx->model_scratch;
    if (cap < 1 || a->active_config < 0 || a->active_config >= a->n_configs) return 0;
    int i = a->active_config;
    out[0].period = i;                                       // config index (warm-start lookup)
    out[0].j = 0; out[0].k = 0;
    out[0].aux[0] = a->cfg_sep[i];                           // the separator digit
    out[0].aux[1] = 0;
    return 1;
}

static void tri_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const TridigitalScratch *a = (const TridigitalScratch *) ctx->model_scratch;
    // Seed from the config's mini-solve warm map, jittered by a few moves for restart diversity.
    // (Fresh-random restarts were tried and rejected: the partition landscape is too rugged for a
    // random start to converge within the climb budget, so they only diluted the warm restarts.)
    int idx = cc->period;
    for (int i = 0; i < TRI_MAP_CELLS; i++) st->key[i] = a->warm[idx][i];
    int jit = rand_int(0, 6);
    for (int s = 0; s < jit; s++) tri_perturb_key(st->key);
    st->key_len = TRI_MAP_CELLS;
}

static void tri_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                        SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    tri_perturb_key(st->key);
    if (frand() < 0.10) tri_perturb_key(st->key);           // occasional coarser move
}

static void tri_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < TRI_MAP_CELLS; i++) dst->key[i] = src->key[i];
    dst->key_len = src->key_len;
}

static void tri_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                             SolverState *st, int *out, double *score_adjust) {
    const TridigitalScratch *a = (const TridigitalScratch *) ctx->model_scratch;
    tri_decode(a, st->key, cc->aux[0], a->beam, out);
    *score_adjust = TRI_WORD_WEIGHT * tri_word_hit_frac(a, out, a->clen);
}

// ===================================================================
//  Reporting
// ===================================================================

// Render the recovered grid (canonical: the 9 letter-columns in digit order, then the
// separator). A column's letters are printed in group order; the row within a column is not
// identifiable (it does not affect the cipher), so this is the equivalent free-row grid.
static void tri_render_grid(const int key[TRI_MAP_CELLS], int sep, char out[]) {
    int digit_of_group[TRI_NGROUPS];
    { int r = 0; for (int d = 0; d < 10; d++) if (d != sep) digit_of_group[r++] = d; }
    int o = 0;
    for (int d = 0; d < 10; d++) {
        if (d == sep) { out[o++] = '/'; continue; }          // separator column marker
        int g = -1; for (int gg = 0; gg < TRI_NGROUPS; gg++) if (digit_of_group[gg] == d) { g = gg; break; }
        int wrote = 0;
        for (int L = 0; L < TRI_MAP_CELLS; L++) if (key[L] == g) { out[o++] = index_to_char(L); wrote++; }
        while (wrote++ < TRI_GROUP_CAP) out[o++] = '.';       // pad short columns for alignment
    }
    out[o] = '\0';
}

static void tri_print_digits(const int digits[], int n) {
    for (int i = 0; i < n; i++) putchar('0' + digits[i]);
}

static void tri_report(const SolverCtx *ctx, const SolverConfig *cc,
                       const SolverState *st, double score, int *decrypted) {
    (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    TridigitalScratch *a = (TridigitalScratch *) ctx->model_scratch;
    int sep = cc->aux[0], C = a->clen;
    int nletters = tri_decode(a, st->key, sep, TRI_FINAL_BEAM, g_tri_dec);
    (void) nletters;

    // Always stash this config's best state so solve_tridigital can score its decode by
    // dictionary coverage and replay the cross-config winner. While `quiet`, that is the only output.
    a->stash_sep = sep; a->stash_score = score; a->stash_declen = C;
    for (int i = 0; i < TRI_MAP_CELLS; i++) a->stash_key[i] = st->key[i];
    for (int i = 0; i < C; i++) a->stash_decode[i] = g_tri_dec[i];
    if (a->quiet) return;

    int n_words_found = 0;
    static char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < C; i++) plaintext_string[i] = index_to_char(g_tri_dec[i]);
    plaintext_string[C] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char gridstr[64]; tri_render_grid(st->key, sep, gridstr);

    printf("\nResult Score: %.2f | Words: %d | separator=%d | grid=%s\n",
        score, n_words_found, sep, gridstr);

    tri_print_digits(a->digits, C);
    printf("\n");
    print_text(g_tri_dec, C);
    printf("\n%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(g_tri_dec, ctx->result->decrypted, C);
        ctx->result->decrypted_len = C;
    }

    // >>> score, [words,] type, grid=, sep=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, grid=%s, sep=%d, ", score, n_words_found, cfg->cipher_type, gridstr, sep);
    else
        printf(">>> %.2f, %d, grid=%s, sep=%d, ", score, cfg->cipher_type, gridstr, sep);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    tri_print_digits(a->digits, C);
    printf(", ");
    print_text(g_tri_dec, C);
    printf("\n");
}

static const CipherModel TRIDIGITAL_MODEL = {
    .name = "tridigital", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = tri_enumerate, .key_len = NULL,
    .seed = tri_seed, .perturb = tri_perturb, .copy_state = tri_copy,
    .decrypt = tri_decrypt_hook, .report = tri_report,
};

// ===================================================================
//  Entry point
// ===================================================================

static int tri_parse_digits(const char *s, int out[], int cap) {
    int n = 0;
    for (int i = 0; s[i] && n < cap; i++)
        if (s[i] >= '0' && s[i] <= '9') out[n++] = s[i] - '0';
    return n;
}

void solve_tridigital(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_len; (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != DEFAULT_ALPHABET_SIZE) {
        printf("\n\nERROR: Tridigital needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    int clen = tri_parse_digits(ciphertext_str, g_tri_digits, MAX_CIPHER_LENGTH);
    if (clen < 8) {
        printf("\n\nERROR: parsed only %d ciphertext digits; need >= 8. The Tridigital "
               "ciphertext must be a stream of digits 0-9.\n\n", clen);
        return;
    }
    // Tridigital's decode carries space sentinels; ensure ngram_score compacts them (Fact B).
    g_score_no_sentinel = false;

    static TridigitalScratch scratch;
    scratch.digits = g_tri_digits;
    scratch.clen = clen;
    scratch.ngram_size = cfg->ngram_size;
    scratch.ngram = shared->ngram_data;
    scratch.beam = TRI_BEAM;

    int keep = (cfg->n_primers > 0) ? cfg->n_primers : TRI_KEEP;
    tri_rank_configs(&scratch, keep, cfg->verbose);

    scratch.wordset = (cfg->dictionary_present && shared->dict)
                    ? tri_build_wordset(shared->dict, shared->n_dict_words, &scratch.wordset_cap)
                    : NULL;
    if (!scratch.wordset) scratch.wordset_cap = 0;
    scratch.max_word_len = shared->max_dict_word_len;

    if (cfg->verbose)
        printf("\ntridigital: %d ciphertext digits, %d config(s), selector=%s\n",
               clen, scratch.n_configs, scratch.wordset ? "dictionary-coverage" : "n-gram (no dictionary)");

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, clen, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    // Drive the engine once per kept config; choose the winner by the EXACT whole-word-hit
    // fraction of its decode (coverage then n-gram as tiebreaks). This is the load-bearing
    // gaming defence: with the TRUE separator the word runs are real word boundaries, so a
    // well-solved grid makes most runs whole dictionary words; a WRONG separator splits the
    // stream at arbitrary positions whose runs are almost never whole words -- so its gibberish
    // cannot win even though lenient substring coverage (or raw n-gram) can be gamed high.
    int win_sep = -1, win_key[TRI_MAP_CELLS], win_declen = 0;
    static int win_decode[MAX_CIPHER_LENGTH];
    double win_metric = -1e18, win_cov = -1.0, win_score = 0.0;
    for (int i = 0; i < scratch.n_configs; i++) {
        scratch.active_config = i;
        scratch.quiet = true;
        run_solver(&TRIDIGITAL_MODEL, &ctx);
        double whf = tri_word_hit_frac(&scratch, scratch.stash_decode, scratch.stash_declen);
        double cov = tri_coverage(&scratch, scratch.stash_decode, scratch.stash_declen);
        // Primary: whole-word-hit fraction (0 when no dictionary => n-gram fallback via tiebreaks).
        double metric = scratch.wordset ? whf : scratch.stash_score;
        bool better = (metric > win_metric) ||
                      (metric == win_metric && cov > win_cov) ||
                      (metric == win_metric && cov == win_cov && scratch.stash_score > win_score);
        if (better) {
            win_metric = metric; win_cov = cov; win_score = scratch.stash_score;
            win_sep = scratch.stash_sep; win_declen = scratch.stash_declen;
            for (int k = 0; k < TRI_MAP_CELLS; k++) win_key[k] = scratch.stash_key[k];
            for (int k = 0; k < win_declen; k++) win_decode[k] = scratch.stash_decode[k];
        }
        if (cfg->verbose)
            printf("  sep=%d: n-gram=%.4f wordhit=%.3f coverage=%.3f\n",
                   scratch.stash_sep, scratch.stash_score, whf, cov);
    }

    // Replay the winner through tri_report (non-quiet) for the canonical human + >>> output.
    if (win_sep >= 0) {
        (void) win_decode; (void) win_declen;
        SolverState ws; ws.key_len = TRI_MAP_CELLS;
        for (int k = 0; k < TRI_MAP_CELLS; k++) ws.key[k] = win_key[k];
        SolverConfig wc; wc.period = 0; wc.j = 0; wc.k = 0; wc.aux[0] = win_sep; wc.aux[1] = 0;
        scratch.quiet = false;
        tri_report(&ctx, &wc, &ws, win_score, NULL);
    }

    if (scratch.wordset) free(scratch.wordset);
}
