#include "monome_dinome_solver.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Monome-Dinome solver -- see monome_dinome_solver.h
// =====================================================================

// Decoded-plaintext + tiled scratch (single-threaded), off the stack.
static _Thread_local int g_md_decode[MAX_CIPHER_LENGTH];
static _Thread_local int g_md_tiled[MAX_CIPHER_LENGTH];
// Parsed ciphertext digit stream, filled by solve_monome_dinome.
static int g_md_digits[MAX_CIPHER_LENGTH];

typedef struct {
    int *digits;                             // ciphertext digit stream (length clen)
    int  clen;                               // digit-stream length (== engine cipher_len)
    int  ngram_size;
    float *ngram;
    int  n_configs;                          // number of validity-kept indicator-pair configs
    int  cfg_r0[MD_N_CONFIGS];               // config indicator pairs (r0 < r1), validity-kept
    int  cfg_r1[MD_N_CONFIGS];
    int  warm[MD_N_CONFIGS][MD_MAP_CELLS];   // mini-solve warm-start map per kept config

    // Per-config engine driving + gaming-resistant cross-config selection (see the header):
    // the engine is run once per kept config; md_report stashes that config's best state here
    // (no printing while `quiet`), and solve_monome_dinome selects across configs by DICTIONARY
    // COVERAGE of the best decode -- a signal n-gram GAMING cannot fake (a mis-segmented free
    // substitution reaches a high n-gram score but has few real words).
    int   active_config;                     // the single config md_enumerate exposes to run_solver
    bool  quiet;                             // suppress md_report printing during the per-config pass
    int   stash_r0, stash_r1;                // md_report stashes the reported config + state here
    int   stash_map[MD_MAP_CELLS];
    double stash_score;
    int   stash_decode[MAX_CIPHER_LENGTH];
    int   stash_declen;
    // Dictionary word-set for coverage scoring (built once from shared->dict); NULL => none.
    char **wordset;                          // open-addressing hash table of dict words (len >= 4)
    int   wordset_cap;                       // table capacity (power of two), 0 if no dictionary
    int   max_word_len;
} MonomeDinomeScratch;

// ===================================================================
//  Gaming-resistant selector: dictionary coverage of a decode
// ===================================================================

static unsigned long md_str_hash(const char *s) {           // FNV-1a
    unsigned long h = 1469598103934665603UL;
    for (; *s; s++) { h ^= (unsigned char) *s; h *= 1099511628211UL; }
    return h;
}

// Build an open-addressing hash set of the dictionary words of length >= 4 (the discriminating
// lengths; 3-letter words are too easily hit by chance). Returns the table (caller frees) or
// NULL, and sets *cap. Words are stored as borrowed pointers into `dict` (no copy).
static char **md_build_wordset(char **dict, int n_dict, int *cap_out) {
    if (!dict || n_dict <= 0) { *cap_out = 0; return NULL; }
    int cap = 1; while (cap < 4 * n_dict) cap <<= 1;         // <=50% load
    char **tab = calloc((size_t) cap, sizeof(char *));
    if (!tab) { *cap_out = 0; return NULL; }
    for (int i = 0; i < n_dict; i++) {
        if (strlen(dict[i]) < 4) continue;
        unsigned long h = md_str_hash(dict[i]) & (unsigned long) (cap - 1);
        while (tab[h]) h = (h + 1) & (unsigned long) (cap - 1);
        tab[h] = dict[i];
    }
    *cap_out = cap;
    return tab;
}

static bool md_wordset_has(char **tab, int cap, const char *frag) {
    unsigned long h = md_str_hash(frag) & (unsigned long) (cap - 1);
    while (tab[h]) { if (strcmp(tab[h], frag) == 0) return true; h = (h + 1) & (unsigned long) (cap - 1); }
    return false;
}

// Fraction of the decode covered by dictionary words (greedy longest non-overlapping match,
// min length 4). Real plaintext covers ~0.6-0.8; an n-gram-gamed mis-segmentation ~0.2-0.3.
static double md_coverage(const MonomeDinomeScratch *a, const int *decode, int m) {
    if (!a->wordset || m <= 0) return 0.0;
    char buf[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < m; i++) buf[i] = index_to_char(decode[i]);
    buf[m] = '\0';
    int maxL = a->max_word_len; if (maxL > 24) maxL = 24;
    int covered = 0, i = 0;
    while (i < m) {
        int hit = 0;
        for (int L = (m - i < maxL) ? m - i : maxL; L >= 4; L--) {
            char c = buf[i + L]; buf[i + L] = '\0';
            bool in = md_wordset_has(a->wordset, a->wordset_cap, buf + i);
            buf[i + L] = c;
            if (in) { covered += L; i += L; hit = 1; break; }
        }
        if (!hit) i++;
    }
    return (double) covered / m;
}

// ===================================================================
//  Code indexing + map decode (the solver's board-free decrypt)
// ===================================================================
//
// Within a config {r0,r1} the 24 token codes are indexed: the 8 non-indicator monome digits
// -> 0..7 (via single_slot[]); the (r0,d2) dinomes -> 8+single_slot[d2]; the (r1,d2) dinomes
// -> 16+single_slot[d2]. map[code] in {0..23}. Every valid token decodes to a letter -- there
// is no figure-shift and no null cell (the Monome-Dinome box is a clean 24-letter fill). This
// matches the primitive's board decrypt exactly (map[code(token)] = board.letter_at_cell[cell]).

static void md_single_slot(int r0, int r1, int single_slot[10]) {
    int s = 0;
    for (int d = 0; d < 10; d++) single_slot[d] = (d == r0 || d == r1) ? -1 : s++;
}

static int md_decode_map(const int *digits, int clen, int r0, int r1, const int single_slot[10],
                         const int map[MD_MAP_CELLS], int out[], int filler,
                         int *n_tokens, int *n_valid) {
    int o = 0, nt = 0, nv = 0, i = 0;
    while (i < clen) {
        int g = digits[i], code;
        if (g < 0 || g > 9) { i++; continue; }
        if (g == r0 || g == r1) {                                            // a dinome
            if (i + 1 >= clen) { out[o++] = filler; nt++; i++; continue; }   // truncated
            int g2 = digits[i + 1];
            i += 2;
            if (g2 < 0 || g2 > 9 || g2 == r0 || g2 == r1) {                  // 2nd digit is an
                out[o++] = filler; nt++; continue;                          //   indicator -> invalid
            }
            code = ((g == r0) ? 8 : 16) + single_slot[g2];
        } else {                                                             // a monome
            i += 1;
            code = single_slot[g];                                           // 0..7
        }
        nt++; nv++;
        out[o++] = map[code];                                                // a letter 0..23
    }
    if (n_tokens) *n_tokens = nt;
    if (n_valid)  *n_valid  = nv;
    return o;
}

// Score a decoded map: mean n-gram of the decode tiled to `clen` (length-fair) + validity.
static double md_score_map(const MonomeDinomeScratch *a, int r0, int r1, const int single_slot[10],
                           const int map[MD_MAP_CELLS]) {
    int nt = 0, nv = 0;
    int m = md_decode_map(a->digits, a->clen, r0, r1, single_slot, map,
                          g_md_decode, MD_FILLER, &nt, &nv);
    if (m <= 0) return -1e9;
    for (int i = 0; i < a->clen; i++) g_md_tiled[i] = g_md_decode[i % m];
    double s = ngram_score(g_md_tiled, a->clen, a->ngram, a->ngram_size);
    s += MD_VALID_WEIGHT * ((nt > 0) ? (double) nv / nt : 0.0);
    return s;
}

// ===================================================================
//  Config pre-pass: a short substitution hill-climb per indicator pair (the real ranking)
// ===================================================================

// One simulated-annealing swap climb of the code->letter map from a random start; returns the
// best score reached and its map. Metropolis accept-worse (geometric T0 -> Tmin) escapes the
// gibberish local optima a greedy climb sticks in, so a WRONG config cannot masquerade as the
// true one and the true config reliably ranks near-optimal.
static double md_anneal_map(const MonomeDinomeScratch *a, int r0, int r1, const int single_slot[10],
                            int iters, double t0, double tmin, int best_map_out[MD_MAP_CELLS]) {
    int map[MD_MAP_CELLS];
    for (int i = 0; i < MD_MAP_CELLS; i++) map[i] = i;   // identity permutation of 0..23
    shuffle(map, MD_MAP_CELLS);
    double cur = md_score_map(a, r0, r1, single_slot, map);
    double best = cur;
    for (int i = 0; i < MD_MAP_CELLS; i++) best_map_out[i] = map[i];

    double ratio = (iters > 1) ? pow(tmin / t0, 1.0 / (iters - 1)) : 1.0, temp = t0;
    for (int it = 0; it < iters; it++, temp *= ratio) {
        int i = rand_int(0, MD_MAP_CELLS), j = rand_int(0, MD_MAP_CELLS);
        if (i == j) continue;
        int t = map[i]; map[i] = map[j]; map[j] = t;
        double s = md_score_map(a, r0, r1, single_slot, map);
        double d = s - cur;
        if (d >= 0 || frand() < exp(d / temp)) {
            cur = s;
            if (s > best) { best = s; for (int k = 0; k < MD_MAP_CELLS; k++) best_map_out[k] = map[k]; }
        } else {
            t = map[i]; map[i] = map[j]; map[j] = t;                          // revert
        }
    }
    return best;
}

// Rank-and-warm-start a config: a few independent SA restarts, keep the best map + score.
static double md_mini_solve(const MonomeDinomeScratch *a, int r0, int r1, int iters,
                            int best_map_out[MD_MAP_CELLS]) {
    int single_slot[10];
    md_single_slot(r0, r1, single_slot);
    int nrst = 3, per = iters / nrst;
    double best = -1e18;
    int mp[MD_MAP_CELLS];
    for (int r = 0; r < nrst; r++) {
        double s = md_anneal_map(a, r0, r1, single_slot, per, 0.30, 0.02, mp);
        if (s > best) { best = s; for (int k = 0; k < MD_MAP_CELLS; k++) best_map_out[k] = mp[k]; }
    }
    return best;
}

// Structural token validity of a config {a,b} -- MAP-INDEPENDENT and purely structural: the
// fraction of tokens legal under indicator pair {a,b} (a dinome whose second digit is itself
// an indicator, or a truncated final dinome, is illegal). The TRUE pair is always exactly 1.0
// (the encryptor emits only legal tokens); a wrong pair drops below 1.0 as soon as a dinome's
// second digit lands on an indicator. So validity is a cheap, deterministic pre-filter that is
// GUARANTEED to retain the true config -- the "decoupling reward" pattern (a statistic that
// depends only on the segmentation half, giving it a gradient flat in the substitution half).
static double md_config_validity(const int *digits, int clen, int a, int b) {
    int i = 0, nt = 0, nv = 0;
    while (i < clen) {
        int g = digits[i];
        if (g < 0 || g > 9) { i++; continue; }
        if (g == a || g == b) {
            if (i + 1 >= clen) { nt++; i++; continue; }                 // truncated dinome
            int g2 = digits[i + 1]; i += 2;
            if (g2 < 0 || g2 > 9 || g2 == a || g2 == b) { nt++; continue; }   // 2nd digit an indicator
            nt++; nv++;
        } else { i++; nt++; nv++; }                                      // monome (always legal)
    }
    return (nt > 0) ? (double) nv / nt : 0.0;
}

// Pre-filter the 45 indicator pairs by structural validity, then mini-solve the kept set for
// warm-start maps. The kept set is the FULLY-VALID configs (validity within MD_VALID_EPS of the
// max -- always includes the true pair), padded to MD_MIN_CANDIDATES by the next-most-valid
// configs (insurance for a slightly indel-corrupted ciphertext whose true validity < 1.0) and
// capped at `keep`. The mini-solve only warm-starts (the engine anneal does the real solving),
// so its n-gram ranking merely orders the kept configs -- the validity filter, not the n-gram
// score, is what guarantees the true pair survives.
static void md_rank_configs(MonomeDinomeScratch *a, int keep, bool verbose) {
    int r0s[MD_N_CONFIGS], r1s[MD_N_CONFIGS];
    double val[MD_N_CONFIGS];
    int n = 0;
    double maxv = 0.0;
    for (int r0 = 0; r0 < 10; r0++)
        for (int r1 = r0 + 1; r1 < 10; r1++) {
            r0s[n] = r0; r1s[n] = r1;
            val[n] = md_config_validity(a->digits, a->clen, r0, r1);
            if (val[n] > maxv) maxv = val[n];
            n++;
        }
    // Order all configs by validity DESC (selection sort, carrying r0/r1).
    for (int i = 0; i < n; i++) {
        int best = i;
        for (int j = i + 1; j < n; j++) if (val[j] > val[best]) best = j;
        double tv = val[i]; val[i] = val[best]; val[best] = tv;
        int t0 = r0s[i]; r0s[i] = r0s[best]; r0s[best] = t0;
        int t1 = r1s[i]; r1s[i] = r1s[best]; r1s[best] = t1;
    }
    // Candidate set: the fully-valid prefix, padded to MD_MIN_CANDIDATES, capped at `keep`.
    int ncand = 0;
    while (ncand < n && val[ncand] >= maxv - MD_VALID_EPS) ncand++;
    if (ncand < MD_MIN_CANDIDATES) ncand = (n < MD_MIN_CANDIDATES) ? n : MD_MIN_CANDIDATES;
    if (ncand > keep) ncand = keep;

    // Mini-solve each kept config for a warm map + an n-gram score, then order by that score.
    static int wm[MD_N_CONFIGS][MD_MAP_CELLS];
    double sc[MD_N_CONFIGS];
    for (int i = 0; i < ncand; i++)
        sc[i] = md_mini_solve(a, r0s[i], r1s[i], MD_PREPASS_ITERS, wm[i]);
    for (int i = 0; i < ncand; i++) {
        int best = i;
        for (int j = i + 1; j < ncand; j++) if (sc[j] > sc[best]) best = j;
        double ts = sc[i]; sc[i] = sc[best]; sc[best] = ts;
        double tv = val[i]; val[i] = val[best]; val[best] = tv;
        int t0 = r0s[i]; r0s[i] = r0s[best]; r0s[best] = t0;
        int t1 = r1s[i]; r1s[i] = r1s[best]; r1s[best] = t1;
        for (int k = 0; k < MD_MAP_CELLS; k++) { int t = wm[i][k]; wm[i][k] = wm[best][k]; wm[best][k] = t; }
    }
    a->n_configs = ncand;
    for (int i = 0; i < ncand; i++) {
        a->cfg_r0[i] = r0s[i]; a->cfg_r1[i] = r1s[i];
        for (int k = 0; k < MD_MAP_CELLS; k++) a->warm[i][k] = wm[i][k];
    }
    if (verbose) {
        printf("\nmonome-dinome pre-pass: %d config(s) kept of 45 by validity (max %.4f):\n", ncand, maxv);
        for (int i = 0; i < ncand && i < 10; i++)
            printf("  {%d,%d} valid=%.4f nscore=%.4f\n", r0s[i], r1s[i], val[i], sc[i]);
    }
}

// ===================================================================
//  CipherModel hooks
// ===================================================================

// Expose the SINGLE active config to run_solver: solve_monome_dinome drives the engine once
// per kept config so the cross-config winner can be chosen by dictionary coverage (not by the
// gameable n-gram score run_solver would compare on).
static int md_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const MonomeDinomeScratch *a = (const MonomeDinomeScratch *) ctx->model_scratch;
    if (cap < 1 || a->active_config < 0 || a->active_config >= a->n_configs) return 0;
    int i = a->active_config;
    out[0].period = i;                              // config index (for the warm-start lookup)
    out[0].j = 0; out[0].k = 0;
    out[0].aux[0] = a->cfg_r0[i];
    out[0].aux[1] = a->cfg_r1[i];
    return 1;
}

// Seed from the config's mini-solve warm map, jittered by a few random swaps for diversity.
static void md_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const MonomeDinomeScratch *a = (const MonomeDinomeScratch *) ctx->model_scratch;
    int idx = cc->period;
    for (int i = 0; i < MD_MAP_CELLS; i++) st->key[i] = a->warm[idx][i];
    int jit = rand_int(0, 7);
    for (int s = 0; s < jit; s++) {
        int i = rand_int(0, MD_MAP_CELLS), j = rand_int(0, MD_MAP_CELLS);
        int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
    }
    st->key_len = MD_MAP_CELLS;
}

static void md_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                       SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int i = rand_int(0, MD_MAP_CELLS), j = rand_int(0, MD_MAP_CELLS);
    int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
    if (frand() < 0.10) {                            // occasional second swap (coarser move)
        int p = rand_int(0, MD_MAP_CELLS), q = rand_int(0, MD_MAP_CELLS);
        t = st->key[p]; st->key[p] = st->key[q]; st->key[q] = t;
    }
}

static void md_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < MD_MAP_CELLS; i++) dst->key[i] = src->key[i];
    dst->key_len = src->key_len;
}

static void md_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                            SolverState *st, int *out, double *score_adjust) {
    const MonomeDinomeScratch *a = (const MonomeDinomeScratch *) ctx->model_scratch;
    int single_slot[10];
    md_single_slot(cc->aux[0], cc->aux[1], single_slot);
    int C = a->clen, nt = 0, nv = 0;
    int m = md_decode_map(a->digits, C, cc->aux[0], cc->aux[1], single_slot, st->key,
                          g_md_decode, MD_FILLER, &nt, &nv);
    if (m <= 0) {
        for (int i = 0; i < C; i++) out[i] = MD_FILLER;
        *score_adjust = 0.0;
        return;
    }
    for (int i = 0; i < C; i++) out[i] = g_md_decode[i % m];
    *score_adjust = MD_VALID_WEIGHT * ((nt > 0) ? (double) nv / nt : 0.0);
}

// ===================================================================
//  Reporting
// ===================================================================

// Render the recovered box (canonical 0..9 labels): top row (8 non-indicator monomes),
// then the r0 and r1 rows. Keyed labels are not separately identifiable, so this is the
// equivalent fixed-label box.
static void md_render_board(int r0, int r1, const int map[MD_MAP_CELLS], char out[]) {
    int single_slot[10];
    md_single_slot(r0, r1, single_slot);
    int o = 0;
    for (int d = 0; d < 10; d++) if (single_slot[d] >= 0) out[o++] = index_to_char(map[single_slot[d]]);
    for (int d = 0; d < 8; d++) out[o++] = index_to_char(map[8 + d]);
    for (int d = 0; d < 8; d++) out[o++] = index_to_char(map[16 + d]);
    out[o] = '\0';
}

static void md_print_digits(const int digits[], int n) {
    for (int i = 0; i < n; i++) putchar('0' + digits[i]);
}

static void md_report(const SolverCtx *ctx, const SolverConfig *cc,
                      const SolverState *st, double score, int *decrypted) {
    (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    MonomeDinomeScratch *a = (MonomeDinomeScratch *) ctx->model_scratch;
    int r0 = cc->aux[0], r1 = cc->aux[1], C = a->clen, nt = 0, nv = 0;
    int single_slot[10];
    md_single_slot(r0, r1, single_slot);
    int m = md_decode_map(a->digits, C, r0, r1, single_slot, st->key,
                          g_md_decode, MD_FILLER, &nt, &nv);
    if (m > MAX_CIPHER_LENGTH) m = MAX_CIPHER_LENGTH;

    // Always stash this config's best state so solve_monome_dinome can score its decode by
    // dictionary coverage and replay the cross-config winner. During the per-config pass
    // (`quiet`) that stash is the ONLY output -- no printing, no result population.
    a->stash_r0 = r0; a->stash_r1 = r1; a->stash_score = score; a->stash_declen = m;
    for (int i = 0; i < MD_MAP_CELLS; i++) a->stash_map[i] = st->key[i];
    for (int i = 0; i < m; i++) a->stash_decode[i] = g_md_decode[i];
    if (a->quiet) return;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < m; i++) plaintext_string[i] = index_to_char(g_md_decode[i]);
    plaintext_string[m] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char boardstr[MD_MAP_CELLS + 1]; md_render_board(r0, r1, st->key, boardstr);

    printf("\nResult Score: %.2f | Words: %d | valid=%d/%d | indicators=%d%d | board=%s\n",
        score, n_words_found, nv, nt, r0, r1, boardstr);

    md_print_digits(a->digits, C);
    printf("\n");
    print_text(g_md_decode, m);
    printf("\n%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(g_md_decode, ctx->result->decrypted, m);
        ctx->result->decrypted_len = m;
    }

    // >>> score, [words,] type, board=, ind=, valid=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, board=%s, ind=%d%d, valid=%d/%d, ",
            score, n_words_found, cfg->cipher_type, boardstr, r0, r1, nv, nt);
    else
        printf(">>> %.2f, %d, board=%s, ind=%d%d, valid=%d/%d, ",
            score, cfg->cipher_type, boardstr, r0, r1, nv, nt);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    md_print_digits(a->digits, C);
    printf(", ");
    print_text(g_md_decode, m);
    printf("\n");
}

static const CipherModel MONOME_DINOME_MODEL = {
    .name = "monome-dinome", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = md_enumerate, .key_len = NULL,
    .seed = md_seed, .perturb = md_perturb, .copy_state = md_copy,
    .decrypt = md_decrypt_hook, .report = md_report,
};

// ===================================================================
//  Entry point
// ===================================================================

static int md_parse_digits(const char *s, int out[], int cap) {
    int n = 0;
    for (int i = 0; s[i] && n < cap; i++)
        if (s[i] >= '0' && s[i] <= '9') out[n++] = s[i] - '0';
    return n;
}

void solve_monome_dinome(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_len; (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != MD_NALPHA) {
        printf("\n\nERROR: Monome-Dinome needs the 24-letter alphabet (A..Z with J->I, Z->Y, "
               "got %d). Run -type md so the alphabet is forced.\n\n", g_alpha);
        return;
    }

    int clen = md_parse_digits(ciphertext_str, g_md_digits, MAX_CIPHER_LENGTH);
    if (clen < 4) {
        printf("\n\nERROR: parsed only %d ciphertext digits; need >= 4. The Monome-Dinome "
               "ciphertext must be a stream of digits 0-9.\n\n", clen);
        return;
    }

    static MonomeDinomeScratch scratch;
    scratch.digits = g_md_digits;
    scratch.clen = clen;
    scratch.ngram_size = cfg->ngram_size;
    scratch.ngram = shared->ngram_data;

    int keep = (cfg->n_primers > 0) ? cfg->n_primers : MD_KEEP;
    md_rank_configs(&scratch, keep, cfg->verbose);

    // Dictionary word-set for the gaming-resistant cross-config selector (built once).
    scratch.wordset = (cfg->dictionary_present && shared->dict)
                    ? md_build_wordset(shared->dict, shared->n_dict_words, &scratch.wordset_cap)
                    : NULL;
    if (!scratch.wordset) scratch.wordset_cap = 0;
    scratch.max_word_len = shared->max_dict_word_len;

    if (cfg->verbose)
        printf("\nmonome-dinome: %d ciphertext digits, %d-letter alphabet %s, %d config(s), "
               "selector=%s\n", clen, g_alpha, g_idx_to_char_arr, scratch.n_configs,
               scratch.wordset ? "dictionary-coverage" : "n-gram (no dictionary)");

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, clen, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    // Drive the engine ONCE PER kept config; choose the winner by dictionary coverage of its
    // best decode (falling back to the n-gram score when no dictionary is loaded). The n-gram
    // score alone is gameable across configs -- a mis-segmented free substitution can out-score
    // the true segmentation -- but dictionary coverage cannot be faked (see the header).
    int  win_r0 = -1, win_r1 = -1, win_map[MD_MAP_CELLS], win_declen = 0;
    static int win_decode[MAX_CIPHER_LENGTH];
    double win_metric = -1e18, win_score = 0.0;
    for (int i = 0; i < scratch.n_configs; i++) {
        scratch.active_config = i;
        scratch.quiet = true;
        run_solver(&MONOME_DINOME_MODEL, &ctx);
        double cov = md_coverage(&scratch, scratch.stash_decode, scratch.stash_declen);
        // Primary key: coverage (0 when no dictionary => pure n-gram tiebreak); tiebreak: n-gram.
        double metric = scratch.wordset ? cov : scratch.stash_score;
        double tie    = scratch.stash_score;
        if (metric > win_metric || (metric == win_metric && tie > win_score)) {
            win_metric = metric; win_score = scratch.stash_score;
            win_r0 = scratch.stash_r0; win_r1 = scratch.stash_r1; win_declen = scratch.stash_declen;
            for (int k = 0; k < MD_MAP_CELLS; k++) win_map[k] = scratch.stash_map[k];
            for (int k = 0; k < win_declen; k++) win_decode[k] = scratch.stash_decode[k];
        }
        if (cfg->verbose)
            printf("  config {%d,%d}: n-gram=%.4f coverage=%.3f\n",
                   scratch.stash_r0, scratch.stash_r1, scratch.stash_score, cov);
    }

    // Replay the winner through md_report (non-quiet) for the canonical human + >>> output.
    if (win_r0 >= 0) {
        SolverState ws; ws.key_len = MD_MAP_CELLS;
        for (int k = 0; k < MD_MAP_CELLS; k++) ws.key[k] = win_map[k];
        SolverConfig wc; wc.period = 0; wc.j = 0; wc.k = 0; wc.aux[0] = win_r0; wc.aux[1] = win_r1;
        scratch.quiet = false;
        md_report(&ctx, &wc, &ws, win_score, win_decode);
    }

    if (scratch.wordset) free(scratch.wordset);
}
