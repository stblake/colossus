#include "trifid_solver.h"
#include "engine.h"
#include "trans_common.h"

// =====================================================================
//  Trifid solver (TYPE trifid)
// =====================================================================
//
// Trifid (Delastelle) is Bifid lifted into three dimensions: each letter fractionates
// into three Polybius-CUBE coordinates (layer, row, column) instead of two, and the
// coordinate stream is re-grouped block-by-block (block size = the PERIOD) into triples.
// Breaking it is the same two coupled problems as Bifid -- recover the period and recover
// the keyed cube -- attacked the same way. The cube search is a hill-climb / anneal over
// the cube itself (a permutation of 0..n-1 carried in st->key) with n-gram scoring; the
// move set is a dominant single-cell swap plus the cube's structured moves (swap two
// whole planes along an axis, reflect along an axis), which jump the basins a single
// cell swap cannot escape. The period is recovered by the same columnar index-of-
// coincidence test as Bifid (trifid_estimate_periods): at the correct period every
// within-block position is a fixed coordinate class sharing one elevated distribution,
// so the columnar IoC stands out (and also peaks at multiples, so the true period is in
// the top-K but not always rank 1). We anneal the estimator's top-K periods (one engine
// config each) and the n-gram score picks the winner. -period pins a single period.
//
// No anti-collapse penalty is needed (a cube is a bijection), so score_adjust stays 0.

#define TRIFID_MAX_PERIODS 64

typedef struct {
    int side;                        // cube side (3 for the 27-letter default)
    int side2;                       // side*side
    int cube_size;                   // side^3 == g_alpha
    int n_periods;                   // number of candidate periods
    int periods[TRIFID_MAX_PERIODS]; // the candidate periods (config order)
} TrifidScratch;

// --- period estimation --------------------------------------------------------
//
// Rank trial periods in [min_p .. max_p] by the columnar Index of Coincidence
// (mean_ioc, perioc.c) -- the same statistic Bifid and the Vigenere period estimator
// use, which peaks at the trifid period because each within-block position becomes a
// coordinate-pure column. Returns the top n_want periods by raw IoC into out[]
// (descending).
int trifid_estimate_periods(int cipher[], int len, int min_p, int max_p,
                            int n_want, int out[], bool verbose) {
    static int col[MAX_CIPHER_LENGTH];
    static double ioc[MAX_CIPHER_LENGTH];     // ioc[p] for p in [min_p..max_p]
    if (max_p > len / 2) max_p = len / 2;
    if (max_p < min_p) max_p = min_p;
    if (n_want < 1) n_want = 1;
    if (n_want > TRIFID_MAX_PERIODS) n_want = TRIFID_MAX_PERIODS;

    for (int p = min_p; p <= max_p; p++)
        ioc[p] = mean_ioc(cipher, len, p, col);

    if (verbose) {
        printf("\nTrifid period estimate (columnar IoC):\n  period\tIoC\n");
        for (int p = min_p; p <= max_p; p++) printf("  %d\t%.4f\n", p, ioc[p]);
    }

    // Selection by repeated max (n_want is small); mark taken periods with -1.
    int n = 0;
    for (; n < n_want; n++) {
        int best_p = -1;
        double best_ioc = -1.0;
        for (int p = min_p; p <= max_p; p++)
            if (ioc[p] >= 0.0 && ioc[p] > best_ioc) { best_ioc = ioc[p]; best_p = p; }
        if (best_p < 0) break;
        out[n] = best_p;
        ioc[best_p] = -1.0;
    }
    if (verbose) {
        printf("  -> annealing periods:");
        for (int i = 0; i < n; i++) printf(" %d", out[i]);
        printf("\n");
    }
    return n;
}

// One config per candidate period; period carries the block size, key the cube.
static int trifid_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const TrifidScratch *t = (const TrifidScratch *) ctx->model_scratch;
    int n = t->n_periods;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = t->periods[i];
        out[i].j = 0; out[i].k = 0; out[i].aux[0] = 0; out[i].aux[1] = 0;
    }
    return n;
}

// Seed: a uniformly random cube (Fisher-Yates shuffle of 0..cube_size-1).
static void trifid_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    const TrifidScratch *t = (const TrifidScratch *) ctx->model_scratch;
    int n = t->cube_size;
    for (int i = 0; i < n; i++) st->key[i] = i;
    for (int i = n - 1; i > 0; i--) {
        int j = rand_int(0, i + 1);
        int tmp = st->key[i]; st->key[i] = st->key[j]; st->key[j] = tmp;
    }
    st->key_len = n;
}

// Swap the two planes a, b along coordinate axis `ax` (0 layer / 1 row / 2 column):
// every cell whose `ax` coordinate is a is exchanged with the matching b cell.
static void trifid_swap_plane(int key[], int side, int side2, int ax, int a, int b) {
    if (a == b) return;
    for (int u = 0; u < side; u++)
        for (int v = 0; v < side; v++) {
            int pa, pb;
            if (ax == 0)      { pa = a * side2 + u * side + v;  pb = b * side2 + u * side + v; }
            else if (ax == 1) { pa = u * side2 + a * side + v;  pb = u * side2 + b * side + v; }
            else              { pa = u * side2 + v * side + a;  pb = u * side2 + v * side + b; }
            int tmp = key[pa]; key[pa] = key[pb]; key[pb] = tmp;
        }
}

// Reflect the cube along coordinate axis `ax` (mirror coordinate w <-> side-1-w).
static void trifid_reflect(int key[], int side, int side2, int ax) {
    for (int w = 0; w < side / 2; w++)
        trifid_swap_plane(key, side, side2, ax, w, side - 1 - w);
}

// Neighbour move on the cube (the 3D analogue of the Bifid square move): 82% swap two
// cells; 12% swap two whole planes along a random axis (4% each axis); 6% reflect along
// a random axis (2% each axis). The structured moves jump the basins a single cell swap
// cannot escape, exactly as the Bifid row/column swaps and reflections do.
static void trifid_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                           SolverState *st, bool *force_primary) {
    const TrifidScratch *t = (const TrifidScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int s = t->side, s2 = t->side2, n = t->cube_size;
    double r = frand();
    if (r < 0.82) {                              // swap two cells
        int a = rand_int(0, n), c = rand_int(0, n);
        int tmp = st->key[a]; st->key[a] = st->key[c]; st->key[c] = tmp;
    } else if (r < 0.94) {                       // swap two planes along a random axis
        int ax = rand_int(0, 3);
        int a = rand_int(0, s), b = rand_int(0, s);
        trifid_swap_plane(st->key, s, s2, ax, a, b);
    } else {                                     // reflect along a random axis
        int ax = rand_int(0, 3);
        trifid_reflect(st->key, s, s2, ax);
    }
}

static void trifid_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    // cube_size cells; copy a fixed upper bound so the engine need not know the side.
    for (int i = 0; i < TRIFID_MAX_CELLS; i++) dst->key[i] = src->key[i];
}

static void trifid_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                SolverState *st, int *out, double *score_adjust) {
    const TrifidScratch *t = (const TrifidScratch *) ctx->model_scratch;
    trifid_decrypt(ctx->cipher, ctx->cipher_len, st->key, t->side, cc->period, out);
    *score_adjust = 0.0;
}

// Render the cube as `side` indented layers of letters.
static void trifid_print_cube(const int cube[], int side) {
    int side2 = side * side;
    for (int l = 0; l < side; l++) {
        printf("    layer %d:\n", l + 1);
        for (int r = 0; r < side; r++) {
            printf("      ");
            for (int c = 0; c < side; c++) printf("%c ", index_to_char(cube[l * side2 + r * side + c]));
            printf("\n");
        }
    }
}

static void trifid_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const TrifidScratch *t = (const TrifidScratch *) ctx->model_scratch;
    printf("\n  period %d, cube:\n", cc->period);
    trifid_print_cube(st->key, t->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "trifid");
}

static void trifid_report(const SolverCtx *ctx, const SolverConfig *cc,
                          const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const TrifidScratch *t = (const TrifidScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, side = t->side, n = t->cube_size;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The recovered cube, read cell-major (layer/row/column order).
    char cubestr[TRIFID_MAX_CELLS + 1];
    for (int i = 0; i < n; i++) cubestr[i] = index_to_char(st->key[i]);
    cubestr[n] = '\0';

    printf("\nResult Score: %.2f | Words: %d | period=%d | cube=%s\n",
        score, n_words_found, cc->period, cubestr);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_solution_check(decrypted, len);
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered %dx%dx%d cube (cell major):\n", side, side, side);
    trifid_print_cube(st->key, side);

    // Publish the recovered solution for callers that pass a SolveResult.
    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = cc->period;   // report the recovered period here
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, period=, cube=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, period=%d, cube=%s, ",
            score, n_words_found, cfg->cipher_type, cc->period, cubestr);
    else
        printf(">>> %.2f, %d, period=%d, cube=%s, ",
            score, cfg->cipher_type, cc->period, cubestr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel TRIFID_MODEL = {
    .name = "trifid", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = trifid_enumerate, .key_len = NULL,
    .seed = trifid_seed, .perturb = trifid_perturb, .copy_state = trifid_copy,
    .decrypt = trifid_decrypt_hook, .report = trifid_report,
    .report_verbose = trifid_report_verbose,
};

// Integer cube root: returns s if n == s^3 for some s >= 1, else -1.
static int exact_icbrt(int n) {
    for (int s = 1; s * s * s <= n; s++)
        if (s * s * s == n) return s;
    return -1;
}

void solve_trifid(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Trifid needs a perfect-cube (e.g. 3x3x3 -> 27) alphabet; the binary forces 27 for
    // the default, but guard in case solve_trifid is driven directly.
    int side = exact_icbrt(g_alpha);
    if (side < 2 || side > TRIFID_MAX_SIDE || side * side * side != g_alpha) {
        printf("\n\nERROR: Trifid needs a perfect-cube alphabet (e.g. 27; got %d). "
               "Use the default -type trifid (27 symbols: A..Z + '%c').\n\n",
               g_alpha, TRIFID_EXTRA_CHAR);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Trifid solve.\n\n");
        return;
    }
    // Trifid letters must be solid: reject any sentinel (space/punctuation).
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a symbol outside the 27-symbol cube "
                   "alphabet (A..Z + '%c') at position %d; remove stray characters.\n\n",
                   TRIFID_EXTRA_CHAR, i);
            return;
        }

    TrifidScratch scratch;
    scratch.side = side;
    scratch.side2 = side * side;
    scratch.cube_size = g_alpha;

    // Candidate periods: pinned, or the estimator's top-K over [2 .. max_period].
    if (cfg->period_present) {
        scratch.periods[0] = cfg->period;
        scratch.n_periods = 1;
        if (cfg->verbose) printf("\ntrifid: period pinned to %d\n", cfg->period);
    } else {
        int max_p = (cfg->max_period > 0) ? cfg->max_period : 20;
        if (max_p > cipher_len / 2) max_p = cipher_len / 2;
        int n_want = cfg->n_periods;
        scratch.n_periods = trifid_estimate_periods(cipher_indices, cipher_len, 2, max_p,
            n_want, scratch.periods, cfg->verbose);
        if (scratch.n_periods < 1) { scratch.periods[0] = 2; scratch.n_periods = 1; }
    }

    if (cfg->verbose)
        printf("\ntrifid: %d positions, %d-symbol alphabet %s, %d candidate period(s)\n",
            cipher_len, g_alpha, g_idx_to_char_arr, scratch.n_periods);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // trifid_report fills it (may be NULL for CLI use)

    run_solver(&TRIFID_MODEL, &ctx);
}
