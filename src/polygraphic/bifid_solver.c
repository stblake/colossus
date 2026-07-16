#include "bifid_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Bifid solver (TYPE bifid)
// =====================================================================
//
// Bifid (Delastelle) fractionates each letter into Polybius-square coordinates, then
// re-pairs the coordinate stream block-by-block (block size = the PERIOD). Breaking it
// is two coupled problems: recover the period and recover the keyed square. The square
// search is identical to Playfair's -- the state is the square itself (a permutation of
// 0..n-1 carried in st->key), hill-climbed / annealed with n-gram scoring, the move set
// a dominant cell swap plus row/column swaps and reflections. The period is recovered by
// an index-of-coincidence test (bifid_estimate_periods): for the correct period every
// within-block position is a fixed coordinate-class (row-row or col-col letters) sharing
// one elevated distribution, so the columnar IoC at the true period stands out. We anneal
// the estimator's top-K periods (one engine config each) and the n-gram score picks the
// winner -- a wrong period yields gibberish and loses. -period pins a single period.
//
// No anti-collapse penalty is needed (a square is a bijection), so score_adjust stays 0.

#define BIFID_MAX_PERIODS 64

typedef struct {
    int side;                       // grid side (5 for the 25-letter default, up to 6)
    int grid_size;                  // side*side == g_alpha
    int n_periods;                  // number of candidate periods
    int periods[BIFID_MAX_PERIODS]; // the candidate periods (config order)
} BifidScratch;

// --- period estimation --------------------------------------------------------
//
// Rank trial periods in [min_p .. max_p] by the columnar Index of Coincidence
// (mean_ioc, perioc.c) -- the same statistic the Vigenere period estimator uses, which
// peaks at the bifid period because each within-block position becomes a coordinate-
// pure column. Returns the top n_want periods by raw IoC into out[] (descending).
int bifid_estimate_periods(int cipher[], int len, int min_p, int max_p,
                           int n_want, int out[], bool verbose) {
    static int col[MAX_CIPHER_LENGTH];
    static double ioc[MAX_CIPHER_LENGTH];     // ioc[p] for p in [min_p..max_p]
    if (max_p > len / 2) max_p = len / 2;
    if (max_p < min_p) max_p = min_p;
    if (n_want < 1) n_want = 1;
    if (n_want > BIFID_MAX_PERIODS) n_want = BIFID_MAX_PERIODS;

    for (int p = min_p; p <= max_p; p++)
        ioc[p] = mean_ioc(cipher, len, p, col);

    if (verbose) {
        printf("\nBifid period estimate (columnar IoC):\n  period\tIoC\n");
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

// One config per candidate period; period carries the block size, key the square.
static int bifid_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const BifidScratch *b = (const BifidScratch *) ctx->model_scratch;
    int n = b->n_periods;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = b->periods[i];
        out[i].j = 0; out[i].k = 0; out[i].aux[0] = 0; out[i].aux[1] = 0;
    }
    return n;
}

// Seed: a uniformly random square (Fisher-Yates shuffle of 0..grid_size-1).
static void bifid_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    const BifidScratch *b = (const BifidScratch *) ctx->model_scratch;
    int n = b->grid_size;
    for (int i = 0; i < n; i++) st->key[i] = i;
    for (int i = n - 1; i > 0; i--) {
        int j = rand_int(0, i + 1);
        int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
    }
    st->key_len = n;
}

// Neighbour move on the square (identical to Playfair's): 80% swap two cells; 8% swap
// two rows; 8% swap two columns; 2% reverse the whole grid; 1% flip rows; 1% flip
// columns. The larger moves jump the basins a single cell swap cannot escape.
static void bifid_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                          SolverState *st, bool *force_primary) {
    const BifidScratch *b = (const BifidScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int s = b->side, n = b->grid_size;
    double r = frand();
    if (r < 0.80) {                              // swap two cells
        int a = rand_int(0, n), c = rand_int(0, n);
        int t = st->key[a]; st->key[a] = st->key[c]; st->key[c] = t;
    } else if (r < 0.88) {                       // swap two rows
        int r1 = rand_int(0, s), r2 = rand_int(0, s);
        for (int c = 0; c < s; c++) {
            int t = st->key[r1 * s + c]; st->key[r1 * s + c] = st->key[r2 * s + c]; st->key[r2 * s + c] = t;
        }
    } else if (r < 0.96) {                       // swap two columns
        int c1 = rand_int(0, s), c2 = rand_int(0, s);
        for (int rr = 0; rr < s; rr++) {
            int t = st->key[rr * s + c1]; st->key[rr * s + c1] = st->key[rr * s + c2]; st->key[rr * s + c2] = t;
        }
    } else if (r < 0.98) {                       // reverse the whole grid (rotate 180)
        for (int i = 0, j = n - 1; i < j; i++, j--) {
            int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
        }
    } else if (r < 0.99) {                       // flip rows top<->bottom
        for (int r1 = 0, r2 = s - 1; r1 < r2; r1++, r2--)
            for (int c = 0; c < s; c++) {
                int t = st->key[r1 * s + c]; st->key[r1 * s + c] = st->key[r2 * s + c]; st->key[r2 * s + c] = t;
            }
    } else {                                     // flip columns left<->right
        for (int c1 = 0, c2 = s - 1; c1 < c2; c1++, c2--)
            for (int rr = 0; rr < s; rr++) {
                int t = st->key[rr * s + c1]; st->key[rr * s + c1] = st->key[rr * s + c2]; st->key[rr * s + c2] = t;
            }
    }
}

static void bifid_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    // grid_size cells; copy a fixed upper bound so the engine need not know the side.
    for (int i = 0; i < BIFID_MAX_GRID; i++) dst->key[i] = src->key[i];
}

static void bifid_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                               SolverState *st, int *out, double *score_adjust) {
    const BifidScratch *b = (const BifidScratch *) ctx->model_scratch;
    bifid_decrypt(ctx->cipher, ctx->cipher_len, st->key, b->side, cc->period, out);
    *score_adjust = 0.0;
}

// Render the square as an indented box of letters.
static void bifid_print_grid(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void bifid_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const BifidScratch *b = (const BifidScratch *) ctx->model_scratch;
    printf("\n  period %d, square:\n", cc->period);
    bifid_print_grid(st->key, b->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "bifid");
}

static void bifid_report(const SolverCtx *ctx, const SolverConfig *cc,
                         const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const BifidScratch *b = (const BifidScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, side = b->side, n = b->grid_size;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The recovered square, read row-major (unique only up to cyclic row/column
    // rotation, which all decrypt identically -- this is one representative).
    char gridstr[BIFID_MAX_GRID + 1];
    for (int i = 0; i < n; i++) gridstr[i] = index_to_char(st->key[i]);
    gridstr[n] = '\0';

    printf("\nResult Score: %.2f | Words: %d | period=%d | square=%s\n",
        score, n_words_found, cc->period, gridstr);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_solution_check(decrypted, len);
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered %dx%d square (row major):\n", side, side);
    bifid_print_grid(st->key, side);

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

    // One-liner summary: >>> score, [words,] type, period=, square=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, period=%d, square=%s, ",
            score, n_words_found, cfg->cipher_type, cc->period, gridstr);
    else
        printf(">>> %.2f, %d, period=%d, square=%s, ",
            score, cfg->cipher_type, cc->period, gridstr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel BIFID_MODEL = {
    .name = "bifid", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = bifid_enumerate, .key_len = NULL,
    .seed = bifid_seed, .perturb = bifid_perturb, .copy_state = bifid_copy,
    .decrypt = bifid_decrypt_hook, .report = bifid_report,
    .report_verbose = bifid_report_verbose,
};

void solve_bifid(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Bifid needs a perfect-square (e.g. 5x5 -> 25, 6x6 -> 36) alphabet; the binary
    // forces 25 for the default, but guard in case solve_bifid is driven directly.
    int side = exact_isqrt(g_alpha);
    if (side < 2 || side > BIFID_MAX_SIDE || side * side != g_alpha) {
        printf("\n\nERROR: Bifid needs a perfect-square alphabet (e.g. 25 or 36; got %d). "
               "Exclude a letter so it is 25 (e.g. -excludeletter J).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Bifid solve.\n\n");
        return;
    }
    // Bifid letters must be solid: reject any sentinel (space/punctuation).
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Bifid ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    BifidScratch scratch;
    scratch.side = side;
    scratch.grid_size = g_alpha;

    // Candidate periods: pinned, or the estimator's top-K over [2 .. max_period].
    if (cfg->period_present) {
        scratch.periods[0] = cfg->period;
        scratch.n_periods = 1;
        if (cfg->verbose) printf("\nbifid: period pinned to %d\n", cfg->period);
    } else {
        int max_p = (cfg->max_period > 0) ? cfg->max_period : 20;
        if (max_p > cipher_len / 2) max_p = cipher_len / 2;
        int n_want = cfg->n_periods;
        scratch.n_periods = bifid_estimate_periods(cipher_indices, cipher_len, 2, max_p,
            n_want, scratch.periods, cfg->verbose);
        if (scratch.n_periods < 1) { scratch.periods[0] = 2; scratch.n_periods = 1; }
    }

    if (cfg->verbose)
        printf("\nbifid: %d positions, %d-letter alphabet %s, %d candidate period(s)\n",
            cipher_len, g_alpha, g_idx_to_char_arr, scratch.n_periods);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // bifid_report fills it (may be NULL for CLI use)

    run_solver(&BIFID_MODEL, &ctx);
}
