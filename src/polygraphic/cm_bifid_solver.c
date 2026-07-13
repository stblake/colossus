#include "cm_bifid_solver.h"
#include "bifid_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  CM Bifid solver (TYPE cm-bifid)
// =====================================================================
//
// CM Bifid (Conjugated Matrix Bifid) is plain Bifid fractionation with TWO keyed Polybius
// squares: square 1 fractionates the plaintext into coords, and after the rows-then-cols
// reshape and re-pairing the coordinate pairs are mapped to ciphertext letters through a
// DIFFERENT square 2 (see cm_bifid.c). Breaking it couples three problems: the period and
// the TWO squares.
//
// The square search is the proven joint two-square anneal (cf. Two/Four-Square): the state
// is the pair of squares -- two independent permutations of 0..n-1 packed back-to-back in
// st->key (sq1 = key[0..g-1], sq2 = key[g..2g-1]) -- annealed with n-gram scoring, each
// move perturbing ONE square (chosen uniformly) with the Playfair move set (cell-swap-
// dominant + row/column swaps + reflections). There is NO square-independent decoupling
// reward (both squares are entangled in the n-gram fitness -- unlike ADFGVX's transposition-
// only IoC or Nihilist-sub's additive-only validity), so this is a genuine joint search and,
// like every square type, effectively needs -logprob. No anti-collapse penalty is needed
// (each square is a bijection, so the whole digraph map is), so score_adjust stays 0.
//
// The period is recovered exactly as Bifid's -- the columnar Index of Coincidence is square-
// AGNOSTIC (square 2 only relabels the coordinate pairs, and IoC is relabel-invariant), so
// bifid_estimate_periods is reused unchanged: we anneal its top-K periods (one engine config
// each) and the n-gram score picks the winner. -period pins a single period.

#define CM_BIFID_MAX_PERIODS 64

typedef struct {
    int side;                          // grid side (5 for the 25-letter default, up to 6)
    int grid_size;                     // side*side == g_alpha == cells per square
    int n_periods;                     // number of candidate periods
    int periods[CM_BIFID_MAX_PERIODS]; // the candidate periods (config order)
} CMBifidScratch;

// One config per candidate period; period carries the block size, key the two squares
// (always 2*grid_size cells -- the per-square geometry is in the scratch).
static int cm_bifid_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const CMBifidScratch *b = (const CMBifidScratch *) ctx->model_scratch;
    int n = b->n_periods;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = b->periods[i];
        out[i].j = 0; out[i].k = 0; out[i].aux[0] = 0; out[i].aux[1] = 0;
    }
    return n;
}

// Seed: two independent uniformly-random squares (a Fisher-Yates shuffle per block).
static void cm_bifid_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    const CMBifidScratch *b = (const CMBifidScratch *) ctx->model_scratch;
    int g = b->grid_size;
    for (int s = 0; s < 2; s++) {
        int *blk = st->key + s * g;
        for (int i = 0; i < g; i++) blk[i] = i;
        for (int i = g - 1; i > 0; i--) {
            int j = rand_int(0, i + 1);
            int t = blk[i]; blk[i] = blk[j]; blk[j] = t;
        }
    }
    st->key_len = 2 * g;
}

// Apply one Playfair-style move to a single side x side square `blk` (n = side*side cells):
// 80% swap two cells; 8% swap two rows; 8% swap two columns; 2% reverse (rotate 180);
// 1% flip rows; 1% flip columns. The larger moves jump the basins a cell swap cannot escape.
static void cm_bifid_perturb_block(int *blk, int side, int n) {
    double r = frand();
    if (r < 0.80) {                              // swap two cells
        int a = rand_int(0, n), b = rand_int(0, n);
        int t = blk[a]; blk[a] = blk[b]; blk[b] = t;
    } else if (r < 0.88) {                       // swap two rows
        int r1 = rand_int(0, side), r2 = rand_int(0, side);
        for (int c = 0; c < side; c++) {
            int t = blk[r1 * side + c]; blk[r1 * side + c] = blk[r2 * side + c]; blk[r2 * side + c] = t;
        }
    } else if (r < 0.96) {                       // swap two columns
        int c1 = rand_int(0, side), c2 = rand_int(0, side);
        for (int rr = 0; rr < side; rr++) {
            int t = blk[rr * side + c1]; blk[rr * side + c1] = blk[rr * side + c2]; blk[rr * side + c2] = t;
        }
    } else if (r < 0.98) {                       // reverse the whole square (rotate 180)
        for (int i = 0, j = n - 1; i < j; i++, j--) {
            int t = blk[i]; blk[i] = blk[j]; blk[j] = t;
        }
    } else if (r < 0.99) {                       // flip rows top<->bottom
        for (int r1 = 0, r2 = side - 1; r1 < r2; r1++, r2--)
            for (int c = 0; c < side; c++) {
                int t = blk[r1 * side + c]; blk[r1 * side + c] = blk[r2 * side + c]; blk[r2 * side + c] = t;
            }
    } else {                                     // flip columns left<->right
        for (int c1 = 0, c2 = side - 1; c1 < c2; c1++, c2--)
            for (int rr = 0; rr < side; rr++) {
                int t = blk[rr * side + c1]; blk[rr * side + c1] = blk[rr * side + c2]; blk[rr * side + c2] = t;
            }
    }
}

// Neighbour move: perturb ONE of the two squares, chosen uniformly.
static void cm_bifid_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                             SolverState *st, bool *force_primary) {
    const CMBifidScratch *b = (const CMBifidScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int g = b->grid_size;
    int s = rand_int(0, 2);
    cm_bifid_perturb_block(st->key + s * g, b->side, g);
}

static void cm_bifid_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    // Two squares; copy a fixed upper bound (cc->period is the BLOCK size, not the state
    // length) so the engine need not know the side. 2*BIFID_MAX_GRID = 72 cells.
    for (int i = 0; i < 2 * BIFID_MAX_GRID; i++) dst->key[i] = src->key[i];
}

static void cm_bifid_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                  SolverState *st, int *out, double *score_adjust) {
    const CMBifidScratch *b = (const CMBifidScratch *) ctx->model_scratch;
    cm_bifid_decrypt(ctx->cipher, ctx->cipher_len, st->key, st->key + b->grid_size,
                     b->side, cc->period, out);
    *score_adjust = 0.0;
}

// Render one side x side square as an indented box of letters.
static void cm_bifid_print_grid(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void cm_bifid_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const CMBifidScratch *b = (const CMBifidScratch *) ctx->model_scratch;
    printf("\n  period %d, square 1 (fractionation):\n", cc->period);
    cm_bifid_print_grid(st->key, b->side);
    printf("  square 2 (recombination):\n");
    cm_bifid_print_grid(st->key + b->grid_size, b->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "cm-bifid");
}

static void cm_bifid_report(const SolverCtx *ctx, const SolverConfig *cc,
                            const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const CMBifidScratch *b = (const CMBifidScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, side = b->side, n = b->grid_size;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The two recovered squares, read row-major. Each square is unique only up to a cyclic
    // row/column rotation (which all re-decrypt identically); the recovered PLAINTEXT is unique.
    char sq1[BIFID_MAX_GRID + 1], sq2[BIFID_MAX_GRID + 1];
    for (int i = 0; i < n; i++) { sq1[i] = index_to_char(st->key[i]); sq2[i] = index_to_char(st->key[n + i]); }
    sq1[n] = '\0'; sq2[n] = '\0';

    printf("\nResult Score: %.2f | Words: %d | period=%d | sq1=%s | sq2=%s\n",
        score, n_words_found, cc->period, sq1, sq2);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered %dx%d square 1 (fractionation, row major):\n", side, side);
    cm_bifid_print_grid(st->key, side);
    printf("recovered %dx%d square 2 (recombination, row major):\n", side, side);
    cm_bifid_print_grid(st->key + n, side);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = cc->period;   // report the recovered period here
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, period=, sq1=, sq2=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, period=%d, sq1=%s, sq2=%s, ",
            score, n_words_found, cfg->cipher_type, cc->period, sq1, sq2);
    else
        printf(">>> %.2f, %d, period=%d, sq1=%s, sq2=%s, ",
            score, cfg->cipher_type, cc->period, sq1, sq2);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel CM_BIFID_MODEL = {
    .name = "cm-bifid", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = cm_bifid_enumerate, .key_len = NULL,
    .seed = cm_bifid_seed, .perturb = cm_bifid_perturb, .copy_state = cm_bifid_copy,
    .decrypt = cm_bifid_decrypt_hook, .report = cm_bifid_report,
    .report_verbose = cm_bifid_report_verbose,
};

void solve_cm_bifid(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // CM Bifid needs a perfect-square (e.g. 5x5 -> 25, 6x6 -> 36) alphabet; the binary
    // forces 25 for the default, but guard in case solve_cm_bifid is driven directly.
    int side = exact_isqrt(g_alpha);
    if (side < 2 || side > BIFID_MAX_SIDE || side * side != g_alpha) {
        printf("\n\nERROR: CM Bifid needs a perfect-square alphabet (e.g. 25 or 36; got %d). "
               "Exclude a letter so it is 25 (e.g. -excludeletter J).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a CM Bifid solve.\n\n");
        return;
    }
    // CM Bifid letters must be solid: reject any sentinel (space/punctuation).
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "CM Bifid ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    CMBifidScratch scratch;
    scratch.side = side;
    scratch.grid_size = g_alpha;

    // Candidate periods: pinned, or the (Bifid) estimator's top-K over [2 .. max_period].
    if (cfg->period_present) {
        scratch.periods[0] = cfg->period;
        scratch.n_periods = 1;
        if (cfg->verbose) printf("\ncm-bifid: period pinned to %d\n", cfg->period);
    } else {
        int max_p = (cfg->max_period > 0) ? cfg->max_period : 20;
        if (max_p > cipher_len / 2) max_p = cipher_len / 2;
        int n_want = cfg->n_periods;
        scratch.n_periods = bifid_estimate_periods(cipher_indices, cipher_len, 2, max_p,
            n_want, scratch.periods, cfg->verbose);
        if (scratch.n_periods < 1) { scratch.periods[0] = 2; scratch.n_periods = 1; }
    }

    if (cfg->verbose)
        printf("\ncm-bifid: %d positions, %d-letter alphabet %s, %d candidate period(s)\n",
            cipher_len, g_alpha, g_idx_to_char_arr, scratch.n_periods);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // cm_bifid_report fills it (may be NULL for CLI use)

    run_solver(&CM_BIFID_MODEL, &ctx);
}
