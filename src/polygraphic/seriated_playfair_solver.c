#include "seriated_playfair_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Seriated Playfair solver (TYPE seriated-playfair)
// =====================================================================
//
// Seriated Playfair (ACA) is plain Playfair over a single 5x5 keyed square, but the
// digraphs are the VERTICAL PAIRS of a two-row seriated layout of period P: in each 2P
// block, pair j couples block-letter j (top) with j+P (bottom) (seriated_playfair.c).
//
// Unlike Portax / Slidefair there is NO per-column independence -- one square enciphers
// every pair -- so this is NOT a per-column key search but a single-grid attack, exactly
// Playfair's: the state is the grid (a permutation of 0..24 in st->key), seeded random
// and annealed with the Playfair move set (a cell swap dominant, plus row/column swaps
// and whole-grid reflections; no anti-collapse penalty -- a grid is a bijection). The ONE
// addition over Playfair is the seriation period P: IoC estimation is useless through the
// vertical pairing, so P is SWEPT (one engine config per period, P in st->cc->period) and
// the n-gram score picks the true one -- a wrong P pairs the wrong letters and decrypts to
// gibberish (only the exact P pairs correctly; a multiple does not). Each config is a full
// grid anneal, so the sweep multiplies the Playfair cost; like every square type it
// effectively needs -logprob. Cribs are not used (the prepare null-insertion shifts
// positions, as in Playfair). The grid is unique only up to a cyclic row/column rotation
// (all re-decrypt identically); the recovered plaintext is unique.

#define SERPF_DEFAULT_MAXP 15   // default top of the period sweep when -maxcols is at its default

typedef struct {
    int grid_size;     // == g_alpha (PLAYFAIR_GRID), the cells climbed
    int side;          // == PLAYFAIR_SIDE
    int minP, maxP;    // seriation-period sweep bounds
} SeriatedPlayfairScratch;

// One engine config per swept seriation period P (period carries P; the grid is always
// grid_size cells). j/k/aux unused.
static int serpf_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const SeriatedPlayfairScratch *a = (const SeriatedPlayfairScratch *) ctx->model_scratch;
    int c = 0;
    for (int P = a->minP; P <= a->maxP && c < cap; P++) {
        out[c].period = P; out[c].j = 0; out[c].k = 0;
        out[c].aux[0] = 0; out[c].aux[1] = 0;
        c++;
    }
    return c;
}

// Seed: a uniformly random grid (Fisher-Yates of 0..grid_size-1). aux[0] caches the grid
// size for copy_state (which gets only the SolverConfig, whose period is the seriation P).
static void serpf_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    const SeriatedPlayfairScratch *a = (const SeriatedPlayfairScratch *) ctx->model_scratch;
    int n = a->grid_size;
    for (int i = 0; i < n; i++) st->key[i] = i;
    for (int i = n - 1; i > 0; i--) {
        int j = rand_int(0, i + 1);
        int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
    }
    st->key_len = n;
    st->aux[0] = n;
}

// Neighbour move (Playfair's): 80% swap two cells; 8% swap two rows; 8% swap two columns;
// 2% reverse the whole grid; 1% flip rows; 1% flip columns. The row/column swaps and
// reflections are genuine jumps (cyclic ROTATIONS are excluded -- they re-decrypt
// identically). Operates on the grid (grid_size cells), independent of the seriation P.
static void serpf_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                          SolverState *st, bool *force_primary) {
    const SeriatedPlayfairScratch *a = (const SeriatedPlayfairScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int s = a->side, n = a->grid_size;
    double r = frand();
    if (r < 0.80) {                              // swap two cells
        int x = rand_int(0, n), y = rand_int(0, n);
        int t = st->key[x]; st->key[x] = st->key[y]; st->key[y] = t;
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

static void serpf_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int n = src->aux[0];
    for (int i = 0; i < n; i++) dst->key[i] = src->key[i];
    dst->key_len = n;
    dst->aux[0] = n;
}

static void serpf_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                               SolverState *st, int *out, double *score_adjust) {
    seriated_playfair_decrypt(ctx->cipher, ctx->cipher_len, st->key, cc->period, out);
    *score_adjust = 0.0;
}

// --- reporting --------------------------------------------------------------------

static void serpf_print_grid(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void serpf_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const SeriatedPlayfairScratch *a = (const SeriatedPlayfairScratch *) ctx->model_scratch;
    printf("\n  P=%d grid:\n", cc->period);
    serpf_print_grid(st->key, a->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "seriated-playfair");
}

static void serpf_report(const SolverCtx *ctx, const SolverConfig *cc,
                         const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const SeriatedPlayfairScratch *a = (const SeriatedPlayfairScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, side = a->side, n = a->grid_size, P = cc->period;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The recovered key square, read row-major (unique only up to cyclic row/column
    // rotation, which all decrypt identically -- this is one representative).
    char gridstr[PLAYFAIR_GRID + 1];
    for (int i = 0; i < n; i++) gridstr[i] = index_to_char(st->key[i]);
    gridstr[n] = '\0';

    printf("\nResult Score: %.2f | Words: %d | P=%d | grid=%s\n", score, n_words_found, P, gridstr);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered 5x5 grid (row major), seriation period %d:\n", P);
    serpf_print_grid(st->key, side);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = P;          // the seriation period (tests read this)
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, P=..., grid=..., file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, P=%d, grid=%s, ", score, n_words_found, cfg->cipher_type, P, gridstr);
    else
        printf(">>> %.2f, %d, P=%d, grid=%s, ", score, cfg->cipher_type, P, gridstr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel SERIATED_PLAYFAIR_MODEL = {
    .name = "seriated-playfair", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = serpf_enumerate, .key_len = NULL,
    .seed = serpf_seed, .perturb = serpf_perturb, .copy_state = serpf_copy,
    .decrypt = serpf_decrypt_hook, .report = serpf_report,
    .report_verbose = serpf_report_verbose,
};

// =====================================================================
//  Entry point
// =====================================================================

void solve_seriated_playfair(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Needs the 25-letter (5x5) grid; main() forces it (J->I) for this type, but guard
    // in case solve_seriated_playfair is driven directly.
    int side = (int) (sqrt((double) g_alpha) + 0.5);
    if (side != PLAYFAIR_SIDE || side * side != g_alpha) {
        printf("\n\nERROR: Seriated Playfair needs a 25-letter alphabet (got %d). "
               "Exclude one letter so it is 25 (e.g. -excludeletter J).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Seriated Playfair solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Seriated Playfair ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    // Period sweep: -period pins; else -mincols..-maxcols, defaulting the top to
    // SERPF_DEFAULT_MAXP when -maxcols is left at its global default (periods are small).
    // A period needs at least one full 2P block, so clamp maxP to cipher_len/2.
    SeriatedPlayfairScratch a;
    a.grid_size = g_alpha;
    a.side = side;
    if (cfg->period_present) { a.minP = a.maxP = cfg->period; }
    else {
        a.minP = cfg->min_cols < 1 ? 1 : cfg->min_cols;
        if (a.minP < 1) a.minP = 1;
        a.maxP = (cfg->max_cols == 30) ? SERPF_DEFAULT_MAXP : cfg->max_cols;
    }
    if (a.maxP > MAX_COLS) a.maxP = MAX_COLS;
    if (a.maxP > cipher_len / 2) a.maxP = cipher_len / 2;
    if (a.maxP < a.minP) a.maxP = a.minP;

    if (cfg->verbose)
        printf("\nseriated-playfair: %d positions, single-grid anneal, seriation P %d..%d\n",
            cipher_len, a.minP, a.maxP);

    // Cribs are not used (the prepare null-insertion shifts plaintext positions); the crib
    // args are still threaded through make_solver_ctx for a uniform interface.
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &a;
    ctx.result = result;

    run_solver(&SERIATED_PLAYFAIR_MODEL, &ctx);
}
