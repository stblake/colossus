#include "phillips_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Phillips solver (TYPE phillips / phillips-c / phillips-rc)
// =====================================================================
//
// Phillips is a periodic monographic substitution over a 5x5 keyed grid (J->I, so the
// binary forces g_alpha == 25). The plaintext is split into blocks of 5 letters; block b
// is enciphered with one of nsq = 8 squares derived from the base grid (square b mod 8),
// each letter going to its down-right neighbour (overall period 40). See phillips.c for
// the square-generation variants (PHILLIPS_ROW / _COL / _ROWCOL, chosen by cipher_type).
//
// The only unknown is the base grid -- a permutation of 0..24 carried in st->key -- so the
// search is identical to Playfair's (the classic SA square break): the move set is a single
// cell swap (dominant) plus row/column swaps and whole-grid reflections, which jump the
// local optima a cell swap gets stuck in. There is no period to estimate (block size and the
// 8-square cycle are fixed -> one config) and no anti-collapse penalty is needed (every
// derived square, and hence the whole map, is a bijection), so the model leaves score_adjust
// at 0 and rides the generic state_score (n-gram + crib). Like Playfair it effectively needs
// -logprob with higher-order n-grams.

typedef struct {
    int grid_size;     // == g_alpha (PHILLIPS_GRID)
    int side;          // == PHILLIPS_SIDE
    int variant;       // PHILLIPS_ROW / PHILLIPS_COL / PHILLIPS_ROWCOL
} PhillipsScratch;

// Map a Phillips cipher type to its square-generation variant.
static int phillips_variant_for_type(int cipher_type) {
    switch (cipher_type) {
        case PHILLIPS_C:  return PHILLIPS_COL;
        case PHILLIPS_RC: return PHILLIPS_ROWCOL;
        default:          return PHILLIPS_ROW;
    }
}

static const char *phillips_variant_name(int variant) {
    switch (variant) {
        case PHILLIPS_COL:    return "column";
        case PHILLIPS_ROWCOL: return "row-column";
        default:              return "row";
    }
}

// One config: the whole grid is climbed at once. period carries the key length.
static int phillips_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const PhillipsScratch *p = (const PhillipsScratch *) ctx->model_scratch;
    if (cap < 1) return 0;
    out[0].period = p->grid_size;
    out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Seed: a uniformly random grid (Fisher-Yates shuffle of 0..n-1 via the engine RNG).
static void phillips_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    int n = cc->period;
    for (int i = 0; i < n; i++) st->key[i] = i;
    for (int i = n - 1; i > 0; i--) {
        int j = rand_int(0, i + 1);
        int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
    }
    st->key_len = n;
}

// Neighbour move (identical to Playfair's): 80% swap two cells; 8% swap two rows; 8% swap
// two columns; 2% reverse the whole grid; 1% flip rows; 1% flip columns. The row/column
// swaps and reflections are genuine jumps (they change the derived squares, so they decrypt
// differently) and let annealing escape the basins a single cell swap cannot leave.
static void phillips_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                             SolverState *st, bool *force_primary) {
    const PhillipsScratch *p = (const PhillipsScratch *) ctx->model_scratch;
    (void) force_primary;
    int s = p->side, n = cc->period;
    double r = frand();
    if (r < 0.80) {                              // swap two cells
        int a = rand_int(0, n), b = rand_int(0, n);
        int t = st->key[a]; st->key[a] = st->key[b]; st->key[b] = t;
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

static void phillips_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}

static void phillips_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                  SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    const PhillipsScratch *p = (const PhillipsScratch *) ctx->model_scratch;
    phillips_decrypt(ctx->cipher, ctx->cipher_len, st->key, p->side, p->variant, out);
    *score_adjust = 0.0;
}

// Render the recovered 5x5 grid as an indented box of letters.
static void phillips_print_grid(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void phillips_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc;
    const PhillipsScratch *p = (const PhillipsScratch *) ctx->model_scratch;
    printf("\n  grid:\n");
    phillips_print_grid(st->key, p->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "phillips");
}

static void phillips_report(const SolverCtx *ctx, const SolverConfig *cc,
                            const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const PhillipsScratch *p = (const PhillipsScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, side = p->side, n = p->grid_size;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The recovered base square, read row-major (it is unique only up to a cyclic column
    // rotation, which re-enciphers identically -- this is one representative).
    char gridstr[PHILLIPS_GRID + 1];
    for (int i = 0; i < n; i++) gridstr[i] = index_to_char(st->key[i]);
    gridstr[n] = '\0';

    printf("\nResult Score: %.2f | Words: %d | variant=%s | grid=%s\n",
        score, n_words_found, phillips_variant_name(p->variant), gridstr);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_solution_check(decrypted, len);
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered 5x5 grid (row major):\n");
    phillips_print_grid(st->key, side);

    // Publish the recovered solution for callers that pass a SolveResult (the in-process
    // tests inspect it instead of scraping stdout).
    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, grid=..., file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, grid=%s, ", score, n_words_found, cfg->cipher_type, gridstr);
    else
        printf(">>> %.2f, %d, grid=%s, ", score, cfg->cipher_type, gridstr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel PHILLIPS_MODEL = {
    .name = "phillips", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = phillips_enumerate, .key_len = NULL,
    .seed = phillips_seed, .perturb = phillips_perturb, .copy_state = phillips_copy,
    .decrypt = phillips_decrypt_hook, .report = phillips_report,
    .report_verbose = phillips_report_verbose,
};

void solve_phillips(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Phillips needs a 25-letter (5x5) grid; the binary's main forces this by excluding a
    // letter, but guard in case solve_phillips is driven directly.
    int side = (int) (sqrt((double) g_alpha) + 0.5);
    if (side != PHILLIPS_SIDE || side * side != g_alpha) {
        printf("\n\nERROR: Phillips needs a 25-letter alphabet (got %d). "
               "Exclude one letter so it is 25 (e.g. -excludeletter J).\n\n", g_alpha);
        return;
    }
    if (cipher_len < side) {
        printf("\n\nERROR: ciphertext too short for a Phillips solve.\n\n");
        return;
    }
    // Phillips letters must be solid: reject any sentinel (space/punctuation).
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Phillips ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    PhillipsScratch scratch;
    scratch.grid_size = g_alpha;
    scratch.side = side;
    scratch.variant = phillips_variant_for_type(cfg->cipher_type);

    if (cfg->verbose)
        printf("\nphillips: %d positions, %d-letter alphabet %s, %s variant\n",
            cipher_len, g_alpha, g_idx_to_char_arr, phillips_variant_name(scratch.variant));

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // phillips_report fills it (may be NULL for CLI use)

    run_solver(&PHILLIPS_MODEL, &ctx);
}
