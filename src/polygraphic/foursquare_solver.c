#include "foursquare_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Four-Square solver (TYPE foursquare)
// =====================================================================
//
// Four-Square is a digraphic substitution over four 5x5 squares (J->I, so the binary forces
// g_alpha == 25): the upper-left and lower-right are the FIXED standard square, the upper-
// right (UR) and lower-left (LL) are the keyed unknowns. The search state is therefore the
// pair (UR, LL) -- two independent permutations of 0..24 packed back-to-back in st->key
// (UR = key[0..24], LL = key[25..49]) -- hill-climbed / annealed with n-gram scoring, like
// the single-square Playfair break but over twice the cells. Each move perturbs ONE of the
// two squares (chosen uniformly) with the classic Playfair move set: a single cell swap
// (dominant) plus row/column swaps and whole-square reflections. No anti-collapse penalty is
// needed (every square is a bijection, so the whole digraph map is), so the model leaves
// score_adjust at 0 and rides the generic state_score (n-gram + crib). Like Playfair it
// effectively needs -logprob; the 50-cell state needs more text than single-square Playfair.

typedef struct {
    int grid_size;     // cells per square (== g_alpha, SQUARE_GRID)
    int side;          // square side (== SQUARE_SIDE)
} FourSquareScratch;

// One config: both keyed squares are climbed at once. period carries the TOTAL state length
// (2 * grid_size), which seed/copy_state use; the per-square geometry is in the scratch.
static int foursquare_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const FourSquareScratch *p = (const FourSquareScratch *) ctx->model_scratch;
    if (cap < 1) return 0;
    out[0].period = 2 * p->grid_size;
    out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Seed: two independent uniformly-random keyed squares (a Fisher-Yates shuffle per block).
static void foursquare_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const FourSquareScratch *p = (const FourSquareScratch *) ctx->model_scratch;
    (void) cc;
    int g = p->grid_size;
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
// 1% flip rows; 1% flip columns. (Cyclic row/column ROTATIONS are excluded -- they leave
// the cipher unchanged.)
static void foursquare_perturb_block(int *blk, int side, int n) {
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

// Neighbour move: perturb ONE of the two keyed squares, chosen uniformly.
static void foursquare_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                               SolverState *st, bool *force_primary) {
    const FourSquareScratch *p = (const FourSquareScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int g = p->grid_size;
    int s = rand_int(0, 2);
    foursquare_perturb_block(st->key + s * g, p->side, g);
}

static void foursquare_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}

static void foursquare_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                    SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    const FourSquareScratch *p = (const FourSquareScratch *) ctx->model_scratch;
    foursquare_decrypt(ctx->cipher, ctx->cipher_len, st->key, st->key + p->grid_size,
                       p->side, out);
    *score_adjust = 0.0;
}

// Render one side x side square as an indented box of letters.
static void foursquare_print_grid(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void foursquare_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc;
    const FourSquareScratch *p = (const FourSquareScratch *) ctx->model_scratch;
    printf("\n  upper-right square:\n");
    foursquare_print_grid(st->key, p->side);
    printf("  lower-left square:\n");
    foursquare_print_grid(st->key + p->grid_size, p->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "foursquare");
}

static void foursquare_report(const SolverCtx *ctx, const SolverConfig *cc,
                              const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const FourSquareScratch *p = (const FourSquareScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, side = p->side, n = p->grid_size;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The two recovered keyed squares (upper-right, lower-left), read row-major. The
    // plaintext squares (upper-left, lower-right) are the fixed standard alphabet.
    char ur[SQUARE_GRID + 1], ll[SQUARE_GRID + 1];
    for (int i = 0; i < n; i++) { ur[i] = index_to_char(st->key[i]); ll[i] = index_to_char(st->key[n + i]); }
    ur[n] = '\0'; ll[n] = '\0';

    printf("\nResult Score: %.2f | Words: %d | UR=%s | LL=%s\n",
        score, n_words_found, ur, ll);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered upper-right square (row major):\n");
    foursquare_print_grid(st->key, side);
    printf("recovered lower-left square (row major):\n");
    foursquare_print_grid(st->key + n, side);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, UR=..., LL=..., file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, UR=%s, LL=%s, ", score, n_words_found, cfg->cipher_type, ur, ll);
    else
        printf(">>> %.2f, %d, UR=%s, LL=%s, ", score, cfg->cipher_type, ur, ll);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel FOURSQUARE_MODEL = {
    .name = "foursquare", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = foursquare_enumerate, .key_len = NULL,
    .seed = foursquare_seed, .perturb = foursquare_perturb, .copy_state = foursquare_copy,
    .decrypt = foursquare_decrypt_hook, .report = foursquare_report,
    .report_verbose = foursquare_report_verbose,
};

void solve_foursquare(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Four-Square needs a 25-letter (5x5) alphabet; the binary's main forces this by
    // excluding a letter, but guard in case solve_foursquare is driven directly.
    int side = (int) (sqrt((double) g_alpha) + 0.5);
    if (side != SQUARE_SIDE || side * side != g_alpha) {
        printf("\n\nERROR: Four-Square needs a 25-letter alphabet (got %d). "
               "Exclude one letter so it is 25 (e.g. -excludeletter J).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Four-Square solve.\n\n");
        return;
    }
    if (cipher_len % 2 != 0)
        printf("\nWARNING: odd ciphertext length (%d); a Four-Square ciphertext is always "
               "even. The trailing letter is left undecrypted.\n", cipher_len);
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Four-Square ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    FourSquareScratch scratch;
    scratch.grid_size = g_alpha;
    scratch.side = side;

    if (cfg->verbose)
        printf("\nfoursquare: %d positions (%d digraphs), %d-letter alphabet %s\n",
            cipher_len, cipher_len / 2, g_alpha, g_idx_to_char_arr);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // foursquare_report fills it (may be NULL for CLI use)

    run_solver(&FOURSQUARE_MODEL, &ctx);
}
