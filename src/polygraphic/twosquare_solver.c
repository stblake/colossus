#include "twosquare_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Two-Square solver (TYPE twosquare / twosquare-v)
// =====================================================================
//
// Two-Square is a digraphic substitution over TWO keyed 5x5 squares (J->I, so the binary
// forces g_alpha == 25). The search state is the pair of squares -- two independent
// permutations of 0..24 packed back-to-back in st->key (sq1 = key[0..24], sq2 = key[25..49])
// -- hill-climbed / annealed with n-gram scoring, exactly like the single-square Playfair
// break but over twice the cells. Each move perturbs ONE of the two squares (chosen
// uniformly) with the classic Playfair move set: a single cell swap (dominant) plus
// row/column swaps and whole-square reflections, which jump the local optima a cell swap
// cannot escape. No anti-collapse penalty is needed (every square is a bijection, so the
// whole digraph map is), so the model leaves score_adjust at 0 and rides the generic
// state_score (n-gram + crib). Like Playfair it effectively needs -logprob.
//
// The `variant` (cipher_type) selects the arrangement -- horizontal (ACA) or vertical
// (Wikipedia, self-inverse) -- which only changes the decrypt rule (see twosquare.c).

typedef struct {
    int grid_size;     // cells per square (== g_alpha, SQUARE_GRID)
    int side;          // square side (== SQUARE_SIDE)
    int variant;       // TWO_SQ_HORIZONTAL / TWO_SQ_VERTICAL
} TwoSquareScratch;

static int twosquare_variant_for_type(int cipher_type) {
    return (cipher_type == TWO_SQUARE_V) ? TWO_SQ_VERTICAL : TWO_SQ_HORIZONTAL;
}

static const char *twosquare_variant_name(int variant) {
    return (variant == TWO_SQ_VERTICAL) ? "vertical" : "horizontal";
}

// One config: both squares are climbed at once. period carries the TOTAL state length
// (2 * grid_size), which seed/copy_state use; the per-square geometry is in the scratch.
static int twosquare_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const TwoSquareScratch *p = (const TwoSquareScratch *) ctx->model_scratch;
    if (cap < 1) return 0;
    out[0].period = 2 * p->grid_size;
    out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Seed: two independent uniformly-random squares (a Fisher-Yates shuffle per block).
static void twosquare_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const TwoSquareScratch *p = (const TwoSquareScratch *) ctx->model_scratch;
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
static void twosquare_perturb_block(int *blk, int side, int n) {
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
static void twosquare_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                              SolverState *st, bool *force_primary) {
    const TwoSquareScratch *p = (const TwoSquareScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int g = p->grid_size;
    int s = rand_int(0, 2);
    twosquare_perturb_block(st->key + s * g, p->side, g);
}

static void twosquare_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}

static void twosquare_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                   SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    const TwoSquareScratch *p = (const TwoSquareScratch *) ctx->model_scratch;
    twosquare_decrypt(ctx->cipher, ctx->cipher_len, st->key, st->key + p->grid_size,
                      p->side, p->variant, out);
    *score_adjust = 0.0;
}

// Render one side x side square as an indented box of letters.
static void twosquare_print_grid(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void twosquare_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc;
    const TwoSquareScratch *p = (const TwoSquareScratch *) ctx->model_scratch;
    printf("\n  square 1:\n");
    twosquare_print_grid(st->key, p->side);
    printf("  square 2:\n");
    twosquare_print_grid(st->key + p->grid_size, p->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "twosquare");
}

static void twosquare_report(const SolverCtx *ctx, const SolverConfig *cc,
                             const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const TwoSquareScratch *p = (const TwoSquareScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, side = p->side, n = p->grid_size;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The two recovered squares, read row-major. The pair is unique only up to the
    // arrangement's symmetries (e.g. for the vertical type a cyclic column rotation of both
    // squares re-enciphers identically); the recovered PLAINTEXT is unique.
    char sq1[SQUARE_GRID + 1], sq2[SQUARE_GRID + 1];
    for (int i = 0; i < n; i++) { sq1[i] = index_to_char(st->key[i]); sq2[i] = index_to_char(st->key[n + i]); }
    sq1[n] = '\0'; sq2[n] = '\0';

    printf("\nResult Score: %.2f | Words: %d | variant=%s | sq1=%s | sq2=%s\n",
        score, n_words_found, twosquare_variant_name(p->variant), sq1, sq2);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_solution_check(decrypted, len);
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered square 1 (row major):\n");
    twosquare_print_grid(st->key, side);
    printf("recovered square 2 (row major):\n");
    twosquare_print_grid(st->key + n, side);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, sq1=..., sq2=..., file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, sq1=%s, sq2=%s, ", score, n_words_found, cfg->cipher_type, sq1, sq2);
    else
        printf(">>> %.2f, %d, sq1=%s, sq2=%s, ", score, cfg->cipher_type, sq1, sq2);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel TWOSQUARE_MODEL = {
    .name = "twosquare", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = twosquare_enumerate, .key_len = NULL,
    .seed = twosquare_seed, .perturb = twosquare_perturb, .copy_state = twosquare_copy,
    .decrypt = twosquare_decrypt_hook, .report = twosquare_report,
    .report_verbose = twosquare_report_verbose,
};

void solve_twosquare(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Two-Square needs a 25-letter (5x5) alphabet; the binary's main forces this by
    // excluding a letter, but guard in case solve_twosquare is driven directly.
    int side = (int) (sqrt((double) g_alpha) + 0.5);
    if (side != SQUARE_SIDE || side * side != g_alpha) {
        printf("\n\nERROR: Two-Square needs a 25-letter alphabet (got %d). "
               "Exclude one letter so it is 25 (e.g. -excludeletter J).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Two-Square solve.\n\n");
        return;
    }
    if (cipher_len % 2 != 0)
        printf("\nWARNING: odd ciphertext length (%d); a Two-Square ciphertext is always "
               "even. The trailing letter is left undecrypted.\n", cipher_len);
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Two-Square ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    TwoSquareScratch scratch;
    scratch.grid_size = g_alpha;
    scratch.side = side;
    scratch.variant = twosquare_variant_for_type(cfg->cipher_type);

    if (cfg->verbose)
        printf("\ntwosquare: %d positions (%d digraphs), %d-letter alphabet %s, %s arrangement\n",
            cipher_len, cipher_len / 2, g_alpha, g_idx_to_char_arr,
            twosquare_variant_name(scratch.variant));

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // twosquare_report fills it (may be NULL for CLI use)

    run_solver(&TWOSQUARE_MODEL, &ctx);
}
