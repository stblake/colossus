#include "playfair_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Playfair solver (TYPE playfair)
// =====================================================================
//
// Playfair is a digraphic substitution over a 5x5 keyed grid of the 25-letter
// alphabet (the binary forces g_alpha == 25 for this type by excluding one letter,
// J->I by default). The search state is the grid itself -- a permutation of 0..24
// carried in st->key -- hill-climbed / annealed with n-gram scoring (the classic
// simulated-annealing attack). The move set is a single cell swap (dominant) plus
// row/column swaps and whole-grid reflections; the larger moves jump across the
// local optima that pure cell swaps get stuck in. No anti-collapse penalty is needed
// (a grid is a bijection -- it cannot pile the alphabet onto a few letters), so the
// model leaves score_adjust at 0 and rides the generic state_score (n-gram + crib).

typedef struct {
    int grid_size;     // == g_alpha (PLAYFAIR_GRID)
    int side;          // == PLAYFAIR_SIDE
} PlayfairScratch;

// One config: the whole grid is climbed at once. period carries the key length.
static int playfair_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const PlayfairScratch *p = (const PlayfairScratch *) ctx->model_scratch;
    if (cap < 1) return 0;
    out[0].period = p->grid_size;
    out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Seed: a uniformly random grid (Fisher-Yates shuffle of 0..n-1 via the engine RNG).
static void playfair_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    int n = cc->period;
    for (int i = 0; i < n; i++) st->key[i] = i;
    for (int i = n - 1; i > 0; i--) {
        int j = rand_int(0, i + 1);
        int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
    }
    st->key_len = n;
}

// Neighbour move. 80% swap two cells (the workhorse); 8% swap two rows; 8% swap two
// columns; 2% reverse the whole grid; 1% flip rows (top<->bottom); 1% flip columns
// (left<->right). The row/column swaps and reflections are genuine jumps (they change
// the down/up-neighbour structure, so they decrypt differently) and let annealing
// escape the basins a single cell swap cannot leave. (Cyclic row/column ROTATIONS are
// deliberately not in the set -- they leave the cipher unchanged, so they would just
// burn iterations.)
static void playfair_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                             SolverState *st, bool *force_primary) {
    const PlayfairScratch *p = (const PlayfairScratch *) ctx->model_scratch;
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

static void playfair_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}

static void playfair_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                  SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    playfair_decrypt(ctx->cipher, ctx->cipher_len, st->key, out);
    *score_adjust = 0.0;
}

// Render the recovered 5x5 grid as an indented box of letters.
static void playfair_print_grid(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void playfair_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc;
    const PlayfairScratch *p = (const PlayfairScratch *) ctx->model_scratch;
    printf("\n  grid:\n");
    playfair_print_grid(st->key, p->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "playfair");
}

static void playfair_report(const SolverCtx *ctx, const SolverConfig *cc,
                            const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const PlayfairScratch *p = (const PlayfairScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, side = p->side, n = p->grid_size;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The recovered key square, read row-major (the grid is unique only up to cyclic
    // row/column rotation, which all decrypt identically -- this is one representative).
    char gridstr[PLAYFAIR_GRID + 1];
    for (int i = 0; i < n; i++) gridstr[i] = index_to_char(st->key[i]);
    gridstr[n] = '\0';

    printf("\nResult Score: %.2f | Words: %d | grid=%s\n", score, n_words_found, gridstr);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered 5x5 grid (row major):\n");
    playfair_print_grid(st->key, side);

    // Publish the recovered solution for callers that pass a SolveResult (the
    // in-process tests inspect it instead of scraping stdout).
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

static const CipherModel PLAYFAIR_MODEL = {
    .name = "playfair", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = playfair_enumerate, .key_len = NULL,
    .seed = playfair_seed, .perturb = playfair_perturb, .copy_state = playfair_copy,
    .decrypt = playfair_decrypt_hook, .report = playfair_report,
    .report_verbose = playfair_report_verbose,
};

void solve_playfair(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Playfair needs a perfect-square (5x5 -> 25) alphabet; the binary's main forces
    // this by excluding a letter, but guard in case solve_playfair is driven directly.
    int side = (int) (sqrt((double) g_alpha) + 0.5);
    if (side != PLAYFAIR_SIDE || side * side != g_alpha) {
        printf("\n\nERROR: Playfair needs a 25-letter alphabet (got %d). "
               "Exclude one letter so it is 25 (e.g. -excludeletter J).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Playfair solve.\n\n");
        return;
    }
    if (cipher_len % 2 != 0)
        printf("\nWARNING: odd ciphertext length (%d); a Playfair ciphertext is always "
               "even. The trailing letter is left undecrypted.\n", cipher_len);
    // Playfair digraphs must be solid letters: reject any sentinel (space/punctuation).
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Playfair ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    if (cfg->verbose)
        printf("\nplayfair: %d positions (%d digraphs), %d-letter alphabet %s\n",
            cipher_len, cipher_len / 2, g_alpha, g_idx_to_char_arr);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);

    PlayfairScratch scratch;
    scratch.grid_size = g_alpha;
    scratch.side = side;
    ctx.model_scratch = &scratch;
    ctx.result = result;          // playfair_report fills it (may be NULL for CLI use)

    run_solver(&PLAYFAIR_MODEL, &ctx);
}


