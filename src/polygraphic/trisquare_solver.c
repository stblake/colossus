#include "trisquare_solver.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Tri-Square solver (TYPE trisquare)
// =====================================================================
//
// Tri-Square is a digraphic substitution over THREE independent keyed 5x5 squares (J->I, so
// the binary forces g_alpha == 25). A plaintext digraph (p1 in sq1, p2 in sq2) becomes a
// ciphertext TRIGRAPH (see trisquare.c) -- a 3:2 length expansion, so an N-letter plaintext
// yields 3N/2 ciphertext letters. The search state is the triple (sq1, sq2, sq3), three
// independent permutations of 0..24 packed back-to-back in st->key (sq1 = key[0..24],
// sq2 = key[25..49], sq3 = key[50..74]) -- the largest square state of the family. It is
// hill-climbed / annealed with n-gram scoring exactly like Four-Square / CM-Bifid: each move
// perturbs ONE of the three squares (chosen uniformly) with the classic Playfair move set (a
// single cell swap, dominant, plus row/column swaps and whole-square reflections). There is
// no square-independent decoupling reward -- every square is a bijection, so any triple
// decrypts to a permuted-but-valid stream and only JOINT correctness yields English n-grams
// -- so the model leaves score_adjust at 0 and rides the generic state_score. Like the other
// square types it effectively needs -logprob.
//
// The 3:2 length change is handled like ADFGVX: the raw ciphertext (a stream of 3M trigraph
// letters) lives in the scratch, the solver passes the PLAINTEXT/scoring length n = 2M to
// make_solver_ctx, and the decrypt hook reads the trigraphs from scratch and emits n plaintext
// symbols that the engine n-gram-scores. Cribs are not used (positions are over the trigraph
// stream and the encode is polyphonic).

typedef struct {
    int  grid_size;    // cells per square (== g_alpha, SQUARE_GRID)
    int  side;         // square side (== SQUARE_SIDE)
    int  len3;         // ciphertext length (== 3M trigraph letters)
    int  m;            // trigraph / digraph count (== len3 / 3)
    int  n;            // plaintext / scoring length (== 2M)
    const int *cipher; // the ciphertext trigraph stream (length len3)
} TriSquareScratch;

// One config: all three keyed squares are climbed at once. period carries the TOTAL state
// length (3 * grid_size); the per-square geometry and the ciphertext are in the scratch.
static int trisquare_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const TriSquareScratch *p = (const TriSquareScratch *) ctx->model_scratch;
    if (cap < 1) return 0;
    out[0].period = 3 * p->grid_size;
    out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Seed: three independent uniformly-random keyed squares (a Fisher-Yates shuffle per block).
static void trisquare_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const TriSquareScratch *p = (const TriSquareScratch *) ctx->model_scratch;
    (void) cc;
    int g = p->grid_size;
    for (int s = 0; s < 3; s++) {
        int *blk = st->key + s * g;
        for (int i = 0; i < g; i++) blk[i] = i;
        for (int i = g - 1; i > 0; i--) {
            int j = rand_int(0, i + 1);
            int t = blk[i]; blk[i] = blk[j]; blk[j] = t;
        }
    }
    st->key_len = 3 * g;
}

// Apply one Playfair-style move to a single side x side square `blk` (n = side*side cells):
// 80% swap two cells; 8% swap two rows; 8% swap two columns; 2% reverse (rotate 180);
// 1% flip rows; 1% flip columns. (Cyclic row/column ROTATIONS are excluded -- they leave
// the cipher unchanged.) Identical to the Four-Square / Bifid / Playfair move set.
static void trisquare_perturb_block(int *blk, int side, int n) {
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

// Neighbour move: perturb ONE of the three keyed squares, chosen uniformly.
static void trisquare_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                              SolverState *st, bool *force_primary) {
    const TriSquareScratch *p = (const TriSquareScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int g = p->grid_size;
    int s = rand_int(0, 3);
    trisquare_perturb_block(st->key + s * g, p->side, g);
}

static void trisquare_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}

static void trisquare_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                   SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    const TriSquareScratch *p = (const TriSquareScratch *) ctx->model_scratch;
    int g = p->grid_size;
    trisquare_decrypt(p->cipher, p->len3, st->key, st->key + g, st->key + 2 * g, p->side, out);
    *score_adjust = 0.0;
}

// Render one side x side square as an indented box of letters.
static void trisquare_print_grid(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void trisquare_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) decrypted;
    const TriSquareScratch *p = (const TriSquareScratch *) ctx->model_scratch;
    int g = p->grid_size;
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  score=%.4f  [%.1fs, %d restarts]\n", score, elapsed, stats->n_restarts);
    printf("  square 1:\n"); trisquare_print_grid(st->key, p->side);
    printf("  square 2:\n"); trisquare_print_grid(st->key + g, p->side);
    printf("  square 3:\n"); trisquare_print_grid(st->key + 2 * g, p->side);
    fflush(stdout);
}

static void trisquare_report(const SolverCtx *ctx, const SolverConfig *cc,
                             const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const TriSquareScratch *p = (const TriSquareScratch *) ctx->model_scratch;
    int n = p->n, side = p->side, g = p->grid_size;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The three recovered keyed squares, read row-major. They are unique only up to the
    // structural symmetries of the cipher (sq1's columns / sq2's rows are unconstrained by
    // the c0/c2 redundancy), but the recovered PLAINTEXT is unique.
    char s1[SQUARE_MAX_GRID + 1], s2[SQUARE_MAX_GRID + 1], s3[SQUARE_MAX_GRID + 1];
    for (int i = 0; i < g; i++) {
        s1[i] = index_to_char(st->key[i]);
        s2[i] = index_to_char(st->key[g + i]);
        s3[i] = index_to_char(st->key[2 * g + i]);
    }
    s1[g] = s2[g] = s3[g] = '\0';

    printf("\nResult Score: %.2f | Words: %d | sq1=%s | sq2=%s | sq3=%s\n",
        score, n_words_found, s1, s2, s3);

    print_cipher((int *) p->cipher, p->len3, NULL);
    printf("\n");
    print_text(decrypted, n);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, n);

    printf("\nrecovered square 1 (row major):\n");
    trisquare_print_grid(st->key, side);
    printf("recovered square 2 (row major):\n");
    trisquare_print_grid(st->key + g, side);
    printf("recovered square 3 (row major):\n");
    trisquare_print_grid(st->key + 2 * g, side);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(decrypted, ctx->result->decrypted, n);
        ctx->result->decrypted_len = n;
    }

    // One-liner summary: >>> score, [words,] type, sq1=..., sq2=..., sq3=..., file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, sq1=%s, sq2=%s, sq3=%s, ",
            score, n_words_found, cfg->cipher_type, s1, s2, s3);
    else
        printf(">>> %.2f, %d, sq1=%s, sq2=%s, sq3=%s, ", score, cfg->cipher_type, s1, s2, s3);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher((int *) p->cipher, p->len3, NULL);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

static const CipherModel TRISQUARE_MODEL = {
    .name = "trisquare", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = trisquare_enumerate, .key_len = NULL,
    .seed = trisquare_seed, .perturb = trisquare_perturb, .copy_state = trisquare_copy,
    .decrypt = trisquare_decrypt_hook, .report = trisquare_report,
    .report_verbose = trisquare_report_verbose,
};

void solve_trisquare(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Tri-Square needs a 25-letter (5x5) alphabet; the binary's main forces this by excluding
    // a letter, but guard in case solve_trisquare is driven directly.
    int side = (int) (sqrt((double) g_alpha) + 0.5);
    if (side != SQUARE_SIDE || side * side != g_alpha) {
        printf("\n\nERROR: Tri-Square needs a 25-letter alphabet (got %d). "
               "Exclude one letter so it is 25 (e.g. -excludeletter J).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 6) {
        printf("\n\nERROR: ciphertext too short for a Tri-Square solve.\n\n");
        return;
    }
    if (cipher_len % 3 != 0) {
        printf("\n\nERROR: Tri-Square ciphertext length must be a multiple of 3 (it is a "
               "stream of trigraphs, 3 cipher letters per plaintext digraph); got %d.\n\n",
               cipher_len);
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Tri-Square ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    int m = cipher_len / 3;    // trigraphs
    int n = 2 * m;             // plaintext / scoring length

    TriSquareScratch scratch;
    scratch.grid_size = g_alpha;
    scratch.side = side;
    scratch.len3 = 3 * m;
    scratch.m = m;
    scratch.n = n;
    scratch.cipher = cipher_indices;

    if (cfg->verbose)
        printf("\ntrisquare: %d ciphertext letters (%d trigraphs, %d plaintext), "
               "%d-letter alphabet %s\n", cipher_len, m, n, g_alpha, g_idx_to_char_arr);

    // Cribs are over the trigraph stream (and the encode is polyphonic), so they do not map to
    // plaintext positions; Tri-Square ignores cribs (like ADFGVX / Bifid).
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    // Pass the PLAINTEXT length n as the engine's scoring length (the decrypt hook reads the
    // trigraphs from the scratch and emits n symbols), mirroring ADFGVX's 2:1 handling.
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, n, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // trisquare_report fills it (may be NULL for CLI use)

    run_solver(&TRISQUARE_MODEL, &ctx);
}
