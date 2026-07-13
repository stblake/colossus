#include "adfgvx_solver.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  ADFGX / ADFGVX solver (TYPE adfgx / adfgvx)
// =====================================================================
//
// ADFGVX is a keyed Polybius-square fractionation (each plaintext symbol -> its two
// cell coordinates) composed with a keyed columnar transposition of the resulting 2N
// coordinate stream (see adfgvx.c). Breaking it is therefore a COUPLED problem: recover
// the square AND the columnar column order at once. The state carries both:
//   key[0 .. grid_size-1]            the keyed square (a permutation of 0..grid_size-1)
//   key[grid_size .. grid_size+K-1]  the columnar read order (a permutation of 0..K-1)
//   aux[0]=K  aux[1]=dir  aux[2]=side  aux[3]=grid_size
// One engine config is enumerated per column count K in [min_cols .. max_cols].
//
// The coupling is the hard part: with a WRONG column order the coordinate stream is
// re-paired across block seams, so the plaintext is gibberish for EVERY square and the
// n-gram score gives the column-order moves no gradient until the square is also right.
// The fix (standard ADFGVX cryptanalysis): after undoing the correct columnar, the
// paired cell ids are a MONOALPHABETIC image of the plaintext, so the decrypted text's
// Index of Coincidence is English (~0.066); a wrong order pairs randomly -> flat
// (~1/grid_size). Crucially the IoC of the decrypt depends ONLY on the column order, not
// on the square (a square just relabels the cells, which IoC is invariant to). Folding a
// structural IoC reward (ADFGVX_IOC_WEIGHT * ioc) into score_adjust thus gives the
// column-order search a gradient independent of the square, decoupling the two halves:
// the climb locks the column order by IoC, then the n-gram score recovers the square.
//
// The square move set is identical to Bifid/Playfair (cell-swap-dominant + row/column
// swaps and reflections); the column-order move set is the columnar solver's (swap
// dominant + short reverse / block-move). No anti-collapse penalty is needed -- the
// square is a bijection and the IoC term is additive reward, not a guard.

// Weight of the structural IoC reward (see above). Chosen so its variation across
// column orders (~0.02-0.04 in IoC) is large enough to drive the column-order search
// over the n-gram noise of a wrong square, while staying a constant offset once the
// column order is fixed (so the n-gram score then decides the square).
#define ADFGVX_IOC_WEIGHT 30.0

typedef struct {
    int  side;        // Polybius square side (5 for ADFGX, 6 for ADFGVX)
    int  grid_size;   // side*side == g_alpha
    int  len2;        // ciphertext coordinate length (== 2N)
    int  n;           // plaintext length (== len2/2), the scoring length
    int  klo, khi;    // column-count (K) sweep range
    int *coords;      // ciphertext as 0..side-1 coordinates (length len2)
} AdfgvxScratch;

// Ciphertext coordinate buffer (single-threaded), filled by solve_adfgvx.
static _Thread_local int g_adfgvx_coords[MAX_CIPHER_LENGTH];

// --- neighbour moves --------------------------------------------------------------

// One move on the square (sq[0..n-1], side x side). Identical to Bifid/Playfair: 80%
// swap two cells; 8% swap two rows; 8% swap two columns; 2% rotate 180; 1% flip rows;
// 1% flip columns.
static void adfgvx_square_move(int sq[], int s, int n) {
    double r = frand();
    if (r < 0.80) {
        int a = rand_int(0, n), c = rand_int(0, n);
        int t = sq[a]; sq[a] = sq[c]; sq[c] = t;
    } else if (r < 0.88) {
        int r1 = rand_int(0, s), r2 = rand_int(0, s);
        for (int c = 0; c < s; c++) {
            int t = sq[r1 * s + c]; sq[r1 * s + c] = sq[r2 * s + c]; sq[r2 * s + c] = t;
        }
    } else if (r < 0.96) {
        int c1 = rand_int(0, s), c2 = rand_int(0, s);
        for (int rr = 0; rr < s; rr++) {
            int t = sq[rr * s + c1]; sq[rr * s + c1] = sq[rr * s + c2]; sq[rr * s + c2] = t;
        }
    } else if (r < 0.98) {
        for (int i = 0, j = n - 1; i < j; i++, j--) { int t = sq[i]; sq[i] = sq[j]; sq[j] = t; }
    } else if (r < 0.99) {
        for (int r1 = 0, r2 = s - 1; r1 < r2; r1++, r2--)
            for (int c = 0; c < s; c++) {
                int t = sq[r1 * s + c]; sq[r1 * s + c] = sq[r2 * s + c]; sq[r2 * s + c] = t;
            }
    } else {
        for (int c1 = 0, c2 = s - 1; c1 < c2; c1++, c2--)
            for (int rr = 0; rr < s; rr++) {
                int t = sq[rr * s + c1]; sq[rr * s + c1] = sq[rr * s + c2]; sq[rr * s + c2] = t;
            }
    }
}

// One move on the column order (order[0..K-1]). The columnar solver's move set: swap
// dominant (70%), short reverse (15%), short block-move (15%); the direction flips only
// when both directions are searched (search_dir == COL_READ_BOTH).
static void adfgvx_order_move(int order[], int K, int *dir, int search_dir) {
    if (search_dir == COL_READ_BOTH && frand() < 0.05) { *dir = 1 - *dir; return; }
    if (K < 2) return;
    double r = frand();
    if (r < 0.70) {
        int a = rand_int(0, K), b = rand_int(0, K);
        int t = order[a]; order[a] = order[b]; order[b] = t;
    } else if (r < 0.85) {
        int max_blk = min(K, 8);
        int blk = rand_int(2, max_blk + 1);
        int s = rand_int(0, K - blk + 1);
        for (int a = s, b = s + blk - 1; a < b; a++, b--) {
            int t = order[a]; order[a] = order[b]; order[b] = t;
        }
    } else {
        int max_blk = min(K, 8);
        int blk = rand_int(1, max_blk + 1);
        int s = rand_int(0, K - blk + 1);
        int d = rand_int(0, K - blk + 1);
        if (d == s) return;
        int tmp[8];
        for (int a = 0; a < blk; a++) tmp[a] = order[s + a];
        if (d < s) { for (int a = s - 1; a >= d; a--) order[a + blk] = order[a]; }
        else       { for (int a = s + blk; a < d + blk; a++) order[a - blk] = order[a]; }
        for (int a = 0; a < blk; a++) order[d + a] = tmp[a];
    }
}

// --- model hooks ------------------------------------------------------------------

static int adfgvx_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const AdfgvxScratch *a = (const AdfgvxScratch *) ctx->model_scratch;
    int n = 0;
    for (int K = a->klo; K <= a->khi && n < cap; K++) {   // one config per column count
        out[n].period = K; out[n].j = 0; out[n].k = 0;
        out[n].aux[0] = 0; out[n].aux[1] = 0;
        n++;
    }
    return n;
}

static void adfgvx_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const AdfgvxScratch *a = (const AdfgvxScratch *) ctx->model_scratch;
    int gs = a->grid_size, K = cc->period;
    for (int i = 0; i < gs; i++) st->key[i] = i;          // random square
    shuffle(st->key, gs);
    int *order = &st->key[gs];                            // random column order
    for (int c = 0; c < K; c++) order[c] = c;
    shuffle(order, K);
    st->aux[0] = K;
    st->aux[1] = (ctx->cfg->read_direction == COL_READ_BOTH)
                     ? rand_int(0, 2) : ctx->cfg->read_direction;
    st->aux[2] = a->side;
    st->aux[3] = gs;
    st->key_len = gs + K;
}

static void adfgvx_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                           SolverState *st, bool *force_primary) {
    (void) cc; (void) force_primary;
    int gs = st->aux[3], K = st->aux[0], side = st->aux[2];
    if (frand() < 0.75) {
        adfgvx_square_move(st->key, side, gs);            // square move
    } else {
        int dir = st->aux[1];
        adfgvx_order_move(&st->key[gs], K, &dir, ctx->cfg->read_direction);
        st->aux[1] = dir;
    }
}

static void adfgvx_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int gs = src->aux[3], K = src->aux[0];
    for (int i = 0; i < gs + K; i++) dst->key[i] = src->key[i];
    for (int i = 0; i < 4; i++) dst->aux[i] = src->aux[i];
    dst->key_len = src->key_len;
}

static void adfgvx_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    const AdfgvxScratch *a = (const AdfgvxScratch *) ctx->model_scratch;
    int side = st->aux[2], K = st->aux[0], dir = st->aux[1];
    int *order = &st->key[st->aux[3]];                    // column order at offset grid_size
    adfgvx_decrypt(a->coords, a->len2, st->key, side, K, order, dir, out);
    // Structural IoC reward: depends only on the column order, not the square -> a
    // gradient for the column-order search that is flat in the n-gram landscape early on.
    *score_adjust = ADFGVX_IOC_WEIGHT * index_of_coincidence(out, a->n);
}

// Render the recovered square as an indented box of symbols.
static void adfgvx_print_square(const int sq[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(sq[r * side + c]));
        printf("\n");
    }
}

// Reconstruct the ciphertext label string from the coordinate stream (coord c -> labels[c]).
static void adfgvx_print_cipher_labels(const int coords[], int len2, int side) {
    const char *labels = adfgvx_labels(side);
    for (int i = 0; i < len2; i++) putchar(labels[coords[i]]);
}

static void adfgvx_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) decrypted;
    const AdfgvxScratch *a = (const AdfgvxScratch *) ctx->model_scratch;
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  K=%d dir=%s score=%.4f ioc=%.4f  [%.1fs, %d restarts]\n",
        st->aux[0], st->aux[1] == COL_READ_BT ? "bt" : "tb", score,
        index_of_coincidence(decrypted, a->n), elapsed, stats->n_restarts);
    adfgvx_print_square(st->key, a->side);
    fflush(stdout);
}

static void adfgvx_report(const SolverCtx *ctx, const SolverConfig *cc,
                          const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const AdfgvxScratch *a = (const AdfgvxScratch *) ctx->model_scratch;
    int n = a->n, side = a->side, gs = a->grid_size;
    int K = st->aux[0], dir = st->aux[1];
    const int *order = &st->key[gs];

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char sqstr[SQUARE_MAX_GRID + 1];
    for (int i = 0; i < gs; i++) sqstr[i] = index_to_char(st->key[i]);
    sqstr[gs] = '\0';

    const char *type_name = (cfg->cipher_type == ADFGVX) ? "adfgvx" : "adfgx";
    printf("\nResult Score: %.2f | Words: %d | K=%d dir=%s | square=%s\n",
        score, n_words_found, K, dir == COL_READ_BT ? "bt" : "tb", sqstr);

    adfgvx_print_cipher_labels(a->coords, a->len2, side);
    printf("\n");
    print_text(decrypted, n);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, n);

    printf("\nrecovered %dx%d square (row major):\n", side, side);
    adfgvx_print_square(st->key, side);
    printf("column order (K=%d, dir=%s):", K, dir == COL_READ_BT ? "bt" : "tb");
    for (int c = 0; c < K; c++) printf(" %d", order[c]);
    printf("\n");

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = K;                  // report the column count here
        vec_copy(decrypted, ctx->result->decrypted, n);
        ctx->result->decrypted_len = n;
    }

    // One-liner summary: >>> score, [words,] type, K=, square=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, K=%d, square=%s, ", score, n_words_found, cfg->cipher_type, K, sqstr);
    else
        printf(">>> %.2f, %d, K=%d, square=%s, ", score, cfg->cipher_type, K, sqstr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    (void) type_name;
    adfgvx_print_cipher_labels(a->coords, a->len2, side);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

static const CipherModel ADFGVX_MODEL = {
    .name = "adfgvx", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = adfgvx_enumerate, .key_len = NULL,
    .seed = adfgvx_seed, .perturb = adfgvx_perturb, .copy_state = adfgvx_copy,
    .decrypt = adfgvx_decrypt_hook, .report = adfgvx_report,
    .report_verbose = adfgvx_report_verbose,
};

void solve_adfgvx(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    int side = (cfg->cipher_type == ADFGVX) ? ADFGVX_SIDE : ADFGX_SIDE;
    int gs = side * side;

    // The alphabet must match the square (forced in main before load_ngrams).
    if (g_alpha != gs) {
        printf("\n\nERROR: ADFG%sX needs a %d-symbol alphabet (got %d). "
               "Run -type %s so the alphabet is forced.\n\n",
               side == 6 ? "V" : "", gs, g_alpha, side == 6 ? "adfgvx" : "adfgx");
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for an ADFGVX solve.\n\n");
        return;
    }
    if (cipher_len % 2 != 0) {
        printf("\n\nERROR: ADFGVX ciphertext length must be even (it is 2x the plaintext "
               "length); got %d.\n\n", cipher_len);
        return;
    }

    // Map the ciphertext label characters to coordinates 0..side-1. Build the inverse of
    // the label table over the active alphabet, then convert + validate each symbol.
    const char *labels = adfgvx_labels(side);
    int idx_to_coord[MAX_ALPHABET_SIZE];
    for (int i = 0; i < MAX_ALPHABET_SIZE; i++) idx_to_coord[i] = -1;
    for (int p = 0; p < side; p++) {
        int li = g_char_to_idx[(unsigned char) labels[p]];
        if (li >= 0 && li < MAX_ALPHABET_SIZE) idx_to_coord[li] = p;
    }
    for (int i = 0; i < cipher_len; i++) {
        int ci = cipher_indices[i];
        if (ci < 0 || ci >= MAX_ALPHABET_SIZE || idx_to_coord[ci] < 0) {
            printf("\n\nERROR: ciphertext has a non-label symbol at position %d; ADFG%sX "
                   "ciphertext must be over {%s}.\n\n", i, side == 6 ? "V" : "", labels);
            return;
        }
        g_adfgvx_coords[i] = idx_to_coord[ci];
    }

    int n = cipher_len / 2;

    // Column-count (K) sweep range, clamped to [2, len2/2] and the array bound.
    int klo = cfg->min_cols, khi = cfg->max_cols;
    int cap = cipher_len / 2;
    if (cap > MAX_COLS) cap = MAX_COLS;
    if (cap < 2) cap = 2;
    if (klo < 2) klo = 2;
    if (khi > cap) khi = cap;
    if (klo > khi) klo = khi;

    AdfgvxScratch scratch;
    scratch.side = side;
    scratch.grid_size = gs;
    scratch.len2 = cipher_len;
    scratch.n = n;
    scratch.klo = klo;
    scratch.khi = khi;
    scratch.coords = g_adfgvx_coords;

    if (cfg->verbose)
        printf("\nadfgvx: %d ciphertext labels (%d plaintext), side %d, %d-symbol alphabet, "
               "K in [%d..%d]\n", cipher_len, n, side, gs, klo, khi);

    // Cribs are over plaintext positions, but the crib machinery aligned them to the 2N
    // ciphertext; they do not correspond to plaintext positions, so ADFGVX ignores cribs.
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, n, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    run_solver(&ADFGVX_MODEL, &ctx);
}
