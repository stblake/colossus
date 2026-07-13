#include "nihilist_sub_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Nihilist Substitution solver (TYPE nihilist-sub / -nc / -m100)
// =====================================================================
//
// Nihilist Substitution is a periodic ADDITIVE cipher over a keyed Polybius square: each
// plaintext letter -> its 2-digit coordinate number, a periodic additive key (its own
// coordinate numbers) is added per position (one of the NIH_ADD_* conventions), and the
// ciphertext is the resulting stream of decimal numbers. Breaking it is a COUPLED problem --
// recover the square AND the additive key -- structurally the twin of ADFGVX.
//
// The state carries both, packed in st->key like ADFGVX:
//   key[0 .. grid_size-1]              the keyed square (a permutation of 0..grid_size-1)
//   key[grid_size .. grid_size+p-1]    the additive key as coordinate cells (0..grid_size-1)
//   aux[0]=p  aux[1]=conv  aux[2]=side  aux[3]=grid_size
// One engine config is enumerated per candidate period p.
//
// The decoupling trick (the crux, mirroring ADFGVX's IoC reward): pt_num = cipher_num (-)
// key_num is INDEPENDENT of the square (the square only maps the resulting coordinate to a
// letter). So the fraction of positions that decrypt to a LEGAL coordinate (both digits in
// 1..side) depends ONLY on the additive key, not the square. Folding that "validity" fraction
// (NIH_VALID_WEIGHT * n_valid/n) into score_adjust gives the additive-key search a gradient
// flat in the square dimension: the climb drives validity to 1 (locking the key), then the
// n-gram score recovers the square as a monoalphabetic map over the decrypted coordinates.
// Once validity saturates the reward is a constant offset, so the n-gram score decides the
// square -- exactly how ADFGVX's IoC term decouples its column order from its square.
//
// The square move set is Bifid/Playfair's (cell-swap dominant + row/column swaps and
// reflections); the additive-key move redraws one column's cell (occasionally swaps two).
// The solver ASSUMES the standard fixed labels (1..side); a keyed-label cipher folds its label
// permutation into the recovered square (it is not separately identifiable ciphertext-only),
// so the same model cracks it as the equivalent relabeled square.

#define NIH_VALID_WEIGHT 14.0    // weight of the square-independent validity reward (tuned)
#define NIH_MAX_PERIODS  64

typedef struct {
    int  side;                       // 5 (25-letter square)
    int  grid_size;                  // side*side == g_alpha
    int  conv;                       // NIH_ADD_CARRY / _NOCARRY / _MOD100
    int  n;                          // plaintext length == number of cipher numbers
    int  n_periods;                  // number of candidate periods
    int  periods[NIH_MAX_PERIODS];   // the candidate periods (config order)
    int  rowlbl[NIHILIST_SUB_MAX_SIDE];
    int  collbl[NIHILIST_SUB_MAX_SIDE];
    int *values;                     // ciphertext numbers (length n)
} NihilistSubScratch;

// Ciphertext number buffer (single-threaded), filled by solve_nihilist_sub.
static int g_nih_values[MAX_CIPHER_LENGTH];

// cipher_type -> addition convention.
static int nih_conv_of_type(int t) {
    if (t == NIHILIST_SUB_NC)   return NIH_ADD_NOCARRY;
    if (t == NIHILIST_SUB_M100) return NIH_ADD_MOD100;
    return NIH_ADD_CARRY;
}

// --- period estimation --------------------------------------------------------
//
// Rank trial periods by the columnar Index of Coincidence over the ciphertext NUMBERS (used
// as symbols; values <= 110 fit a 128-wide table). At the true period each column shares one
// additive, so its number distribution is a fixed shift/relabel of the plaintext coordinate
// distribution -- English, hence elevated IoC; a wrong period mixes additives and flattens it.
// (mean_ioc can't be reused: it tables on [0..MAX_ALPHABET_SIZE), and the numbers exceed that.)
static double nih_columnar_ioc(const int vals[], int n, int period) {
    double sum = 0.0;
    for (int k = 0; k < period; k++) {
        int freq[128];
        for (int v = 0; v < 128; v++) freq[v] = 0;
        int cnt = 0;
        for (int i = k; i < n; i += period) { freq[vals[i] & 127]++; cnt++; }
        if (cnt < 2) continue;
        long num = 0;
        for (int v = 0; v < 128; v++) num += (long) freq[v] * (freq[v] - 1);
        sum += (double) num / ((double) cnt * (cnt - 1));
    }
    return sum / period;
}

int nihilist_sub_estimate_periods(const int values[], int n, int min_p, int max_p,
                                  int n_want, int out[], bool verbose) {
    static double ioc[MAX_CIPHER_LENGTH];
    if (max_p > n / 2) max_p = n / 2;
    if (max_p < min_p) max_p = min_p;
    if (min_p < 1) min_p = 1;
    if (n_want < 1) n_want = 1;
    if (n_want > NIH_MAX_PERIODS) n_want = NIH_MAX_PERIODS;

    for (int p = min_p; p <= max_p; p++)
        ioc[p] = nih_columnar_ioc(values, n, p);

    if (verbose) {
        printf("\nNihilist Substitution period estimate (columnar IoC over ciphertext numbers):\n  period\tIoC\n");
        for (int p = min_p; p <= max_p; p++) printf("  %d\t%.4f\n", p, ioc[p]);
    }

    int cnt = 0;
    for (; cnt < n_want; cnt++) {
        int best_p = -1; double best = -1.0;
        for (int p = min_p; p <= max_p; p++)
            if (ioc[p] >= 0.0 && ioc[p] > best) { best = ioc[p]; best_p = p; }
        if (best_p < 0) break;
        out[cnt] = best_p;
        ioc[best_p] = -1.0;
    }
    if (verbose) {
        printf("  -> annealing periods:");
        for (int i = 0; i < cnt; i++) printf(" %d", out[i]);
        printf("\n");
    }
    return cnt;
}

// --- square neighbour move (identical to Bifid/Playfair/ADFGVX) ---------------
static void nih_square_move(int sq[], int s, int n) {
    double r = frand();
    if (r < 0.80) {                              // swap two cells
        int a = rand_int(0, n), c = rand_int(0, n);
        int t = sq[a]; sq[a] = sq[c]; sq[c] = t;
    } else if (r < 0.88) {                        // swap two rows
        int r1 = rand_int(0, s), r2 = rand_int(0, s);
        for (int c = 0; c < s; c++) {
            int t = sq[r1 * s + c]; sq[r1 * s + c] = sq[r2 * s + c]; sq[r2 * s + c] = t;
        }
    } else if (r < 0.96) {                        // swap two columns
        int c1 = rand_int(0, s), c2 = rand_int(0, s);
        for (int rr = 0; rr < s; rr++) {
            int t = sq[rr * s + c1]; sq[rr * s + c1] = sq[rr * s + c2]; sq[rr * s + c2] = t;
        }
    } else if (r < 0.98) {                        // rotate 180
        for (int i = 0, j = n - 1; i < j; i++, j--) { int t = sq[i]; sq[i] = sq[j]; sq[j] = t; }
    } else if (r < 0.99) {                        // flip rows
        for (int r1 = 0, r2 = s - 1; r1 < r2; r1++, r2--)
            for (int c = 0; c < s; c++) {
                int t = sq[r1 * s + c]; sq[r1 * s + c] = sq[r2 * s + c]; sq[r2 * s + c] = t;
            }
    } else {                                      // flip columns
        for (int c1 = 0, c2 = s - 1; c1 < c2; c1++, c2--)
            for (int rr = 0; rr < s; rr++) {
                int t = sq[rr * s + c1]; sq[rr * s + c1] = sq[rr * s + c2]; sq[rr * s + c2] = t;
            }
    }
}

// --- model hooks --------------------------------------------------------------

static int nih_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const NihilistSubScratch *a = (const NihilistSubScratch *) ctx->model_scratch;
    int n = a->n_periods;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = a->periods[i];
        out[i].j = 0; out[i].k = 0; out[i].aux[0] = 0; out[i].aux[1] = 0;
    }
    return n;
}

static void nih_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const NihilistSubScratch *a = (const NihilistSubScratch *) ctx->model_scratch;
    int gs = a->grid_size, p = cc->period;
    for (int i = 0; i < gs; i++) st->key[i] = i;          // random square
    shuffle(st->key, gs);
    for (int j = 0; j < p; j++) st->key[gs + j] = rand_int(0, gs);   // random additive cells
    st->aux[0] = p;
    st->aux[1] = a->conv;
    st->aux[2] = a->side;
    st->aux[3] = gs;
    st->key_len = gs + p;
}

static void nih_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                        SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int gs = st->aux[3], p = st->aux[0], side = st->aux[2];
    if (frand() < 0.72) {
        nih_square_move(st->key, side, gs);               // square move
    } else {
        int *key = &st->key[gs];
        if (p >= 2 && frand() < 0.15) {                   // swap two additive columns
            int a = rand_int(0, p), b = rand_int(0, p);
            int t = key[a]; key[a] = key[b]; key[b] = t;
        } else {                                          // redraw one additive column's cell
            key[rand_int(0, p)] = rand_int(0, gs);
        }
    }
}

static void nih_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int gs = src->aux[3], p = src->aux[0];
    for (int i = 0; i < gs + p; i++) dst->key[i] = src->key[i];
    for (int i = 0; i < 4; i++) dst->aux[i] = src->aux[i];
    dst->key_len = src->key_len;
}

static void nih_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                             SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    const NihilistSubScratch *a = (const NihilistSubScratch *) ctx->model_scratch;
    int p = st->aux[0], conv = st->aux[1], side = st->aux[2], gs = st->aux[3];
    int *key_cells = &st->key[gs];
    int n_valid = nihilist_sub_decrypt(a->values, a->n, st->key, a->rowlbl, a->collbl,
                                       side, key_cells, p, conv, out);
    // Square-independent validity reward: drives the additive key to all-legal, decoupling it
    // from the square (which the n-gram score then recovers).
    *score_adjust = NIH_VALID_WEIGHT * (double) n_valid / (double) a->n;
}

// --- reporting ----------------------------------------------------------------

static const char *nih_type_name(int t) {
    if (t == NIHILIST_SUB_NC)   return "nihilist-sub-nc";
    if (t == NIHILIST_SUB_M100) return "nihilist-sub-m100";
    return "nihilist-sub";
}

static void nih_print_square(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void nih_print_numbers(const int vals[], int n) {
    for (int i = 0; i < n; i++) printf("%s%d", i ? " " : "", vals[i]);
}

// The additive key cell -> its fixed-label coordinate number (1..side labels).
static int nih_cell_num_fixed(int cell, int side) {
    return (cell / side + 1) * 10 + (cell % side + 1);
}

static void nih_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const NihilistSubScratch *a = (const NihilistSubScratch *) ctx->model_scratch;
    int p = st->aux[0];
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  period %d, score=%.4f  [%.1fs, %d restarts]\n", p, score, elapsed, stats->n_restarts);
    nih_print_square(st->key, a->side);
    (void) cc; (void) decrypted;
    fflush(stdout);
}

static void nih_report(const SolverCtx *ctx, const SolverConfig *cc,
                       const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const NihilistSubScratch *a = (const NihilistSubScratch *) ctx->model_scratch;
    int n = a->n, side = a->side, gs = a->grid_size, p = st->aux[0];
    const int *key_cells = &st->key[gs];

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The recovered square (row major) and the additive key as numbers + reconstructed keyword.
    char gridstr[NIHILIST_SUB_MAX_GRID + 1];
    for (int i = 0; i < gs; i++) gridstr[i] = index_to_char(st->key[i]);
    gridstr[gs] = '\0';

    char keyword[NIH_MAX_PERIODS + 1];
    for (int j = 0; j < p && j < NIH_MAX_PERIODS; j++) keyword[j] = index_to_char(st->key[key_cells[j]]);
    keyword[(p < NIH_MAX_PERIODS) ? p : NIH_MAX_PERIODS] = '\0';

    printf("\nResult Score: %.2f | Words: %d | period=%d | square=%s | key=%s\n",
        score, n_words_found, p, gridstr, keyword);

    nih_print_numbers(a->values, n);
    printf("\n");
    print_text(decrypted, n);
    printf("\n");
    print_solution_check(decrypted, n);
    print_spaces_line(g_spaces_table, decrypted, n);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered %dx%d square (row major):\n", side, side);
    nih_print_square(st->key, side);
    printf("additive key (period %d): keyword=%s, numbers=", p, keyword);
    for (int j = 0; j < p; j++) printf("%s%d", j ? " " : "", nih_cell_num_fixed(key_cells[j], side));
    printf("\n");

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = p;            // report the recovered period here
        vec_copy(decrypted, ctx->result->decrypted, n);
        ctx->result->decrypted_len = n;
    }

    // One-liner summary: >>> score, [words,] type, period=, square=, key=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, period=%d, square=%s, key=%s, ",
            score, n_words_found, cfg->cipher_type, p, gridstr, keyword);
    else
        printf(">>> %.2f, %d, period=%d, square=%s, key=%s, ",
            score, cfg->cipher_type, p, gridstr, keyword);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    nih_print_numbers(a->values, n);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
    (void) cc;
}

static const CipherModel NIHILIST_SUB_MODEL = {
    .name = "nihilist-sub", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = nih_enumerate, .key_len = NULL,
    .seed = nih_seed, .perturb = nih_perturb, .copy_state = nih_copy,
    .decrypt = nih_decrypt_hook, .report = nih_report,
    .report_verbose = nih_report_verbose,
};

// Parse a stream of decimal numbers (any non-digit run is a separator) into out[]; returns n.
static int nih_parse_numbers(const char *s, int out[], int cap) {
    int n = 0, i = 0;
    while (s[i] && n < cap) {
        if (s[i] >= '0' && s[i] <= '9') {
            int v = 0;
            while (s[i] >= '0' && s[i] <= '9') { v = v * 10 + (s[i] - '0'); i++; }
            out[n++] = v;
        } else {
            i++;
        }
    }
    return n;
}

void solve_nihilist_sub(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_len;     // the per-character decode length is irrelevant; we parse numbers

    int side = NIHILIST_SUB_SIDE, gs = side * side;
    if (g_alpha != gs) {
        printf("\n\nERROR: Nihilist Substitution needs a %d-letter alphabet (got %d). "
               "Run -type nihilist-sub so the alphabet is forced (J->I).\n\n", gs, g_alpha);
        return;
    }

    // Parse the ciphertext numbers from the raw string (space/comma/any-separator delimited).
    int n = nih_parse_numbers(ciphertext_str, g_nih_values, MAX_CIPHER_LENGTH);
    if (n < 4) {
        printf("\n\nERROR: parsed only %d ciphertext numbers; need >= 4. The Nihilist "
               "Substitution ciphertext must be space/comma-separated decimal numbers.\n\n", n);
        return;
    }
    int conv = nih_conv_of_type(cfg->cipher_type);

    NihilistSubScratch scratch;
    scratch.side = side;
    scratch.grid_size = gs;
    scratch.conv = conv;
    scratch.n = n;
    scratch.values = g_nih_values;
    nihilist_sub_fixed_labels(scratch.rowlbl, scratch.collbl, side);

    // Candidate periods: pinned, or the estimator's top-K over [1 .. max_period].
    if (cfg->period_present) {
        scratch.periods[0] = cfg->period;
        scratch.n_periods = 1;
        if (cfg->verbose) printf("\nnihilist-sub: period pinned to %d\n", cfg->period);
    } else {
        int max_p = (cfg->max_period > 0) ? cfg->max_period : 15;
        if (max_p > n / 2) max_p = n / 2;
        int n_want = (cfg->n_periods > 0) ? cfg->n_periods : 5;
        scratch.n_periods = nihilist_sub_estimate_periods(g_nih_values, n, 1, max_p,
            n_want, scratch.periods, cfg->verbose);
        if (scratch.n_periods < 1) { scratch.periods[0] = 1; scratch.n_periods = 1; }
    }

    if (cfg->verbose)
        printf("\nnihilist-sub (%s): %d ciphertext numbers, %d-letter alphabet %s, %d candidate period(s)\n",
            nih_type_name(cfg->cipher_type), n, g_alpha, g_idx_to_char_arr, scratch.n_periods);

    // Cribs are over plaintext positions but were aligned to the (per-char) ciphertext;
    // they do not line up with the parsed numbers, so cribs are ignored (like ADFGVX).
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, n, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    run_solver(&NIHILIST_SUB_MODEL, &ctx);
}
