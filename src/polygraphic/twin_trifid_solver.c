#include "twin_trifid_solver.h"
#include "trifid_solver.h"      // trifid_estimate_periods (columnar-IoC estimator, reused)
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Twin Trifid solver (TYPE twin-trifid)
// =====================================================================
//
// The Trifid analogue of Twin Bifid: TWO Trifid messages enciphered under the SAME keyed
// 3x3x3 cube (27 symbols, A..Z + '+') but with DIFFERENT periods; the two plaintexts share a
// common phrase (an ACA >= 16-letter repeat). The shared cube decrypts BOTH messages, so the
// search state is a SINGLE cube (a permutation of 0..26 in st->key) -- exactly Trifid's state
// and move set -- but the objective scores the CONCATENATED decrypt of both ciphertexts.
// Judging every cube move against ~2x the text gives far more n-gram signal than a lone
// Trifid, so recovery works from the short ACA lengths (100-150 letters each).
//
// Structurally identical to the Twin Bifid solver (one perturbed key, score_adjust 0, the
// Tri-Square "scoring length = n1+n2" wiring), swapping Bifid's 5x5 square for Trifid's cube:
// the seed shuffles the cube, the move set is Trifid's (cell swap / plane swap along a random
// axis / axis reflection), and the period estimator is trifid_estimate_periods. The two
// periods are independent unknowns: bifid-style top-K per message, crossed; -period pins
// message 1, -period2 pins message 2. The ACA recommends periods "not both divisible by 3",
// a con-construction guideline -- the sweep stays general (the n-gram score picks).
//
// The second ciphertext is supplied via cfg->twincipher_str (-cipher2 on the CLI; tests set
// it directly) and decoded with the same decode_cipher() the primary uses (27-symbol alphabet).

#define TWIN_TRIFID_MAX_PAIRS 64
#define TWIN_TRIFID_MAX_PER   8    // candidate periods considered per message before crossing

typedef struct {
    int side;                     // cube side (3 for the 27-symbol default)
    int side2;                    // side*side
    int cube_size;                // side^3 == g_alpha
    const int *c1; int n1;        // message 1 ciphertext + length
    const int *c2; int n2;        // message 2 ciphertext + length
    int n_pairs;                  // number of (p1,p2) period pairs to anneal
    int p1[TWIN_TRIFID_MAX_PAIRS];
    int p2[TWIN_TRIFID_MAX_PAIRS];
} TwinTrifidScratch;

// Setup-phase statics (written once on the main thread, read-only during the search): the
// decoded second ciphertext and the combined (c1 ++ c2) buffer the engine sees as ctx->cipher.
static int g_twin_trifid_c2[MAX_CIPHER_LENGTH];
static int g_twin_trifid_combined[MAX_CIPHER_LENGTH];

// Integer cube root: returns s if n == s^3 for some s >= 1, else -1.
static int exact_icbrt(int n) {
    for (int s = 1; s * s * s <= n; s++)
        if (s * s * s == n) return s;
    return -1;
}

// One config per (p1,p2) period pair; period carries message 1's block size, aux[0] message 2's.
static int twin_trifid_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const TwinTrifidScratch *t = (const TwinTrifidScratch *) ctx->model_scratch;
    int n = t->n_pairs;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = t->p1[i];
        out[i].aux[0] = t->p2[i];
        out[i].j = 0; out[i].k = 0; out[i].aux[1] = 0;
    }
    return n;
}

// Seed: one uniformly random cube (Fisher-Yates shuffle of 0..cube_size-1).
static void twin_trifid_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    const TwinTrifidScratch *t = (const TwinTrifidScratch *) ctx->model_scratch;
    int n = t->cube_size;
    for (int i = 0; i < n; i++) st->key[i] = i;
    for (int i = n - 1; i > 0; i--) {
        int j = rand_int(0, i + 1);
        int tmp = st->key[i]; st->key[i] = st->key[j]; st->key[j] = tmp;
    }
    st->key_len = n;
}

// Swap the two planes a, b along coordinate axis `ax` (0 layer / 1 row / 2 column).
static void twin_trifid_swap_plane(int key[], int side, int side2, int ax, int a, int b) {
    if (a == b) return;
    for (int u = 0; u < side; u++)
        for (int v = 0; v < side; v++) {
            int pa, pb;
            if (ax == 0)      { pa = a * side2 + u * side + v;  pb = b * side2 + u * side + v; }
            else if (ax == 1) { pa = u * side2 + a * side + v;  pb = u * side2 + b * side + v; }
            else              { pa = u * side2 + v * side + a;  pb = u * side2 + v * side + b; }
            int tmp = key[pa]; key[pa] = key[pb]; key[pb] = tmp;
        }
}

// Reflect the cube along coordinate axis `ax` (mirror coordinate w <-> side-1-w).
static void twin_trifid_reflect(int key[], int side, int side2, int ax) {
    for (int w = 0; w < side / 2; w++)
        twin_trifid_swap_plane(key, side, side2, ax, w, side - 1 - w);
}

// Neighbour move on the shared cube (Trifid's move set): 82% swap two cells; 12% swap two
// planes along a random axis; 6% reflect along a random axis.
static void twin_trifid_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                                SolverState *st, bool *force_primary) {
    const TwinTrifidScratch *t = (const TwinTrifidScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int s = t->side, s2 = t->side2, n = t->cube_size;
    double r = frand();
    if (r < 0.82) {                              // swap two cells
        int a = rand_int(0, n), c = rand_int(0, n);
        int tmp = st->key[a]; st->key[a] = st->key[c]; st->key[c] = tmp;
    } else if (r < 0.94) {                       // swap two planes along a random axis
        int ax = rand_int(0, 3);
        int a = rand_int(0, s), b = rand_int(0, s);
        twin_trifid_swap_plane(st->key, s, s2, ax, a, b);
    } else {                                     // reflect along a random axis
        int ax = rand_int(0, 3);
        twin_trifid_reflect(st->key, s, s2, ax);
    }
}

static void twin_trifid_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < TRIFID_MAX_CELLS; i++) dst->key[i] = src->key[i];
}

static void twin_trifid_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                     SolverState *st, int *out, double *score_adjust) {
    const TwinTrifidScratch *t = (const TwinTrifidScratch *) ctx->model_scratch;
    twin_trifid_decrypt(t->c1, t->n1, t->c2, t->n2, st->key, t->side,
                        cc->period, cc->aux[0], out);
    *score_adjust = 0.0;
}

// Render the cube as `side` indented layers of letters.
static void twin_trifid_print_cube(const int cube[], int side) {
    int side2 = side * side;
    for (int l = 0; l < side; l++) {
        printf("    layer %d:\n", l + 1);
        for (int r = 0; r < side; r++) {
            printf("      ");
            for (int c = 0; c < side; c++) printf("%c ", index_to_char(cube[l * side2 + r * side + c]));
            printf("\n");
        }
    }
}

static void twin_trifid_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const TwinTrifidScratch *t = (const TwinTrifidScratch *) ctx->model_scratch;
    printf("\n  periods (%d, %d), shared cube:\n", cc->period, cc->aux[0]);
    twin_trifid_print_cube(st->key, t->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "twin-trifid");
}

static void twin_trifid_report(const SolverCtx *ctx, const SolverConfig *cc,
                               const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const TwinTrifidScratch *t = (const TwinTrifidScratch *) ctx->model_scratch;
    int side = t->side, csz = t->cube_size, n1 = t->n1, n2 = t->n2, len = n1 + n2;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The recovered cube, read cell-major (layer/row/column order).
    char cubestr[TRIFID_MAX_CELLS + 1];
    for (int i = 0; i < csz; i++) cubestr[i] = index_to_char(st->key[i]);
    cubestr[csz] = '\0';

    printf("\nResult Score: %.2f | Words: %d | period1=%d | period2=%d | cube=%s\n",
        score, n_words_found, cc->period, cc->aux[0], cubestr);

    printf("cipher 1 (period %d): ", cc->period);
    print_cipher((int *) t->c1, n1, NULL);
    printf("\ncipher 2 (period %d): ", cc->aux[0]);
    print_cipher((int *) t->c2, n2, NULL);
    printf("\nplaintext 1: ");
    print_text(decrypted, n1);
    printf("\nplaintext 2: ");
    print_text(decrypted + n1, n2);
    printf("\n");
    print_solution_check(decrypted, len);

    printf("\nrecovered %dx%dx%d cube (cell major):\n", side, side, side);
    twin_trifid_print_cube(st->key, side);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = cc->period;   // message-1 period (message-2 = cc->aux[0])
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, period1=, period2=, cube=, file,
    // CIPHER1CIPHER2, PLAINTEXT1PLAINTEXT2 (concatenated plaintext is the final field).
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, period1=%d, period2=%d, cube=%s, ",
            score, n_words_found, cfg->cipher_type, cc->period, cc->aux[0], cubestr);
    else
        printf(">>> %.2f, %d, period1=%d, period2=%d, cube=%s, ",
            score, cfg->cipher_type, cc->period, cc->aux[0], cubestr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel TWIN_TRIFID_MODEL = {
    .name = "twin-trifid", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = twin_trifid_enumerate, .key_len = NULL,
    .seed = twin_trifid_seed, .perturb = twin_trifid_perturb, .copy_state = twin_trifid_copy,
    .decrypt = twin_trifid_decrypt_hook, .report = twin_trifid_report,
    .report_verbose = twin_trifid_report_verbose,
};

// Rank a single message's candidate periods: pinned value, or the estimator's top-K.
static int twin_trifid_periods_for(int cipher[], int len, bool pinned, int pinned_p,
                                   const ColossusConfig *cfg, int out[]) {
    if (pinned) { out[0] = pinned_p; return 1; }
    int max_p = (cfg->max_period > 0) ? cfg->max_period : 20;
    if (max_p > len / 2) max_p = len / 2;
    int n_want = cfg->n_periods;
    if (n_want > TWIN_TRIFID_MAX_PER) n_want = TWIN_TRIFID_MAX_PER;
    int n = trifid_estimate_periods(cipher, len, 2, max_p, n_want, out, cfg->verbose);
    if (n < 1) { out[0] = 2; n = 1; }
    return n;
}

void solve_twin_trifid(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    // Twin Trifid needs a perfect-cube (27) alphabet; the binary forces 27.
    int side = exact_icbrt(g_alpha);
    if (side < 2 || side > TRIFID_MAX_SIDE || side * side * side != g_alpha) {
        printf("\n\nERROR: Twin Trifid needs a perfect-cube alphabet (e.g. 27; got %d). "
               "Use the default -type twin-trifid (27 symbols: A..Z + '%c').\n\n",
               g_alpha, TRIFID_EXTRA_CHAR);
        return;
    }

    if (!cfg->twincipher_str || cfg->twincipher_str[0] == '\0') {
        printf("\n\nERROR: Twin Trifid needs a SECOND ciphertext -- supply it with "
               "-cipher2 <file> (the two messages share one keyed cube).\n\n");
        return;
    }

    int n1 = cipher_len;
    SymbolTable symtab2;
    int n2 = decode_cipher(cfg->twincipher_str, cfg, g_twin_trifid_c2, &symtab2);

    if (n1 < 4 || n2 < 4) {
        printf("\n\nERROR: both Twin Trifid ciphertexts must be at least 4 symbols "
               "(got %d and %d).\n\n", n1, n2);
        return;
    }
    if (n1 + n2 > MAX_CIPHER_LENGTH) {
        printf("\n\nERROR: combined Twin Trifid ciphertext length %d exceeds the %d limit.\n\n",
               n1 + n2, MAX_CIPHER_LENGTH);
        return;
    }
    for (int i = 0; i < n1; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext 1 has a symbol outside the 27-symbol cube "
                   "alphabet (A..Z + '%c') at position %d.\n\n", TRIFID_EXTRA_CHAR, i);
            return;
        }
    for (int i = 0; i < n2; i++)
        if (g_twin_trifid_c2[i] < 0 || g_twin_trifid_c2[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext 2 has a symbol outside the 27-symbol cube "
                   "alphabet (A..Z + '%c') at position %d.\n\n", TRIFID_EXTRA_CHAR, i);
            return;
        }

    TwinTrifidScratch scratch;
    scratch.side = side;
    scratch.side2 = side * side;
    scratch.cube_size = g_alpha;
    scratch.c1 = cipher_indices; scratch.n1 = n1;
    scratch.c2 = g_twin_trifid_c2; scratch.n2 = n2;

    int cand1[TWIN_TRIFID_MAX_PER], cand2[TWIN_TRIFID_MAX_PER];
    int nc1 = twin_trifid_periods_for(cipher_indices, n1, cfg->period_present, cfg->period, cfg, cand1);
    int nc2 = twin_trifid_periods_for(g_twin_trifid_c2, n2, cfg->period2_present, cfg->period2, cfg, cand2);
    int np = 0;
    for (int i = 0; i < nc1 && np < TWIN_TRIFID_MAX_PAIRS; i++)
        for (int j = 0; j < nc2 && np < TWIN_TRIFID_MAX_PAIRS; j++) {
            scratch.p1[np] = cand1[i];
            scratch.p2[np] = cand2[j];
            np++;
        }
    scratch.n_pairs = np;

    if (cfg->verbose)
        printf("\ntwin-trifid: msg1 %d symbols, msg2 %d symbols, %d-symbol alphabet %s, "
               "%d candidate period pair(s)\n", n1, n2, g_alpha, g_idx_to_char_arr, np);

    for (int i = 0; i < n1; i++) g_twin_trifid_combined[i] = cipher_indices[i];
    for (int i = 0; i < n2; i++) g_twin_trifid_combined[n1 + i] = g_twin_trifid_c2[i];

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        g_twin_trifid_combined, n1 + n2, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    run_solver(&TWIN_TRIFID_MODEL, &ctx);
}
