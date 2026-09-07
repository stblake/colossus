#include "twin_bifid_solver.h"
#include "bifid_solver.h"       // bifid_estimate_periods (the columnar-IoC estimator, reused)
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"       // exact_isqrt

// =====================================================================
//  Twin Bifid solver (TYPE twin-bifid)
// =====================================================================
//
// A Twin Bifid is TWO Bifid messages enciphered under the SAME keyed 5x5 Polybius square
// (J->I, so the binary forces g_alpha == 25) but with DIFFERENT periods; the two plaintexts
// share a common phrase (an ACA >= 18-letter repeat). The shared square is the whole crack:
// one square decrypts BOTH messages, so the search state is a SINGLE square (a permutation
// of 0..24 in st->key) -- exactly Bifid's state and move set -- but the objective scores the
// CONCATENATED decrypt of both ciphertexts. Judging every square move against ~2x the text
// gives far more n-gram signal per evaluation than a lone Bifid, so recovery works from the
// short ACA lengths (100-150 letters each) where a single Bifid would be near hopeless.
//
// This is the CM-Bifid joint anneal with ONE square instead of two (no anti-collapse penalty
// -- a square is a bijection -- so score_adjust stays 0), plus the Tri-Square "scoring length
// != raw ciphertext length" wiring: both ciphertexts live in the scratch, the engine's
// scoring length is n1+n2, and the decrypt hook emits the two plaintexts back-to-back. The
// single n-gram spanning the concatenation seam is negligible noise (like Tri-Square's).
//
// The two periods are INDEPENDENT unknowns. bifid_estimate_periods (columnar IoC, square-
// agnostic) ranks each message's period separately; we anneal the CROSS PRODUCT of the two
// top-K lists (one engine config per (p1,p2) pair) and the n-gram score picks the winner.
// -period pins message 1's period, -period2 pins message 2's. The ACA recommends periods
// "not both even", but that is a con-construction guideline -- the sweep stays general.
//
// The second ciphertext is supplied via cfg->twincipher_str (-cipher2 on the CLI; the in-
// process tests set it directly). It is decoded here with the same decode_cipher() the
// primary uses, so it honours the active J->I alphabet.

#define TWIN_BIFID_MAX_PAIRS 64
#define TWIN_BIFID_MAX_PER   8    // candidate periods considered per message before crossing

typedef struct {
    int side;                    // grid side (5 for the 25-letter default, up to 6)
    int grid_size;               // side*side == g_alpha
    const int *c1; int n1;       // message 1 ciphertext + length
    const int *c2; int n2;       // message 2 ciphertext + length
    int n_pairs;                 // number of (p1,p2) period pairs to anneal
    int p1[TWIN_BIFID_MAX_PAIRS]; // message-1 period per config
    int p2[TWIN_BIFID_MAX_PAIRS]; // message-2 period per config
} TwinBifidScratch;

// Setup-phase statics (written once on the main thread in solve_twin_bifid, read-only during
// the search -- so safe under -nthreads, like the per-cipher scratch): the decoded second
// ciphertext and the combined (c1 ++ c2) buffer the engine sees as ctx->cipher.
static int g_twin_bifid_c2[MAX_CIPHER_LENGTH];
static int g_twin_bifid_combined[MAX_CIPHER_LENGTH];

// One config per (p1,p2) period pair; period carries message 1's block size, aux[0] carries
// message 2's; st->key carries the single shared square.
static int twin_bifid_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const TwinBifidScratch *b = (const TwinBifidScratch *) ctx->model_scratch;
    int n = b->n_pairs;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = b->p1[i];
        out[i].aux[0] = b->p2[i];
        out[i].j = 0; out[i].k = 0; out[i].aux[1] = 0;
    }
    return n;
}

// Seed: one uniformly random square (Fisher-Yates shuffle of 0..grid_size-1).
static void twin_bifid_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    const TwinBifidScratch *b = (const TwinBifidScratch *) ctx->model_scratch;
    int n = b->grid_size;
    for (int i = 0; i < n; i++) st->key[i] = i;
    for (int i = n - 1; i > 0; i--) {
        int j = rand_int(0, i + 1);
        int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
    }
    st->key_len = n;
}

// Neighbour move on the shared square: the canonical Playfair/Bifid move set -- 80% swap two
// cells; 8% swap two rows; 8% swap two columns; 2% reverse (rotate 180); 1% flip rows; 1%
// flip columns. The larger moves jump the basins a single cell swap cannot escape.
static void twin_bifid_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                               SolverState *st, bool *force_primary) {
    const TwinBifidScratch *b = (const TwinBifidScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int s = b->side, n = b->grid_size;
    double r = frand();
    if (r < 0.80) {                              // swap two cells
        int a = rand_int(0, n), c = rand_int(0, n);
        int t = st->key[a]; st->key[a] = st->key[c]; st->key[c] = t;
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

static void twin_bifid_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < BIFID_MAX_GRID; i++) dst->key[i] = src->key[i];
}

static void twin_bifid_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                    SolverState *st, int *out, double *score_adjust) {
    const TwinBifidScratch *b = (const TwinBifidScratch *) ctx->model_scratch;
    twin_bifid_decrypt(b->c1, b->n1, b->c2, b->n2, st->key, b->side,
                       cc->period, cc->aux[0], out);
    *score_adjust = 0.0;
}

// Render the square as an indented box of letters.
static void twin_bifid_print_grid(const int grid[], int side) {
    for (int r = 0; r < side; r++) {
        printf("    ");
        for (int c = 0; c < side; c++) printf("%c ", index_to_char(grid[r * side + c]));
        printf("\n");
    }
}

static void twin_bifid_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const TwinBifidScratch *b = (const TwinBifidScratch *) ctx->model_scratch;
    printf("\n  periods (%d, %d), shared square:\n", cc->period, cc->aux[0]);
    twin_bifid_print_grid(st->key, b->side);
    report_transposition_verbose(ctx, score, decrypted, stats, "twin-bifid");
}

static void twin_bifid_report(const SolverCtx *ctx, const SolverConfig *cc,
                              const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const TwinBifidScratch *b = (const TwinBifidScratch *) ctx->model_scratch;
    int side = b->side, gsz = b->grid_size, n1 = b->n1, n2 = b->n2, len = n1 + n2;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The recovered square (row-major; unique only up to a cyclic row/column rotation, which
    // all decrypt identically -- one representative). The recovered plaintexts are unique.
    char gridstr[BIFID_MAX_GRID + 1];
    for (int i = 0; i < gsz; i++) gridstr[i] = index_to_char(st->key[i]);
    gridstr[gsz] = '\0';

    printf("\nResult Score: %.2f | Words: %d | period1=%d | period2=%d | square=%s\n",
        score, n_words_found, cc->period, cc->aux[0], gridstr);

    printf("cipher 1 (period %d): ", cc->period);
    print_cipher((int *) b->c1, n1, NULL);
    printf("\ncipher 2 (period %d): ", cc->aux[0]);
    print_cipher((int *) b->c2, n2, NULL);
    printf("\nplaintext 1: ");
    print_text(decrypted, n1);
    printf("\nplaintext 2: ");
    print_text(decrypted + n1, n2);
    printf("\n");
    print_solution_check(decrypted, len);

    printf("\nrecovered %dx%d square (row major):\n", side, side);
    twin_bifid_print_grid(st->key, side);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = cc->period;   // message-1 period (message-2 = cc->aux[0])
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, period1=, period2=, square=, file,
    // CIPHER1CIPHER2, PLAINTEXT1PLAINTEXT2 (the concatenated plaintext is the final field,
    // matching the concatenated .solution the accuracy suite compares against).
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, period1=%d, period2=%d, square=%s, ",
            score, n_words_found, cfg->cipher_type, cc->period, cc->aux[0], gridstr);
    else
        printf(">>> %.2f, %d, period1=%d, period2=%d, square=%s, ",
            score, cfg->cipher_type, cc->period, cc->aux[0], gridstr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel TWIN_BIFID_MODEL = {
    .name = "twin-bifid", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = twin_bifid_enumerate, .key_len = NULL,
    .seed = twin_bifid_seed, .perturb = twin_bifid_perturb, .copy_state = twin_bifid_copy,
    .decrypt = twin_bifid_decrypt_hook, .report = twin_bifid_report,
    .report_verbose = twin_bifid_report_verbose,
};

// Rank a single message's candidate periods: pinned value, or the estimator's top-K.
static int twin_bifid_periods_for(int cipher[], int len, bool pinned, int pinned_p,
                                  const ColossusConfig *cfg, int out[]) {
    if (pinned) { out[0] = pinned_p; return 1; }
    int max_p = (cfg->max_period > 0) ? cfg->max_period : 20;
    if (max_p > len / 2) max_p = len / 2;
    int n_want = cfg->n_periods;
    if (n_want > TWIN_BIFID_MAX_PER) n_want = TWIN_BIFID_MAX_PER;
    int n = bifid_estimate_periods(cipher, len, 2, max_p, n_want, out, cfg->verbose);
    if (n < 1) { out[0] = 2; n = 1; }
    return n;
}

void solve_twin_bifid(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;
    // Cribs do not map cleanly onto the concatenated two-message plaintext, so Twin Bifid
    // ignores them (like Bifid / CM Bifid / Tri-Square).
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    // Twin Bifid needs a perfect-square (25 or 36) alphabet; the binary forces 25.
    int side = exact_isqrt(g_alpha);
    if (side < 2 || side > BIFID_MAX_SIDE || side * side != g_alpha) {
        printf("\n\nERROR: Twin Bifid needs a perfect-square alphabet (e.g. 25 or 36; got %d). "
               "Exclude a letter so it is 25 (e.g. -excludeletter J).\n\n", g_alpha);
        return;
    }

    // The second message arrives via -cipher2 (cfg->twincipher_str); required.
    if (!cfg->twincipher_str || cfg->twincipher_str[0] == '\0') {
        printf("\n\nERROR: Twin Bifid needs a SECOND ciphertext -- supply it with "
               "-cipher2 <file> (the two messages share one keyed square).\n\n");
        return;
    }

    int n1 = cipher_len;
    SymbolTable symtab2;
    int n2 = decode_cipher(cfg->twincipher_str, cfg, g_twin_bifid_c2, &symtab2);

    if (n1 < 4 || n2 < 4) {
        printf("\n\nERROR: both Twin Bifid ciphertexts must be at least 4 letters "
               "(got %d and %d).\n\n", n1, n2);
        return;
    }
    if (n1 + n2 > MAX_CIPHER_LENGTH) {
        printf("\n\nERROR: combined Twin Bifid ciphertext length %d exceeds the %d limit.\n\n",
               n1 + n2, MAX_CIPHER_LENGTH);
        return;
    }
    // Both ciphertexts must be solid letters (reject any sentinel: space/punctuation).
    for (int i = 0; i < n1; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext 1 has a non-alphabet symbol at position %d; "
                   "Twin Bifid ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }
    for (int i = 0; i < n2; i++)
        if (g_twin_bifid_c2[i] < 0 || g_twin_bifid_c2[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext 2 has a non-alphabet symbol at position %d; "
                   "Twin Bifid ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    TwinBifidScratch scratch;
    scratch.side = side;
    scratch.grid_size = g_alpha;
    scratch.c1 = cipher_indices; scratch.n1 = n1;
    scratch.c2 = g_twin_bifid_c2; scratch.n2 = n2;

    // Candidate periods: pinned (-period / -period2) or estimated per message, then crossed.
    int cand1[TWIN_BIFID_MAX_PER], cand2[TWIN_BIFID_MAX_PER];
    int nc1 = twin_bifid_periods_for(cipher_indices, n1, cfg->period_present, cfg->period, cfg, cand1);
    int nc2 = twin_bifid_periods_for(g_twin_bifid_c2, n2, cfg->period2_present, cfg->period2, cfg, cand2);
    int np = 0;
    for (int i = 0; i < nc1 && np < TWIN_BIFID_MAX_PAIRS; i++)
        for (int j = 0; j < nc2 && np < TWIN_BIFID_MAX_PAIRS; j++) {
            scratch.p1[np] = cand1[i];
            scratch.p2[np] = cand2[j];
            np++;
        }
    scratch.n_pairs = np;

    if (cfg->verbose)
        printf("\ntwin-bifid: msg1 %d letters, msg2 %d letters, %d-letter alphabet %s, "
               "%d candidate period pair(s)\n", n1, n2, g_alpha, g_idx_to_char_arr, np);

    // Build the combined ciphertext the engine sees as ctx->cipher (length n1+n2). The
    // decrypt hook splits it back via the scratch; this keeps ctx->cipher fully valid for
    // its stated length and gives print_cipher a sensible stream.
    for (int i = 0; i < n1; i++) g_twin_bifid_combined[i] = cipher_indices[i];
    for (int i = 0; i < n2; i++) g_twin_bifid_combined[n1 + i] = g_twin_bifid_c2[i];

    // Pass the concatenated length n1+n2 as the engine's scoring length (Tri-Square pattern);
    // the decrypt hook emits n1+n2 plaintext symbols the engine n-gram-scores.
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        g_twin_bifid_combined, n1 + n2, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // twin_bifid_report fills it (may be NULL for CLI use)

    run_solver(&TWIN_BIFID_MODEL, &ctx);
}
