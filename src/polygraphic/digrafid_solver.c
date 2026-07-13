#include "digrafid_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Digrafid solver (TYPE digrafid)
// =====================================================================
//
// Digrafid is a digraphic fractionation cipher over TWO independently keyed 27-symbol
// alphabets (A..Z + '#'): a horizontal grid H (3x9) and a vertical grid V (9x3). Breaking
// it is two coupled problems: recover the period and recover both grids.
//
// KEY-STRUCTURED SEARCH (not the full 27!^2 free permutation). ACA Digrafid grids are KEYED
// ALPHABETS -- a keyword (duplicates dropped) followed by the rest of the alphabet ascending
// -- exactly as the generator builds them (digrafid_grid_from_keyword). So instead of
// annealing two free 54-cell permutations (the Two-Square square break, which needs ~700-800
// letters of n-gram signal to pin a grid statistically), we search the SAME space the cipher
// is keyed in: the state is two keyed-alphabet SEQUENCES (H = key[0..26], V = key[27..53]),
// each maintained as "keyword prefix of length kw + ascending tail" by the SAME core moves
// the Quagmire solver uses -- random_keyword() seeds them and perturbate_keyword() perturbs
// them (a 20% in-keyword swap / 80% keyword<->tail swap that re-inserts the displaced letter
// into the tail in sorted order, so the keyed-alphabet invariant always holds). Each move
// perturbs ONE sequence (chosen uniformly); the grids are rebuilt from the sequences at
// decrypt (H row-major, V column-major). The keyword lengths kwH, kwV live in st->aux[0..1]
// and are sampled per restart from [DIGRAFID_KW_MIN .. DIGRAFID_KW_MAX] -- RESTARTS are the
// lever that covers the keyword-length range (the n-gram score across restarts/periods picks
// the winner), the same way Bazeries sweeps its digit count and the period is swept here.
// Collapsing each grid from a free permutation to a short keyword shrinks the keyspace by
// orders of magnitude and makes the structured tail a strong prior, so recovery no longer
// needs ~800 letters -- it tracks the keyword, the way an ACA solver does. No anti-collapse
// penalty is needed (every keyed alphabet is a bijection, so the whole digraph map is), so
// score_adjust stays 0 and it rides the generic state_score; it still effectively needs
// -logprob. (A keyword longer than DIGRAFID_KW_MAX, or a non-keyed/fully-shuffled grid, is
// not representable -- but every ACA Digrafid, and the generator, uses a keyed alphabet.)
//
// The period is recovered by an index-of-coincidence test (digrafid_estimate_periods):
// the ciphertext is laid out in 2P lanes -- one per (digraph-position-in-group, first/
// second-letter role) -- and at the true period each lane is coordinate-constrained, so
// the mean per-lane IoC stands out (it also peaks at multiples, which decrypt to gibberish
// and lose on the n-gram score). We anneal the top-K periods (one engine config each) and
// the n-gram score picks the winner. -period pins a single period.

#define DIGRAFID_MAX_PERIODS 64

// Keyword-length search range for each keyed alphabet. A keyed alphabet of true keyword
// length L is representable iff kw >= L (a too-short kw cannot place the keyword tail in
// sorted order), and a tighter kw (closer to L) is a stronger prior / smaller space. The
// range covers realistic ACA keyword lengths (the test corpus tops out at BERLINCLOCK = 11)
// while keeping a tail of >= DIGRAFID_GRID - DIGRAFID_KW_MAX = 14 cells in fixed ascending
// order -- a genuine structural constraint, never the free permutation (kw == DIGRAFID_GRID).
#define DIGRAFID_KW_MIN 3
#define DIGRAFID_KW_MAX 13

typedef struct {
    int grid_size;                       // cells per grid (== g_alpha == 27)
    int n_periods;                       // number of candidate periods
    int periods[DIGRAFID_MAX_PERIODS];   // the candidate periods (config order)
} DigrafidScratch;

// --- period estimation --------------------------------------------------------
//
// For each trial period P, average the Index of Coincidence over the 2P lanes formed by
// (digraph index mod P, first/second-letter role). At the true period every lane draws a
// coordinate-constrained sub-distribution, so the mean per-lane IoC peaks. Returns the top
// n_want periods (descending) into out[].
int digrafid_estimate_periods(int cipher[], int len, int min_p, int max_p,
                              int n_want, int out[], bool verbose) {
    static double ioc[MAX_CIPHER_LENGTH];   // ioc[p] for p in [min_p..max_p]
    int ndig = len / 2;
    if (max_p > ndig / 2) max_p = ndig / 2;
    if (max_p < min_p) max_p = min_p;
    if (n_want < 1) n_want = 1;
    if (n_want > DIGRAFID_MAX_PERIODS) n_want = DIGRAFID_MAX_PERIODS;

    for (int p = min_p; p <= max_p; p++) {
        double total = 0.0;
        int nlanes = 0;
        for (int lane = 0; lane < 2 * p; lane++) {
            int counts[MAX_ALPHABET_SIZE];
            for (int s = 0; s < g_alpha; s++) counts[s] = 0;
            int N = 0;
            int dpos = lane / 2, role = lane & 1;        // digraph-position-in-group, role
            for (int d = dpos; d < ndig; d += p) {
                int sym = cipher[2 * d + role];
                if (sym >= 0 && sym < g_alpha) { counts[sym]++; N++; }
            }
            if (N > 1) {
                long acc = 0;
                for (int s = 0; s < g_alpha; s++) acc += (long) counts[s] * (counts[s] - 1);
                total += (double) acc / ((double) N * (N - 1));
                nlanes++;
            }
        }
        ioc[p] = (nlanes > 0) ? total / nlanes : 0.0;
    }

    if (verbose) {
        printf("\nDigrafid period estimate (mean per-lane IoC):\n  period\tIoC\n");
        for (int p = min_p; p <= max_p; p++) printf("  %d\t%.4f\n", p, ioc[p]);
    }

    // Selection by repeated max (n_want is small); mark taken periods with -1.
    int n = 0;
    for (; n < n_want; n++) {
        int best_p = -1;
        double best_ioc = -1.0;
        for (int p = min_p; p <= max_p; p++)
            if (ioc[p] >= 0.0 && ioc[p] > best_ioc) { best_ioc = ioc[p]; best_p = p; }
        if (best_p < 0) break;
        out[n] = best_p;
        ioc[best_p] = -1.0;
    }
    if (verbose) {
        printf("  -> annealing periods:");
        for (int i = 0; i < n; i++) printf(" %d", out[i]);
        printf("\n");
    }
    return n;
}

// One config per candidate period; period carries the block size, key the two grids.
static int digrafid_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const DigrafidScratch *d = (const DigrafidScratch *) ctx->model_scratch;
    int n = d->n_periods;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = d->periods[i];
        out[i].j = 0; out[i].k = 0; out[i].aux[0] = 0; out[i].aux[1] = 0;
    }
    return n;
}

// Build the two grids from the two keyed-alphabet SEQUENCES carried in the state. seqH/seqV
// are full DIGRAFID_GRID-length keyed sequences (keyword prefix + ascending tail); H is
// filled row-major, V column-major. Passing the whole sequence as the "keyword" makes
// digrafid_grid_from_keyword a pure sequence->grid placement (every cell present, no tail
// appended), so this is exactly the grid the generator built from the same keyed alphabet.
static void digrafid_build_grids(const int *seqH, const int *seqV, int *gridH, int *gridV) {
    digrafid_grid_from_keyword(seqH, DIGRAFID_GRID, gridH, DIGRAFID_HROWS, DIGRAFID_HCOLS, 0);
    digrafid_grid_from_keyword(seqV, DIGRAFID_GRID, gridV, DIGRAFID_VROWS, DIGRAFID_VCOLS, 1);
}

// Re-establish the keyed-alphabet invariant after the keyword/tail boundary kw moves: keep
// the prefix seq[0..kw-1] as the keyword (any order), then refill seq[kw..n-1] with the
// remaining symbols in ascending order. (n = 27, so the O(n^2) refill is negligible.)
void digrafid_canonicalize(int *seq, int n, int kw) {
    char used[DIGRAFID_GRID];
    for (int i = 0; i < n; i++) used[i] = 0;
    for (int i = 0; i < kw; i++) used[seq[i]] = 1;
    int m = kw;
    for (int s = 0; s < n && m < n; s++) if (!used[s]) seq[m++] = s;
}

// One structure-preserving move on a single keyed alphabet `seq` of n cells whose keyword
// prefix has length *kw. The keyed-alphabet invariant (sorted tail past *kw) always holds on
// exit. Three move classes, tuned for the Digrafid fractionation landscape (coarser than a
// monoalphabetic substitution, so the smooth in-keyword reorder carries more of the fine
// convergence than Quagmire's 20/80 split):
//   ~4%  : change the keyword length by +-1 (re-canonicalising the tail) -- lets ONE restart
//          find the true keyword length instead of relying on the seed to hit it;
//   ~48% : swap a keyword letter with a tail letter and re-sort the tail (COARSE set search --
//          brings new letters into / out of the keyword; reshuffles part of the tail);
//   ~48% : swap two keyword letters (SMOOTH: reorders the keyword, 2 cells, no tail shift) --
//          the fine move that locks the exact order once the keyword set is right.
void digrafid_move_seq(int *seq, int n, int *kw) {
    double r = frand();
    if (r < 0.04) {                                   // grow / shrink the keyword
        int k = *kw;
        if (k <= DIGRAFID_KW_MIN) k++;
        else if (k >= DIGRAFID_KW_MAX) k--;
        else k += (frand() < 0.5) ? 1 : -1;
        *kw = k;
        digrafid_canonicalize(seq, n, k);
    } else if (r < 0.52) {                            // keyword <-> tail (coarse set search)
        int i = rand_int(0, *kw), j = rand_int(*kw, n);
        int t = seq[i]; seq[i] = seq[j]; seq[j] = t;
        digrafid_canonicalize(seq, n, *kw);           // keep the tail sorted
    } else {                                          // in-keyword reorder (smooth)
        int i = rand_int(0, *kw), j = rand_int(0, *kw);
        int t = seq[i]; seq[i] = seq[j]; seq[j] = t;
    }
}

// Seed: two independent random keyed alphabets. random_keyword lays a random keyword of
// length kw (distinct symbols) then the rest of the alphabet ascending; kw is sampled per
// restart from [DIGRAFID_KW_MIN .. DIGRAFID_KW_MAX] and stashed in st->aux. The keyword length
// is then ALSO perturbed within the restart (digrafid_move_seq), so the seed only sets a
// starting length -- restarts + the in-climb length move together cover the range.
static void digrafid_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    const DigrafidScratch *d = (const DigrafidScratch *) ctx->model_scratch;
    int g = d->grid_size;
    for (int s = 0; s < 2; s++) {
        int kw = rand_int(DIGRAFID_KW_MIN, DIGRAFID_KW_MAX + 1);   // [MIN .. MAX]
        random_keyword(st->key + s * g, g, kw);
        st->aux[s] = kw;
    }
    st->key_len = 2 * g;
}

// Neighbour move: perturb ONE of the two keyed alphabets, chosen uniformly, keeping its
// keyword length in st->aux up to date.
static void digrafid_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                             SolverState *st, bool *force_primary) {
    const DigrafidScratch *d = (const DigrafidScratch *) ctx->model_scratch;
    (void) cc; (void) force_primary;
    int g = d->grid_size;
    int s = rand_int(0, 2);
    digrafid_move_seq(st->key + s * g, g, &st->aux[s]);
}

static void digrafid_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    // Both keyed-alphabet sequences (DIGRAFID_STATE cells) plus the two keyword lengths.
    for (int i = 0; i < DIGRAFID_STATE; i++) dst->key[i] = src->key[i];
    dst->aux[0] = src->aux[0];   // kwH
    dst->aux[1] = src->aux[1];   // kwV
}

static void digrafid_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                  SolverState *st, int *out, double *score_adjust) {
    const DigrafidScratch *d = (const DigrafidScratch *) ctx->model_scratch;
    int gridH[DIGRAFID_GRID], gridV[DIGRAFID_GRID];
    digrafid_build_grids(st->key, st->key + d->grid_size, gridH, gridV);
    digrafid_decrypt(ctx->cipher, ctx->cipher_len, gridH, gridV, cc->period, out);
    *score_adjust = 0.0;
}

// Render one rows x cols grid as an indented box of letters.
static void digrafid_print_grid(const int grid[], int rows, int cols) {
    for (int r = 0; r < rows; r++) {
        printf("    ");
        for (int c = 0; c < cols; c++) printf("%c ", index_to_char(grid[r * cols + c]));
        printf("\n");
    }
}

static void digrafid_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const DigrafidScratch *d = (const DigrafidScratch *) ctx->model_scratch;
    int gridH[DIGRAFID_GRID], gridV[DIGRAFID_GRID];
    digrafid_build_grids(st->key, st->key + d->grid_size, gridH, gridV);
    printf("\n  period %d, horizontal grid (3x9):\n", cc->period);
    digrafid_print_grid(gridH, DIGRAFID_HROWS, DIGRAFID_HCOLS);
    printf("  vertical grid (9x3):\n");
    digrafid_print_grid(gridV, DIGRAFID_VROWS, DIGRAFID_VCOLS);
    report_transposition_verbose(ctx, score, decrypted, stats, "digrafid");
}

static void digrafid_report(const SolverCtx *ctx, const SolverConfig *cc,
                            const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const DigrafidScratch *d = (const DigrafidScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, n = d->grid_size;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The two recovered grids, read row-major (H first, then V), rebuilt from the keyed-
    // alphabet sequences carried in the state.
    int gridH[DIGRAFID_GRID], gridV[DIGRAFID_GRID];
    digrafid_build_grids(st->key, st->key + n, gridH, gridV);
    char hstr[DIGRAFID_GRID + 1], vstr[DIGRAFID_GRID + 1];
    for (int i = 0; i < n; i++) { hstr[i] = index_to_char(gridH[i]); vstr[i] = index_to_char(gridV[i]); }
    hstr[n] = '\0'; vstr[n] = '\0';

    printf("\nResult Score: %.2f | Words: %d | period=%d | H=%s | V=%s\n",
        score, n_words_found, cc->period, hstr, vstr);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered horizontal grid (3x9, row major):\n");
    digrafid_print_grid(st->key, DIGRAFID_HROWS, DIGRAFID_HCOLS);
    printf("recovered vertical grid (9x3, row major):\n");
    digrafid_print_grid(st->key + n, DIGRAFID_VROWS, DIGRAFID_VCOLS);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = cc->period;   // report the recovered period here
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, period=, H=, V=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, period=%d, H=%s, V=%s, ",
            score, n_words_found, cfg->cipher_type, cc->period, hstr, vstr);
    else
        printf(">>> %.2f, %d, period=%d, H=%s, V=%s, ",
            score, cfg->cipher_type, cc->period, hstr, vstr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel DIGRAFID_MODEL = {
    .name = "digrafid", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = digrafid_enumerate, .key_len = NULL,
    .seed = digrafid_seed, .perturb = digrafid_perturb, .copy_state = digrafid_copy,
    .decrypt = digrafid_decrypt_hook, .report = digrafid_report,
    .report_verbose = digrafid_report_verbose,
};

void solve_digrafid(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Digrafid needs the 27-symbol alphabet (A..Z + '#'); the binary's main forces it, but
    // guard in case solve_digrafid is driven directly.
    if (g_alpha != DIGRAFID_GRID) {
        printf("\n\nERROR: Digrafid needs the 27-symbol alphabet (A..Z + '%c'; got %d).\n\n",
               DIGRAFID_EXTRA_CHAR, g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Digrafid solve.\n\n");
        return;
    }
    if (cipher_len % 2 != 0)
        printf("\nWARNING: odd ciphertext length (%d); a Digrafid ciphertext is always "
               "even. The trailing letter is left undecrypted.\n", cipher_len);
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Digrafid ciphertext must be solid symbols (A..Z + '%c').\n\n",
                   i, DIGRAFID_EXTRA_CHAR);
            return;
        }

    DigrafidScratch scratch;
    scratch.grid_size = g_alpha;

    // Candidate periods: pinned, or the estimator's top-K over [2 .. max_period].
    if (cfg->period_present) {
        scratch.periods[0] = cfg->period;
        scratch.n_periods = 1;
        if (cfg->verbose) printf("\ndigrafid: period pinned to %d\n", cfg->period);
    } else {
        int max_p = (cfg->max_period > 0) ? cfg->max_period : 20;
        int ndig = cipher_len / 2;
        if (max_p > ndig) max_p = ndig;
        int n_want = cfg->n_periods;
        scratch.n_periods = digrafid_estimate_periods(cipher_indices, cipher_len, 2, max_p,
            n_want, scratch.periods, cfg->verbose);
        if (scratch.n_periods < 1) { scratch.periods[0] = 2; scratch.n_periods = 1; }
    }

    if (cfg->verbose)
        printf("\ndigrafid: %d positions (%d digraphs), %d-symbol alphabet %s, %d candidate period(s)\n",
            cipher_len, cipher_len / 2, g_alpha, g_idx_to_char_arr, scratch.n_periods);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // digrafid_report fills it (may be NULL for CLI use)

    run_solver(&DIGRAFID_MODEL, &ctx);
}
