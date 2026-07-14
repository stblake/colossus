#include "homophonic_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Homophonic substitution solver (TYPE homophonic)
// =====================================================================
//
// A monoalphabetic-in-meaning substitution whose CIPHERTEXT alphabet is larger than
// the plaintext alphabet: each plaintext letter is enciphered by any of several
// distinct ciphertext symbols (its homophones), chosen to flatten the ciphertext
// frequency profile (Zodiac-408 style). The ciphertext is decoded (decode_cipher)
// into a sequence of symbol ids 0..N-1 indexing a SymbolTable; the key is the
// many-to-one map symbol_id -> plaintext letter, so decrypted[i] = key[cipher[i]].
//
// There is no period or transposition to recover -- positions are preserved -- so the
// solver just hill-climbs the N-entry map against the n-gram score, exactly like the
// other CipherModels.
//
// ANTI-COLLAPSE. A homophonic map is free to fold many symbols onto E/T/A... to tile
// high-frequency n-grams -- a low-entropy fixed point that can out-score (or, at high
// symbol-count, tie) the true plaintext on raw n-grams. Two mutually-exclusive terms
// suppress it (pick by symbol-count / redundancy regime):
//   * -weightmono (default 1.0): ADDITIVE chi-squared monogram penalty. Best when the
//     collapse clearly out-scores truth on n-grams (short, few symbols): e.g. Zodiac Z408
//     (54 symbols) needs weightmono ~1.5.
//   * -weightentropy (default 0): AZDecrypt's MULTIPLICATIVE fitness (ngram - floor)*H^w,
//     H = Shannon entropy of the decrypt (see state_score). Best when truth barely beats
//     collapse on n-grams (long, MANY symbols) so an additive penalty just lands
//     frequency-matched near-misses: e.g. Beale cipher 2 (182 symbols, 762 chars) is
//     cracked blind by -weightentropy 1.5 -weightmono 0 (long slow anneal; the entropy
//     score scale is O(10-50), so inittemp must be O(1), not the ~0.02 of the mono scale). It plugs into the shared shotgun/anneal engine (run_solver):
// SHAPE_ANNEAL (Metropolis) acceptance, shotgun restarts, backtracking. Seeds are
// frequency-flattening (symbols drawn from the English monogram distribution, so
// common letters naturally receive more homophones); the move set reassigns one
// symbol's letter (dominant) or swaps two symbols' letters.

typedef struct {
    SymbolTable *tab;     // the interned ciphertext symbols (for display)
    int          n_symbols;

    // --- incremental-scoring caches (see the fast-path hooks below) ---
    // These describe the CURRENT (accepted) decryption, kept live by sync_caches /
    // commit_neighbor so each neighbour is scored as a delta over only the touched
    // cipher positions instead of an O(cipher_len) re-decrypt + n-gram rescan.
    double  scale;                 // n-gram scale factor (matches ngram_score())
    double  ngsum;                 // running raw sum of ngram_data[] over all windows
    int     counts[ALPHABET_SIZE]; // letter histogram of the current decryption
    // Position index: pos[pos_off[s] .. pos_off[s+1]) lists the cipher positions
    // carrying symbol s (so a one-symbol move maps straight to the positions, and
    // thus the n-gram windows, it changes). Built once per solve.
    int    *pos;                   // length cipher_len
    int    *pos_off;               // length n_symbols + 1
    // Per-call scratch for the touched-window set (win_mark stays all-zero between
    // calls; score_neighbor clears every entry it set before returning).
    char   *win_mark;              // length cipher_len (>= n_windows)
    int    *win_list;              // length cipher_len
    // Pending delta stashed by score_neighbor, applied by commit_neighbor.
    double  pend_ngsum;
    int     pend_counts[ALPHABET_SIZE];
    int     pend_nsym;             // number of changed symbols
    int    *pend_sym;              // length n_symbols: which symbols changed
    int    *pend_newc;             // length n_symbols: their new plaintext letters
} HomophonicScratch;

// One config: the whole map is climbed at once. period carries the key length.
static int homophonic_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const HomophonicScratch *h = (const HomophonicScratch *) ctx->model_scratch;
    if (cap < 1) return 0;
    out[0].period = h->n_symbols;
    out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Frequency-flattening seed: draw each symbol's plaintext letter from the English
// monogram distribution. Randomised (so shotgun restarts diversify) yet biased so the
// recovered ciphertext frequencies start out roughly English-shaped.
static void homophonic_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    int n = cc->period;
    double cum[ALPHABET_SIZE], total = 0.;
    for (int c = 0; c < g_alpha; c++) { total += g_monograms[c]; cum[c] = total; }
    for (int s = 0; s < n; s++) {
        double r = frand() * total;
        int c = 0;
        while (c < g_alpha - 1 && r > cum[c]) c++;
        st->key[s] = c;
    }
    st->key_len = n;
}

// The anti-collapse penalty: chi-squared of the decrypted letter-frequency profile
// against English monograms. Unlike a 26->26 substitution (a bijection, which cannot
// pile multiple symbols onto one letter), a homophonic map is free to fold many symbols
// onto E/T/A... to tile high-frequency n-grams -- a fixed point that out-scores the
// true plaintext on raw n-grams alone. Penalising the resulting (wildly non-English)
// monogram distribution removes that fixed point. Returned as a positive quantity to be
// SUBTRACTED from the score.
static double homophonic_penalty(const SolverCtx *ctx, const int *dec) {
    if (ctx->cfg->weight_monogram <= 1.e-9) return 0.0;
    return ctx->cfg->weight_monogram * chi_squared((int *) dec, ctx->cipher_len);
}

// Neighbour move (AZDecrypt-style). The dominant move is a single random symbol
// reassignment -- pick one symbol, give it a random plaintext letter -- accepted or
// rejected by the engine's Metropolis annealing. The old greedy coordinate step
// (try all 26 letters for one symbol, every trial a full re-decrypt + n-gram rescan)
// was the dominant cost AND actively drove the homophonic collapse it then needed the
// chi-squared penalty to undo; a single random reassignment scored incrementally is
// hundreds of times cheaper per iteration and lets annealing find the basin. Two
// auxiliary moves aid exploration: a symbol-pair letter swap, and a letter-class swap
// (exchange the WHOLE homophone classes of two plaintext letters at once, to cross
// equal-frequency ambiguities like W<->M that single-symbol moves cannot, since
// flipping one symbol first makes it worse). Every move stays a small edit to the
// key, so score_neighbor can score it as a delta over the touched positions only.
static void homophonic_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                               SolverState *st, bool *force_primary) {
    (void) ctx; (void) force_primary;
    int n = cc->period;
    if (n < 1) return;
    double r = frand();
    if (r < 0.85) {
        // Single random reassignment (the workhorse move).
        st->key[rand_int(0, n)] = rand_int(0, g_alpha);
    } else if (n >= 2 && r < 0.93) {
        // Swap two symbols' letters (fine-grained exploration).
        int a = rand_int(0, n), b = rand_int(0, n);
        int t = st->key[a]; st->key[a] = st->key[b]; st->key[b] = t;
    } else {
        // Letter-class swap (cross equal-frequency two-letter ambiguities).
        int a = rand_int(0, g_alpha), b = rand_int(0, g_alpha);
        if (a != b)
            for (int s = 0; s < n; s++) {
                if (st->key[s] == a) st->key[s] = b;
                else if (st->key[s] == b) st->key[s] = a;
            }
    }
}

// Packed big-endian base-g_alpha index of the n-gram window at start `w`, reading
// each position's letter through `letter(p)`. Matches ngram_score()'s packing exactly
// (idx = sum_j letter(w+j) * g_alpha^(ng-1-j)), so the same ngram_data[] entry is
// selected -- the incremental ngsum stays bit-identical to the full rescore.
#define HOMO_WINDOW_INDEX(w, ng, letter) ({                 \
    int _idx = 0;                                           \
    for (int _j = 0; _j < (ng); _j++)                       \
        _idx = _idx * g_alpha + letter((w) + _j);           \
    _idx; })

// Rebuild the incremental caches (letter histogram + running n-gram sum) from a full
// decryption `dec`. Called once per restart/backtrack, after the engine has reset the
// current state with a full decrypt -- the per-iteration deltas ride on top of this.
static void homophonic_sync_caches(const SolverCtx *ctx, const SolverConfig *cc, const int *dec) {
    (void) cc;
    HomophonicScratch *h = (HomophonicScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, ng = ctx->cfg->ngram_size;
    tally((int *) dec, len, h->counts, ALPHABET_SIZE);
    h->ngsum = 0.0;
    int n_windows = len - ng + 1;
    #define dec_at(q) (dec[q])
    for (int w = 0; w < n_windows; w++)
        h->ngsum += ctx->ngram_data[HOMO_WINDOW_INDEX(w, ng, dec_at)];
    #undef dec_at
}

// chi-squared of a letter histogram vs English monograms -- the same quantity as
// chi_squared() (utils.c) but from a precomputed count[] (no re-tally).
static double homophonic_chi2_from_counts(const int *counts, int len) {
    double chi2 = 0.0;
    for (int c = 0; c < g_alpha; c++) {
        double f = ((double) counts[c]) / len;
        chi2 += pow(f - g_monograms[c], 2) / g_monograms[c];
    }
    return chi2;
}

// Score the neighbour `loc` (= `cur` after one perturb) incrementally. Diffs loc vs
// cur to find the changed symbols, recomputes only the n-gram windows touching those
// symbols' positions and the two affected histogram bins, and reassembles the score
// to equal decrypt()+state_score()+score_adjust for loc. Stashes the delta (new
// ngsum + new counts + changed symbols) for commit_neighbor.
static double homophonic_score_neighbor(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *cur, const SolverState *loc, const int *cur_dec, double cur_score) {
    (void) cur_score;
    HomophonicScratch *h = (HomophonicScratch *) ctx->model_scratch;
    ColossusConfig *cfg = ctx->cfg;
    int len = ctx->cipher_len, ng = cfg->ngram_size, n = cc->period;
    const int *cipher = ctx->cipher;

    // 1. Changed symbols, and the histogram delta (each symbol's whole class of
    //    positions moves from its old letter to its new one).
    for (int c = 0; c < g_alpha; c++) h->pend_counts[c] = h->counts[c];
    h->pend_nsym = 0;
    for (int s = 0; s < n; s++) {
        if (loc->key[s] == cur->key[s]) continue;
        int oldc = cur->key[s], newc = loc->key[s];
        int m = h->pos_off[s + 1] - h->pos_off[s];
        h->pend_counts[oldc] -= m;
        h->pend_counts[newc] += m;
        h->pend_sym[h->pend_nsym]  = s;
        h->pend_newc[h->pend_nsym] = newc;
        h->pend_nsym++;
    }

    // 2. n-gram delta over the windows touching any changed position. Mark the
    //    distinct touched window starts, then for each re-score old vs new.
    double dsum = 0.0;
    int n_windows = len - ng + 1;
    if (cfg->weight_ngram > 1.e-4 && n_windows > 0) {
        int nlist = 0;
        for (int k = 0; k < h->pend_nsym; k++) {
            int s = h->pend_sym[k];
            for (int pi = h->pos_off[s]; pi < h->pos_off[s + 1]; pi++) {
                int p = h->pos[pi];
                int lo = p - ng + 1; if (lo < 0) lo = 0;
                int hi = p;          if (hi > n_windows - 1) hi = n_windows - 1;
                for (int w = lo; w <= hi; w++)
                    if (!h->win_mark[w]) { h->win_mark[w] = 1; h->win_list[nlist++] = w; }
            }
        }
        // loc's letter at position q (only changed symbols differ from cur_dec).
        #define dec_at(q) (cur_dec[q])
        #define loc_at(q) (loc->key[cipher[q]])
        for (int i = 0; i < nlist; i++) {
            int w = h->win_list[i];
            int old_idx = HOMO_WINDOW_INDEX(w, ng, dec_at);
            int new_idx = HOMO_WINDOW_INDEX(w, ng, loc_at);
            dsum += ctx->ngram_data[new_idx] - ctx->ngram_data[old_idx];
            h->win_mark[w] = 0;            // clear for next call
        }
        #undef loc_at
        #undef dec_at
    }
    h->pend_ngsum = h->ngsum + dsum;

    // 3. Reassemble the score exactly as engine_score (state_score + score_adjust).
    double ng_score = 0.0;
    if (cfg->weight_ngram > 1.e-4 && n_windows > 0)
        ng_score = h->scale * h->pend_ngsum / (len - ng);

    // Mirror state_score's multiplicative entropy term (AZDecrypt fitness) on the fast
    // path: fitness_ngram = (ng_score - floor) * H^weight_entropy, with H the Shannon
    // entropy of the neighbour's letter distribution. Computed from the maintained
    // histogram pend_counts (no re-tally) using the SAME formula as entropy() (utils.c:
    // natural log, over ALPHABET_SIZE bins, divided by len) so the incremental score
    // stays bit-identical to the full state_score. Default weight_entropy == 0 skips it.
    if (cfg->weight_ngram > 1.e-4 && n_windows > 0 && cfg->weight_entropy > 1.e-4) {
        double H = 0.0;
        for (int c = 0; c < ALPHABET_SIZE; c++) {
            if (h->pend_counts[c] > 0) {
                double f = ((double) h->pend_counts[c]) / len;
                H -= f * log(f);
            }
        }
        ng_score = (ng_score - g_ngram_floor) * pow((double) H, (double) cfg->weight_entropy);
    }

    double score;
    if (ctx->n_cribs > 0) {
        double crib = 0.0;
        if (cfg->weight_crib > 1.e-4) {
            for (int i = 0; i < ctx->n_cribs; i++) {
                int got = loc->key[cipher[ctx->crib_positions[i]]];
                int diff = abs(got - ctx->crib_indices[i]);
                crib += (diff == 0) ? 1.0 : 1.0 / (1.0 + diff * diff);
            }
            crib /= (double) ctx->n_cribs;
        }
        score = (cfg->weight_ngram * ng_score + cfg->weight_crib * crib)
                / (cfg->weight_ngram + cfg->weight_crib);
    } else {
        score = ng_score;
    }

    if (cfg->weight_monogram > 1.e-9)
        score -= cfg->weight_monogram * homophonic_chi2_from_counts(h->pend_counts, len);

    return score;
}

// Apply the delta stashed by the most recent score_neighbor: advance cur_dec at the
// changed symbols' positions and adopt the precomputed histogram + n-gram sum.
static void homophonic_commit_neighbor(const SolverCtx *ctx, const SolverConfig *cc, int *cur_dec) {
    (void) cc;
    HomophonicScratch *h = (HomophonicScratch *) ctx->model_scratch;
    for (int k = 0; k < h->pend_nsym; k++) {
        int s = h->pend_sym[k], c = h->pend_newc[k];
        for (int pi = h->pos_off[s]; pi < h->pos_off[s + 1]; pi++)
            cur_dec[h->pos[pi]] = c;
    }
    for (int c = 0; c < g_alpha; c++) h->counts[c] = h->pend_counts[c];
    h->ngsum = h->pend_ngsum;
}

static void homophonic_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}

// Per-thread scratch clone (enables -nthreads on the incremental fast path). Shares
// the read-only fields built once in solve_homophonic (tab, scale, and the position
// index pos/pos_off -- never written during the search) but allocates this worker its
// OWN copies of every buffer the search mutates: the live caches (counts/ngsum are
// inline in the struct, zeroed here and rebuilt by sync_caches each restart) and the
// per-call scratch (win_mark/win_list/pend_sym/pend_newc). Returns NULL on OOM.
static void *homophonic_scratch_clone(const SolverCtx *ctx) {
    const HomophonicScratch *src = (const HomophonicScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, n = src->n_symbols;
    HomophonicScratch *h = malloc(sizeof *h);
    if (!h) return NULL;
    memset(h, 0, sizeof *h);
    // Shared, read-only during the search.
    h->tab = src->tab;
    h->n_symbols = src->n_symbols;
    h->scale = src->scale;
    h->pos = src->pos;          // position index: built once, only read -> shared
    h->pos_off = src->pos_off;
    // Private, mutated during the search.
    h->win_mark  = calloc(len > 0 ? len : 1, sizeof(char));
    h->win_list  = malloc(sizeof(int) * (len > 0 ? len : 1));
    h->pend_sym  = malloc(sizeof(int) * (n > 0 ? n : 1));
    h->pend_newc = malloc(sizeof(int) * (n > 0 ? n : 1));
    if (!h->win_mark || !h->win_list || !h->pend_sym || !h->pend_newc) {
        free(h->win_mark); free(h->win_list); free(h->pend_sym); free(h->pend_newc);
        free(h);
        return NULL;
    }
    return h;
}

// Free a clone: only its private allocations + the struct. pos/pos_off/tab are shared
// with the owner scratch (freed by solve_homophonic), so they are NOT freed here.
static void homophonic_scratch_free(void *scratch) {
    HomophonicScratch *h = (HomophonicScratch *) scratch;
    if (!h) return;
    free(h->win_mark); free(h->win_list); free(h->pend_sym); free(h->pend_newc);
    free(h);
}

static void homophonic_decrypt(const SolverCtx *ctx, const SolverConfig *cc,
                               SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    for (int i = 0; i < ctx->cipher_len; i++) out[i] = st->key[ctx->cipher[i]];
    *score_adjust = -homophonic_penalty(ctx, out);   // engine adds this to state_score
}

static void homophonic_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) st;
    const HomophonicScratch *h = (const HomophonicScratch *) ctx->model_scratch;
    char params[64];
    snprintf(params, sizeof(params), "symbols=%d", h->n_symbols);
    report_transposition_verbose(ctx, score, decrypted, stats, params);
}

static void homophonic_report(const SolverCtx *ctx, const SolverConfig *cc,
                              const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const HomophonicScratch *h = (const HomophonicScratch *) ctx->model_scratch;
    SymbolTable *tab = h->tab;
    int n = h->n_symbols, len = ctx->cipher_len;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    printf("\nResult Score: %.2f | Words: %d | symbols=%d\n", score, n_words_found, n);

    print_cipher(ctx->cipher, len, tab);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    // Recovered homophone classes: each plaintext letter and the symbols decoding to it.
    printf("\nhomophone key (plaintext <- symbols):\n");
    for (int c = 0; c < g_alpha; c++) {
        int any = 0;
        for (int s = 0; s < n; s++) if (st->key[s] == c) {
            if (!any) { printf("  %c <-", index_to_char(c)); any = 1; }
            printf(" %s", tab->tokens[s]);
        }
        if (any) printf("\n");
    }

    // One-liner summary: >>> score, [words,] type, symbols=N, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, symbols=%d, ", score, n_words_found, cfg->cipher_type, n);
    else
        printf(">>> %.2f, %d, symbols=%d, ", score, cfg->cipher_type, n);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, tab);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel HOMOPHONIC_MODEL = {
    .name = "homophonic", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = homophonic_enumerate, .key_len = NULL,
    .seed = homophonic_seed, .perturb = homophonic_perturb, .copy_state = homophonic_copy,
    .decrypt = homophonic_decrypt, .report = homophonic_report,
    .report_verbose = homophonic_report_verbose,
    .sync_caches = homophonic_sync_caches,
    .score_neighbor = homophonic_score_neighbor,
    .commit_neighbor = homophonic_commit_neighbor,
    .scratch_clone = homophonic_scratch_clone,
    .scratch_free = homophonic_scratch_free,
};

void solve_homophonic(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SymbolTable *tab) {

    (void) ciphertext_str;
    if (cipher_len < 4 || tab == NULL || tab->n < 1) {
        printf("\n\nERROR: ciphertext too short for a homophonic solve.\n\n");
        return ;
    }
    if (cfg->verbose)
        printf("\nhomophonic: %d positions, %d distinct symbols\n", cipher_len, tab->n);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);

    int n_sym = tab->n;
    HomophonicScratch scratch;
    memset(&scratch, 0, sizeof scratch);
    scratch.tab = tab;
    scratch.n_symbols = n_sym;
    // n-gram scale factor: must match ngram_score() (legacy g_alpha^ngramsize, or 1
    // in -logprob mode where the table already holds O(1) log10 probabilities).
    scratch.scale = g_ngram_logprob ? 1.0 : pow((double) g_alpha, cfg->ngram_size);

    // Position index: bucket each cipher position by its symbol so a one-symbol move
    // maps straight to the positions (and thus n-gram windows) it changes. Built once.
    scratch.pos     = malloc(sizeof(int) * (cipher_len > 0 ? cipher_len : 1));
    scratch.pos_off = malloc(sizeof(int) * (n_sym + 1));
    scratch.win_mark = calloc(cipher_len > 0 ? cipher_len : 1, sizeof(char));
    scratch.win_list = malloc(sizeof(int) * (cipher_len > 0 ? cipher_len : 1));
    scratch.pend_sym  = malloc(sizeof(int) * (n_sym > 0 ? n_sym : 1));
    scratch.pend_newc = malloc(sizeof(int) * (n_sym > 0 ? n_sym : 1));
    if (!scratch.pos || !scratch.pos_off || !scratch.win_mark ||
        !scratch.win_list || !scratch.pend_sym || !scratch.pend_newc) {
        printf("\n\nERROR: out of memory in homophonic solve.\n\n");
        free(scratch.pos); free(scratch.pos_off); free(scratch.win_mark);
        free(scratch.win_list); free(scratch.pend_sym); free(scratch.pend_newc);
        return ;
    }
    // Counting sort of positions by symbol id into pos[]/pos_off[].
    for (int s = 0; s <= n_sym; s++) scratch.pos_off[s] = 0;
    for (int i = 0; i < cipher_len; i++) scratch.pos_off[cipher_indices[i] + 1]++;
    for (int s = 0; s < n_sym; s++) scratch.pos_off[s + 1] += scratch.pos_off[s];
    {
        int *cursor = malloc(sizeof(int) * (n_sym > 0 ? n_sym : 1));
        for (int s = 0; s < n_sym; s++) cursor[s] = scratch.pos_off[s];
        for (int i = 0; i < cipher_len; i++) {
            int s = cipher_indices[i];
            scratch.pos[cursor[s]++] = i;
        }
        free(cursor);
    }

    ctx.model_scratch = &scratch;
    run_solver(&HOMOPHONIC_MODEL, &ctx);

    free(scratch.pos); free(scratch.pos_off); free(scratch.win_mark);
    free(scratch.win_list); free(scratch.pend_sym); free(scratch.pend_newc);
}


