#include "aristocrat_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Aristocrat / Patristocrat solver (TYPES aristocrat, patristocrat)
// =====================================================================
//
// A simple monoalphabetic substitution: one 26-letter bijection key[] with pt = key[ct]. The
// Aristocrat preserves word divisions (spaces), the Patristocrat writes 5-letter groups; the
// enciphering is identical, so BOTH run this one solver, differing only in the report layout.
//
// SEARCH REPRESENTATION. The engine decodes the ciphertext into cipher_indices[] with letters
// 0..25 and non-letters as negative sentinels (space = -33). This solver derives a LETTERS-ONLY
// search stream (the >=0 entries) and keeps the full original layout for the Aristocrat's spaced
// report. The state is a 26-permutation key[ct]=pt; scoring is straight-through n-gram over the
// letters-only stream (spaces are transparent to the solver -- the "light" design: word divisions
// serve the report / dictionary / keyword recovery, not the solve gradient). Because every position
// is a letter, no sentinel bookkeeping is needed in the hot loop.
//
// ENGINE = the homophonic model made BIJECTIVE. It reuses the homophonic incremental fast path
// (sync_caches / score_neighbor / commit_neighbor + per-thread scratch_clone) so each neighbour is
// scored as a delta over only the touched n-gram windows -- fast and multithreaded. Two changes vs
// homophonic: (1) the move is a SWAP of two ciphertext letters' plaintext images (never a single
// reassignment, which would break the bijection), so a move changes exactly two symbols; (2) there
// is NO chi-squared anti-collapse penalty (a bijection cannot pile letters onto one image), so the
// decrypt hook leaves score_adjust = 0 and no letter histogram is tracked. Cribs are not used (the
// letters-only stream carries no absolute cipher positions), matching Ragbaby.
//
// SEED = per-column monogram warm start: rank the ciphertext letters by stream frequency and map
// the r-th most frequent ciphertext letter to the r-th most frequent English letter, then apply a
// few random swaps for restart diversity. The anneal only has to correct a handful of cells.

typedef struct {
    // Problem instance. `orig` is the full decoded ciphertext (letters + sentinels), kept so the
    // Aristocrat report can reinsert the original word divisions around the recovered letters.
    int    *orig;
    int     orig_len;
    int     n_letters;           // == ctx->cipher_len (the letters-only search length)

    // Warm-start ranks (built once): ciphertext letters by descending stream frequency, and
    // English letters by descending monogram weight. seed maps rank-for-rank.
    int     ct_rank[ALPHABET_SIZE];
    int     eng_rank[ALPHABET_SIZE];

    // --- incremental-scoring caches (homophonic fast path, minus counts/chi2/cribs) ---
    double  scale;               // n-gram scale factor (matches ngram_score())
    double  ngsum;               // running raw sum of ngram_data[] over all windows
    int    *pos;                 // length n_letters: cipher positions grouped by ciphertext letter
    int    *pos_off;             // length ALPHABET_SIZE + 1
    char   *win_mark;            // length n_letters (touched-window set scratch; stays all-zero)
    int    *win_list;            // length n_letters
    double  pend_ngsum;          // delta stashed by score_neighbor, applied by commit_neighbor
    int     pend_nsym;
    int    *pend_sym;            // length ALPHABET_SIZE: changed ciphertext letters
    int    *pend_newc;           // length ALPHABET_SIZE: their new plaintext letters
} AristocratScratch;

// Packed big-endian base-g_alpha index of the n-gram window at start `w`, reading each position's
// plaintext letter through `letter(p)`. Matches ngram_score()'s packing exactly (idx = sum_j
// letter(w+j) * g_alpha^(ng-1-j)), so the incremental ngsum stays bit-identical to a full rescore.
#define ARIST_WINDOW_INDEX(w, ng, letter) ({                \
    int _idx = 0;                                           \
    for (int _j = 0; _j < (ng); _j++)                       \
        _idx = _idx * g_alpha + letter((w) + _j);           \
    _idx; })

// ---------------------------------------------------------------------
//  Move (also used by the solver-test move-invariant check)
// ---------------------------------------------------------------------

// Swap two ciphertext letters' plaintext images: keeps key[] a permutation of 0..alpha-1.
void aristocrat_move(int *key, int alpha) {
    int a = rand_int(0, alpha), b = rand_int(0, alpha);
    int t = key[a]; key[a] = key[b]; key[b] = t;
}

// ---------------------------------------------------------------------
//  CipherModel hooks
// ---------------------------------------------------------------------

static int aristocrat_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    (void) ctx;
    if (cap < 1) return 0;
    out[0].period = g_alpha; out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Frequency rank-match warm start + a few random swaps for restart diversity.
static void aristocrat_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    const AristocratScratch *a = (const AristocratScratch *) ctx->model_scratch;
    for (int r = 0; r < g_alpha; r++) st->key[a->ct_rank[r]] = a->eng_rank[r];
    st->key_len = g_alpha;
    int kicks = rand_int(0, 9);
    for (int k = 0; k < kicks; k++) aristocrat_move(st->key, g_alpha);
}

static void aristocrat_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                               SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    aristocrat_move(st->key, g_alpha);
    if (frand() < 0.10) aristocrat_move(st->key, g_alpha);   // occasional double swap (exploration)
}

static void aristocrat_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < g_alpha; i++) dst->key[i] = src->key[i];
    dst->key_len = src->key_len;
}

static void aristocrat_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                    SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    for (int i = 0; i < ctx->cipher_len; i++) out[i] = st->key[ctx->cipher[i]];
    *score_adjust = 0.0;         // a bijection cannot collapse letters -> no penalty term
}

// ---------------------------------------------------------------------
//  Incremental fast path (homophonic twin, n-gram only)
// ---------------------------------------------------------------------

static void aristocrat_sync_caches(const SolverCtx *ctx, const SolverConfig *cc, const int *dec) {
    (void) cc;
    AristocratScratch *a = (AristocratScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, ng = ctx->cfg->ngram_size;
    a->ngsum = 0.0;
    int n_windows = len - ng + 1;
    #define dec_at(q) (dec[q])
    for (int w = 0; w < n_windows; w++)
        a->ngsum += ctx->ngram_data[ARIST_WINDOW_INDEX(w, ng, dec_at)];
    #undef dec_at
}

static double aristocrat_score_neighbor(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *cur, const SolverState *loc, const int *cur_dec, double cur_score) {
    (void) cur_score;
    AristocratScratch *a = (AristocratScratch *) ctx->model_scratch;
    ColossusConfig *cfg = ctx->cfg;
    int len = ctx->cipher_len, ng = cfg->ngram_size, n = cc->period;
    const int *cipher = ctx->cipher;

    // 1. Which ciphertext letters changed image (a swap touches exactly two).
    a->pend_nsym = 0;
    for (int s = 0; s < n; s++) {
        if (loc->key[s] == cur->key[s]) continue;
        a->pend_sym[a->pend_nsym]  = s;
        a->pend_newc[a->pend_nsym] = loc->key[s];
        a->pend_nsym++;
    }

    // 2. n-gram delta over the windows touching any changed position.
    double dsum = 0.0;
    int n_windows = len - ng + 1;
    if (cfg->weight_ngram > 1.e-4 && n_windows > 0) {
        int nlist = 0;
        for (int k = 0; k < a->pend_nsym; k++) {
            int s = a->pend_sym[k];
            for (int pi = a->pos_off[s]; pi < a->pos_off[s + 1]; pi++) {
                int p = a->pos[pi];
                int lo = p - ng + 1; if (lo < 0) lo = 0;
                int hi = p;          if (hi > n_windows - 1) hi = n_windows - 1;
                for (int w = lo; w <= hi; w++)
                    if (!a->win_mark[w]) { a->win_mark[w] = 1; a->win_list[nlist++] = w; }
            }
        }
        #define dec_at(q) (cur_dec[q])
        #define loc_at(q) (loc->key[cipher[q]])
        for (int i = 0; i < nlist; i++) {
            int w = a->win_list[i];
            int old_idx = ARIST_WINDOW_INDEX(w, ng, dec_at);
            int new_idx = ARIST_WINDOW_INDEX(w, ng, loc_at);
            dsum += ctx->ngram_data[new_idx] - ctx->ngram_data[old_idx];
            a->win_mark[w] = 0;            // clear for next call
        }
        #undef loc_at
        #undef dec_at
    }
    a->pend_ngsum = a->ngsum + dsum;

    // 3. Reassemble the score exactly as ngram_score (n-gram only; no cribs, no penalty).
    if (cfg->weight_ngram > 1.e-4 && n_windows > 0)
        return a->scale * a->pend_ngsum / (len - ng);
    return 0.0;
}

static void aristocrat_commit_neighbor(const SolverCtx *ctx, const SolverConfig *cc, int *cur_dec) {
    (void) cc;
    AristocratScratch *a = (AristocratScratch *) ctx->model_scratch;
    for (int k = 0; k < a->pend_nsym; k++) {
        int s = a->pend_sym[k], c = a->pend_newc[k];
        for (int pi = a->pos_off[s]; pi < a->pos_off[s + 1]; pi++)
            cur_dec[a->pos[pi]] = c;
    }
    a->ngsum = a->pend_ngsum;
}

static void *aristocrat_scratch_clone(const SolverCtx *ctx) {
    const AristocratScratch *src = (const AristocratScratch *) ctx->model_scratch;
    int len = src->n_letters;
    AristocratScratch *a = malloc(sizeof *a);
    if (!a) return NULL;
    memcpy(a, src, sizeof *a);      // share read-only fields (orig, ranks, scale, pos/pos_off)
    // Private buffers the search mutates.
    a->win_mark  = calloc(len > 0 ? len : 1, sizeof(char));
    a->win_list  = malloc(sizeof(int) * (len > 0 ? len : 1));
    a->pend_sym  = malloc(sizeof(int) * ALPHABET_SIZE);
    a->pend_newc = malloc(sizeof(int) * ALPHABET_SIZE);
    if (!a->win_mark || !a->win_list || !a->pend_sym || !a->pend_newc) {
        free(a->win_mark); free(a->win_list); free(a->pend_sym); free(a->pend_newc);
        free(a);
        return NULL;
    }
    return a;
}

static void aristocrat_scratch_free(void *scratch) {
    AristocratScratch *a = (AristocratScratch *) scratch;
    if (!a) return;
    free(a->win_mark); free(a->win_list); free(a->pend_sym); free(a->pend_newc);
    free(a);
}

// ---------------------------------------------------------------------
//  Reporting
// ---------------------------------------------------------------------

// Recovered ciphertext alphabet (the row under plaintext A..Z): ct = key_inv[pt].
static void aristocrat_alpha_string(const int key[], char out[]) {
    int inv[ALPHABET_SIZE];
    for (int ct = 0; ct < g_alpha; ct++) inv[key[ct]] = ct;
    for (int p = 0; p < g_alpha; p++) out[p] = index_to_char(inv[p]);
    out[g_alpha] = '\0';
}

// Aristocrat: walk the original layout, emitting a space/punctuation sentinel verbatim and each
// letter as its recovered plaintext letter -- reproducing the ciphertext's word divisions.
static void aristocrat_print_spaced(const int orig[], int orig_len, const int decrypted[]) {
    int c = 0;
    for (int i = 0; i < orig_len; i++)
        putchar(orig[i] >= 0 ? index_to_char(decrypted[c++]) : index_to_char(orig[i]));
}

// Patristocrat: the recovered letters in 5-letter groups.
static void aristocrat_print_grouped(const int decrypted[], int len) {
    for (int i = 0; i < len; i++) {
        if (i > 0 && i % 5 == 0) putchar(' ');
        putchar(index_to_char(decrypted[i]));
    }
}

static void aristocrat_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) st;
    report_transposition_verbose(ctx, score, decrypted, stats, "monoalphabetic");
}

static void aristocrat_report(const SolverCtx *ctx, const SolverConfig *cc,
                              const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const AristocratScratch *a = (const AristocratScratch *) ctx->model_scratch;
    int len = ctx->cipher_len;
    bool patristocrat = (cfg->cipher_type == PATRISTOCRAT);

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char alpha[ALPHABET_SIZE + 1]; aristocrat_alpha_string(st->key, alpha);

    printf("\nResult Score: %.2f | Words: %d | ct-alphabet=%s\n", score, n_words_found, alpha);

    // Cipher and plaintext, laid out the way this variant presents them.
    if (patristocrat) {
        // Cipher: the letters-only search stream in 5-letter groups; plaintext below it.
        aristocrat_print_grouped(ctx->cipher, len); printf("\n");
        aristocrat_print_grouped(decrypted, len);   printf("\n");
    } else {
        // Cipher and plaintext with the original word divisions restored.
        aristocrat_print_spaced(a->orig, a->orig_len, ctx->cipher); printf("\n");
        aristocrat_print_spaced(a->orig, a->orig_len, decrypted);   printf("\n");
    }
    printf("%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, ct-alphabet=, file, CIPHER, PLAINTEXT
    // (letters-only, so the accuracy suite's whitespace-stripped comparison matches).
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, ct-alphabet=%s, ", score, n_words_found, cfg->cipher_type, alpha);
    else
        printf(">>> %.2f, %d, ct-alphabet=%s, ", score, cfg->cipher_type, alpha);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(ctx->cipher, len);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel ARISTOCRAT_MODEL = {
    .name = "aristocrat", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = aristocrat_enumerate, .key_len = NULL,
    .seed = aristocrat_seed, .perturb = aristocrat_perturb, .copy_state = aristocrat_copy,
    .decrypt = aristocrat_decrypt_hook, .report = aristocrat_report,
    .report_verbose = aristocrat_report_verbose,
    .sync_caches = aristocrat_sync_caches,
    .score_neighbor = aristocrat_score_neighbor,
    .commit_neighbor = aristocrat_commit_neighbor,
    .scratch_clone = aristocrat_scratch_clone,
    .scratch_free = aristocrat_scratch_free,
};

// ---------------------------------------------------------------------
//  Entry point
// ---------------------------------------------------------------------

// The letters-only search stream + the original layout, off the stack.
static int g_arist_letters[MAX_CIPHER_LENGTH];
static int g_arist_orig[MAX_CIPHER_LENGTH];

void solve_aristocrat(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Derive the letters-only stream from the decoded ciphertext; keep the full layout so the
    // Aristocrat report can restore the word divisions. Sentinels (spaces/punctuation) are dropped
    // from the search (they are transparent to scoring anyway) and, for the Patristocrat, are just
    // the meaningless 5-letter grouping.
    int nl = 0;
    for (int i = 0; i < cipher_len; i++) {
        g_arist_orig[i] = cipher_indices[i];
        if (cipher_indices[i] >= 0) g_arist_letters[nl++] = cipher_indices[i];
    }
    if (nl < 4) {
        printf("\n\nERROR: ciphertext too short for an %s solve (%d letters).\n\n",
            cfg->cipher_type == PATRISTOCRAT ? "Patristocrat" : "Aristocrat", nl);
        return;
    }

    if (cfg->verbose)
        printf("\n%s: %d letters, 26-permutation anneal (n-gram only)\n",
            cfg->cipher_type == PATRISTOCRAT ? "patristocrat" : "aristocrat", nl);

    // Cribs are not used (the letters-only stream carries no absolute cipher positions).
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        g_arist_letters, nl, crib_indices, crib_positions, 0);
    ctx.result = result;

    AristocratScratch scratch;
    memset(&scratch, 0, sizeof scratch);
    scratch.orig = g_arist_orig;
    scratch.orig_len = cipher_len;
    scratch.n_letters = nl;
    scratch.scale = g_ngram_logprob ? 1.0 : pow((double) g_alpha, cfg->ngram_size);

    // English letters by descending monogram weight (warm-start target ranking).
    for (int i = 0; i < g_alpha; i++) scratch.eng_rank[i] = i;
    for (int x = 0; x < g_alpha; x++)
        for (int y = x + 1; y < g_alpha; y++)
            if (g_monograms[scratch.eng_rank[y]] > g_monograms[scratch.eng_rank[x]]) {
                int t = scratch.eng_rank[x]; scratch.eng_rank[x] = scratch.eng_rank[y]; scratch.eng_rank[y] = t;
            }
    // Ciphertext letters by descending stream frequency.
    int hist[ALPHABET_SIZE];
    for (int c = 0; c < g_alpha; c++) { hist[c] = 0; scratch.ct_rank[c] = c; }
    for (int i = 0; i < nl; i++) hist[g_arist_letters[i]]++;
    for (int x = 0; x < g_alpha; x++)
        for (int y = x + 1; y < g_alpha; y++)
            if (hist[scratch.ct_rank[y]] > hist[scratch.ct_rank[x]]) {
                int t = scratch.ct_rank[x]; scratch.ct_rank[x] = scratch.ct_rank[y]; scratch.ct_rank[y] = t;
            }

    // Position index: bucket each cipher position by its ciphertext letter so a swap maps straight
    // to the positions (and thus n-gram windows) it changes. Built once (counting sort).
    scratch.pos     = malloc(sizeof(int) * (nl > 0 ? nl : 1));
    scratch.pos_off = malloc(sizeof(int) * (ALPHABET_SIZE + 1));
    scratch.win_mark = calloc(nl > 0 ? nl : 1, sizeof(char));
    scratch.win_list = malloc(sizeof(int) * (nl > 0 ? nl : 1));
    scratch.pend_sym  = malloc(sizeof(int) * ALPHABET_SIZE);
    scratch.pend_newc = malloc(sizeof(int) * ALPHABET_SIZE);
    if (!scratch.pos || !scratch.pos_off || !scratch.win_mark ||
        !scratch.win_list || !scratch.pend_sym || !scratch.pend_newc) {
        printf("\n\nERROR: out of memory in aristocrat solve.\n\n");
        free(scratch.pos); free(scratch.pos_off); free(scratch.win_mark);
        free(scratch.win_list); free(scratch.pend_sym); free(scratch.pend_newc);
        return;
    }
    for (int s = 0; s <= g_alpha; s++) scratch.pos_off[s] = 0;
    for (int i = 0; i < nl; i++) scratch.pos_off[g_arist_letters[i] + 1]++;
    for (int s = 0; s < g_alpha; s++) scratch.pos_off[s + 1] += scratch.pos_off[s];
    {
        int cursor[ALPHABET_SIZE];
        for (int s = 0; s < g_alpha; s++) cursor[s] = scratch.pos_off[s];
        for (int i = 0; i < nl; i++) scratch.pos[cursor[g_arist_letters[i]]++] = i;
    }

    ctx.model_scratch = &scratch;
    run_solver(&ARISTOCRAT_MODEL, &ctx);

    free(scratch.pos); free(scratch.pos_off); free(scratch.win_mark);
    free(scratch.win_list); free(scratch.pend_sym); free(scratch.pend_newc);
}
