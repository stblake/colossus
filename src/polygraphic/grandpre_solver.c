#include "grandpre_solver.h"
#include "grandpre.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Grandpre solver (TYPE grandpre) -- see grandpre_solver.h
// =====================================================================
//
// Grandpre is a HOMOPHONIC substitution over 2-digit numeric codes: decoding a code is a
// unique lookup into the N x N square, but a plaintext letter has one code per cell it
// occupies. The solver parses the digit stream into 2-digit codes, interns the distinct
// codes into symbol ids, and hill-climbs the many-to-one map key[symbol] -> letter exactly
// as homophonic_solver.c does -- same monogram-flattening seed, same single-symbol
// reassignment move set, same chi-squared anti-collapse penalty (a homophonic map is free
// to fold many codes onto E/T/A to tile frequent n-grams), and the same incremental
// score_neighbor / commit_neighbor / per-thread scratch_clone fast path (each move rescans
// only the n-gram windows its changed codes touch). This file is a self-contained twin of
// the homophonic model; the only Grandpre-specific parts are the code parsing (below) and
// the report, which lays the recovered map back into the N x N square.
//
// CAPABILITY. A ~64-code (8x8) homophonic map over the ACA 150-200 letter range is
// undersampled (~2-3 samples/code), so blind recovery is partial and high-variance at the
// short end and reliable only at longer lengths -- characterized (a length/accuracy curve),
// not asserted at 99%, in test_grandpre_solver.c. Best with -logprob (+ quintgrams).

#define GP_MAX_CODES 100                 // <= GRANDPRE_MAX_GRID distinct codes (10x10)

// Parsed ciphertext: the interned symbol-id stream + the symbol->code map, filled by
// solve_grandpre and read-only during the search.
static int g_gp_stream[MAX_CIPHER_LENGTH];   // position -> symbol id (0..n_symbols-1)

typedef struct {
    int    n_symbols;                // distinct observed codes
    int    N;                        // inferred square side (for the report grid)
    int    code_of_sym[GP_MAX_CODES];// symbol id -> its 2-digit code

    // --- incremental-scoring caches (homophonic fast path; identical semantics) ---
    double  scale;                   // n-gram scale factor (matches ngram_score())
    double  ngsum;                   // running raw sum of ngram_data[] over all windows
    int     counts[ALPHABET_SIZE];   // letter histogram of the current decryption
    int    *pos;                     // length cipher_len: positions grouped by symbol
    int    *pos_off;                 // length n_symbols + 1
    char   *win_mark;                // length cipher_len (touched-window scratch, stays 0)
    int    *win_list;                // length cipher_len
    double  pend_ngsum;
    int     pend_counts[ALPHABET_SIZE];
    int     pend_nsym;
    int    *pend_sym;                // length n_symbols: changed symbols
    int    *pend_newc;               // length n_symbols: their new plaintext letters
} GrandpreScratch;

// One config: the whole map is climbed at once; period carries the symbol count.
static int gp_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const GrandpreScratch *h = (const GrandpreScratch *) ctx->model_scratch;
    if (cap < 1) return 0;
    out[0].period = h->n_symbols;
    out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Frequency-flattening seed: draw each symbol's letter from the English monogram
// distribution (common letters naturally get more homophones).
static void gp_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
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

static double gp_penalty(const SolverCtx *ctx, const int *dec) {
    if (ctx->cfg->weight_monogram <= 1.e-9) return 0.0;
    return ctx->cfg->weight_monogram * chi_squared((int *) dec, ctx->cipher_len);
}

// Neighbour move: dominant single random reassignment + a symbol-pair swap + a whole
// letter-class swap (identical to homophonic_perturb).
static void gp_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                       SolverState *st, bool *force_primary) {
    (void) ctx; (void) force_primary;
    int n = cc->period;
    if (n < 1) return;
    double r = frand();
    if (r < 0.85) {
        st->key[rand_int(0, n)] = rand_int(0, g_alpha);
    } else if (n >= 2 && r < 0.93) {
        int a = rand_int(0, n), b = rand_int(0, n);
        int t = st->key[a]; st->key[a] = st->key[b]; st->key[b] = t;
    } else {
        int a = rand_int(0, g_alpha), b = rand_int(0, g_alpha);
        if (a != b)
            for (int s = 0; s < n; s++) {
                if (st->key[s] == a) st->key[s] = b;
                else if (st->key[s] == b) st->key[s] = a;
            }
    }
}

// Packed big-endian base-g_alpha n-gram window index (matches ngram_score()).
#define GP_WINDOW_INDEX(w, ng, letter) ({                   \
    int _idx = 0;                                           \
    for (int _j = 0; _j < (ng); _j++)                       \
        _idx = _idx * g_alpha + letter((w) + _j);           \
    _idx; })

static void gp_sync_caches(const SolverCtx *ctx, const SolverConfig *cc, const int *dec) {
    (void) cc;
    GrandpreScratch *h = (GrandpreScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, ng = ctx->cfg->ngram_size;
    tally((int *) dec, len, h->counts, ALPHABET_SIZE);
    h->ngsum = 0.0;
    int n_windows = len - ng + 1;
    #define dec_at(q) (dec[q])
    for (int w = 0; w < n_windows; w++)
        h->ngsum += ctx->ngram_data[GP_WINDOW_INDEX(w, ng, dec_at)];
    #undef dec_at
}

static double gp_chi2_from_counts(const int *counts, int len) {
    double chi2 = 0.0;
    for (int c = 0; c < g_alpha; c++) {
        double f = ((double) counts[c]) / len;
        chi2 += pow(f - g_monograms[c], 2) / g_monograms[c];
    }
    return chi2;
}

static double gp_score_neighbor(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *cur, const SolverState *loc, const int *cur_dec, double cur_score) {
    (void) cur_score;
    GrandpreScratch *h = (GrandpreScratch *) ctx->model_scratch;
    ColossusConfig *cfg = ctx->cfg;
    int len = ctx->cipher_len, ng = cfg->ngram_size, n = cc->period;
    const int *cipher = ctx->cipher;

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
        #define dec_at(q) (cur_dec[q])
        #define loc_at(q) (loc->key[cipher[q]])
        for (int i = 0; i < nlist; i++) {
            int w = h->win_list[i];
            int old_idx = GP_WINDOW_INDEX(w, ng, dec_at);
            int new_idx = GP_WINDOW_INDEX(w, ng, loc_at);
            dsum += ctx->ngram_data[new_idx] - ctx->ngram_data[old_idx];
            h->win_mark[w] = 0;
        }
        #undef loc_at
        #undef dec_at
    }
    h->pend_ngsum = h->ngsum + dsum;

    double ng_score = 0.0;
    if (cfg->weight_ngram > 1.e-4 && n_windows > 0)
        ng_score = h->scale * h->pend_ngsum / (len - ng);

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
        score -= cfg->weight_monogram * gp_chi2_from_counts(h->pend_counts, len);

    return score;
}

static void gp_commit_neighbor(const SolverCtx *ctx, const SolverConfig *cc, int *cur_dec) {
    (void) cc;
    GrandpreScratch *h = (GrandpreScratch *) ctx->model_scratch;
    for (int k = 0; k < h->pend_nsym; k++) {
        int s = h->pend_sym[k], c = h->pend_newc[k];
        for (int pi = h->pos_off[s]; pi < h->pos_off[s + 1]; pi++)
            cur_dec[h->pos[pi]] = c;
    }
    for (int c = 0; c < g_alpha; c++) h->counts[c] = h->pend_counts[c];
    h->ngsum = h->pend_ngsum;
}

static void gp_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}

static void *gp_scratch_clone(const SolverCtx *ctx) {
    const GrandpreScratch *src = (const GrandpreScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, n = src->n_symbols;
    GrandpreScratch *h = malloc(sizeof *h);
    if (!h) return NULL;
    memcpy(h, src, sizeof *h);          // share read-only fields (code map, scale, pos index)
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

static void gp_scratch_free(void *scratch) {
    GrandpreScratch *h = (GrandpreScratch *) scratch;
    if (!h) return;
    free(h->win_mark); free(h->win_list); free(h->pend_sym); free(h->pend_newc);
    free(h);
}

static void gp_decrypt(const SolverCtx *ctx, const SolverConfig *cc,
                       SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    for (int i = 0; i < ctx->cipher_len; i++) out[i] = st->key[ctx->cipher[i]];
    *score_adjust = -gp_penalty(ctx, out);
}

static void gp_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) st;
    const GrandpreScratch *h = (const GrandpreScratch *) ctx->model_scratch;
    char params[64];
    snprintf(params, sizeof(params), "%dx%d codes=%d", h->N, h->N, h->n_symbols);
    report_transposition_verbose(ctx, score, decrypted, stats, params);
}

// Lay the recovered code -> letter map into the N x N square: cell (r, c) holds the letter
// of code grandpre_code(r, c, N) if that code was observed (else '?').
static void gp_print_square(const GrandpreScratch *h, const SolverState *st) {
    int N = h->N;
    int letter_of_code[100];
    for (int i = 0; i < 100; i++) letter_of_code[i] = -1;
    for (int s = 0; s < h->n_symbols; s++) letter_of_code[h->code_of_sym[s]] = st->key[s];
    printf("     ");
    for (int c = 0; c < N; c++) printf("%d ", grandpre_label(c, N));
    printf("\n");
    for (int r = 0; r < N; r++) {
        printf("  %d  ", grandpre_label(r, N));
        for (int c = 0; c < N; c++) {
            int L = letter_of_code[grandpre_code(r, c, N)];
            printf("%c ", (L >= 0) ? index_to_char(L) : '?');
        }
        printf("\n");
    }
}

static void gp_report(const SolverCtx *ctx, const SolverConfig *cc,
                      const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const GrandpreScratch *h = (const GrandpreScratch *) ctx->model_scratch;
    int n = h->n_symbols, len = ctx->cipher_len;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    printf("\nResult Score: %.2f | Words: %d | %dx%d square, %d codes\n",
        score, n_words_found, h->N, h->N, n);

    print_text(decrypted, len);
    printf("\n%s\n", ctx->cribtext);

    printf("\nrecovered %dx%d square (row/col labels are the grid indices):\n", h->N, h->N);
    gp_print_square(h, st);

    // Homophone key: each plaintext letter and the codes decoding to it.
    printf("code key (plaintext <- codes):\n");
    for (int c = 0; c < g_alpha; c++) {
        int any = 0;
        for (int s = 0; s < n; s++) if (st->key[s] == c) {
            if (!any) { printf("  %c <-", index_to_char(c)); any = 1; }
            printf(" %02d", h->code_of_sym[s]);
        }
        if (any) printf("\n");
    }

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner: >>> score, [words,] type, codes=N, file, PLAINTEXT (letters-only).
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, codes=%d, ", score, n_words_found, cfg->cipher_type, n);
    else
        printf(">>> %.2f, %d, codes=%d, ", score, cfg->cipher_type, n);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel GRANDPRE_MODEL = {
    .name = "grandpre", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = gp_enumerate, .key_len = NULL,
    .seed = gp_seed, .perturb = gp_perturb, .copy_state = gp_copy,
    .decrypt = gp_decrypt, .report = gp_report, .report_verbose = gp_report_verbose,
    .sync_caches = gp_sync_caches, .score_neighbor = gp_score_neighbor,
    .commit_neighbor = gp_commit_neighbor,
    .scratch_clone = gp_scratch_clone, .scratch_free = gp_scratch_free,
};

// Parse the ciphertext digit stream into 2-digit codes; return the code count.
static int gp_parse_codes(const char *s, int out[], int cap) {
    int digits[2 * MAX_CIPHER_LENGTH], nd = 0;
    for (int i = 0; s[i] && nd < 2 * MAX_CIPHER_LENGTH; i++)
        if (s[i] >= '0' && s[i] <= '9') digits[nd++] = s[i] - '0';
    int n = 0;
    for (int i = 0; i + 1 < nd && n < cap; i += 2) out[n++] = digits[i] * 10 + digits[i + 1];
    return n;
}

void solve_grandpre(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_indices; (void) cipher_len;
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != DEFAULT_ALPHABET_SIZE) {
        printf("\n\nERROR: Grandpre needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }

    static int codes[MAX_CIPHER_LENGTH];
    int len = gp_parse_codes(ciphertext_str, codes, MAX_CIPHER_LENGTH);
    if (len < 4) {
        printf("\n\nERROR: parsed only %d Grandpre codes; need >= 4. The ciphertext must be a "
               "stream of 2-digit numbers.\n\n", len);
        return;
    }

    // Intern the distinct codes into symbol ids; infer the square side from the digit range
    // (a '0' digit => 10x10 with labels 0..9, else N = the largest digit seen).
    int sym_of_code[100];
    for (int i = 0; i < 100; i++) sym_of_code[i] = -1;
    GrandpreScratch scratch;
    memset(&scratch, 0, sizeof scratch);
    int n_sym = 0, maxd = 0, has_zero = 0;
    for (int i = 0; i < len; i++) {
        int code = codes[i];
        int d0 = code / 10, d1 = code % 10;
        if (d0 > maxd) maxd = d0;
        if (d1 > maxd) maxd = d1;
        if (d0 == 0 || d1 == 0) has_zero = 1;
        if (sym_of_code[code] < 0) {
            if (n_sym >= GP_MAX_CODES) { printf("\n\nERROR: too many distinct codes.\n\n"); return; }
            sym_of_code[code] = n_sym;
            scratch.code_of_sym[n_sym] = code;
            n_sym++;
        }
        g_gp_stream[i] = sym_of_code[code];
    }
    scratch.N = has_zero ? 10 : (maxd < GRANDPRE_MIN_SIDE ? GRANDPRE_MIN_SIDE : maxd);
    scratch.n_symbols = n_sym;

    if (cfg->verbose)
        printf("\ngrandpre: %d codes, %d distinct, inferred %dx%d square\n",
               len, n_sym, scratch.N, scratch.N);

    scratch.scale = g_ngram_logprob ? 1.0 : pow((double) g_alpha, cfg->ngram_size);

    // Cribs are not used: main() indexes crib positions against the letter-decoded ciphertext,
    // which does not line up with the interned 2-digit-code stream (like the other digit-stream
    // types). Pass 0 to avoid an out-of-range crib index.
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        g_gp_stream, len, NULL, NULL, 0);
    ctx.result = result;

    // Position index: bucket positions by symbol for the incremental fast path.
    scratch.pos     = malloc(sizeof(int) * (len > 0 ? len : 1));
    scratch.pos_off = malloc(sizeof(int) * (n_sym + 1));
    scratch.win_mark = calloc(len > 0 ? len : 1, sizeof(char));
    scratch.win_list = malloc(sizeof(int) * (len > 0 ? len : 1));
    scratch.pend_sym  = malloc(sizeof(int) * (n_sym > 0 ? n_sym : 1));
    scratch.pend_newc = malloc(sizeof(int) * (n_sym > 0 ? n_sym : 1));
    if (!scratch.pos || !scratch.pos_off || !scratch.win_mark ||
        !scratch.win_list || !scratch.pend_sym || !scratch.pend_newc) {
        printf("\n\nERROR: out of memory in grandpre solve.\n\n");
        free(scratch.pos); free(scratch.pos_off); free(scratch.win_mark);
        free(scratch.win_list); free(scratch.pend_sym); free(scratch.pend_newc);
        return;
    }
    for (int s = 0; s <= n_sym; s++) scratch.pos_off[s] = 0;
    for (int i = 0; i < len; i++) scratch.pos_off[g_gp_stream[i] + 1]++;
    for (int s = 0; s < n_sym; s++) scratch.pos_off[s + 1] += scratch.pos_off[s];
    {
        int *cursor = malloc(sizeof(int) * (n_sym > 0 ? n_sym : 1));
        for (int s = 0; s < n_sym; s++) cursor[s] = scratch.pos_off[s];
        for (int i = 0; i < len; i++) scratch.pos[cursor[g_gp_stream[i]]++] = i;
        free(cursor);
    }

    ctx.model_scratch = &scratch;
    run_solver(&GRANDPRE_MODEL, &ctx);

    free(scratch.pos); free(scratch.pos_off); free(scratch.win_mark);
    free(scratch.win_list); free(scratch.pend_sym); free(scratch.pend_newc);
}
