#include "running_key_solver.h"
#include "running_key.h"
#include "engine.h"
#include "scoring.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// =====================================================================
//  Running Key solver  (TYPE running-key / runningkey / rk)
// =====================================================================
//
// A Running Key enciphers a message against a key stream as long as the message, so it
// has no period -- the whole key is the unknown. See running_key.c for the primitive
// and running_key_solver.h for the three modes. The cryptanalytic crux (blind) is that
// fixing the key stream K determines the plaintext P[i] = decode(CT[i], K[i]); scoring
// ONLY P is hopeless (any target P is producible by some gibberish K), so the objective
// scores BOTH streams as English at once. That joint term is folded in through the
// engine's *score_adjust hook (state_score(P) supplies the plaintext half, the hook
// adds the key half), so the engine needs no change.
//
// State: SolverState.key[0..N-1] holds the running-key stream K (key_len = N). One
// CipherModel serves all three modes:
//   * blind (self-keyed / independent): key_len() returns N so the engine climbs K;
//     each restart is SEEDED from a deterministic beam warm start (below) plus a random
//     jitter, then annealed (SHAPE_ANNEAL). Self-keyed adds a seam reward at the K|P
//     join (the two halves are one contiguous passage) and reports the full K||P (2N).
//   * known-key (-runningkeyfile): key_len() returns 0 so each family config is a single
//     SWEEP candidate -- seed loads the known key, decrypt+score once, best family kept.
//
// The four families (Vigenere/Beaufort/Variant/Porta) are swept as engine configs; the
// n-gram picks the winner (like a period sweep). -variant / -beaufort pin one family.
//
// Warm start (the capability + efficiency lever): 26^N is far too large for a cold
// climb, but the additive coupling gives a strong left-to-right signal. rk_beam_warmstart
// runs a beam over positions 0..N-1: state = the key prefix, each of the 26 extensions
// scored by the n-gram window (in BOTH streams) it closes; top-B kept. This descends from
// the "deterministic warm start + n-gram climb" pattern (derive_optimal_cycleword's
// per-column monogram warm start does not transfer -- one char per column when the period
// equals the length -- so a DP replaces it). The beam is O(N^2 * B) in the prefix copies,
// which is fine because Running Keys are short by construction (ACA 40-50 plaintext).

#define RK_BEAM_DEFAULT   96      // beam width for the warm start
#define RK_BEAM_MAX      400
#define RK_BEAM_BUDGET (192 * 1024)  // cap B*N (bounds warm-start memory + O(N^2 B) time)
#define RK_JITTER_MAX    0.18     // max fraction of warm-start positions randomised per restart
#define RK_MAXNG           8      // n-gram order clamp for the window buffers

// Per-config invariant instance handed to the hooks via ctx->model_scratch.
typedef struct {
    int   n;                              // cipher length N
    bool  self_keyed;                     // add the seam term + report K||P (2N)
    bool  known_key;                      // -runningkeyfile: SWEEP, no climb
    int   n_families;                     // number of family configs
    int   families[RK_N_FAMILIES];        // family value per config index
    int  *warm;                           // [n_families * N] per-family warm-start keys (blind)
    int  *known_key_stream;               // [N] loaded key (known-key mode), else NULL
    // Crib ANCHORING. A known plaintext letter at a ciphertext position FORCES the running-key
    // letter there (rk_key_from_pt), so crib positions are pinned to their exact key value and
    // held fixed while the search fills the free positions. This is far stronger than a soft
    // crib reward (which the gaming basins swamp). crib_pt[p] = the crib plaintext letter at pt
    // position p, or -1 if free; free_pos[0..n_free-1] lists the searchable positions.
    int  *crib_pt;                        // [N], letter 0..25 at pinned positions, else -1
    int  *free_pos;                       // [N], the n_free searchable positions
    int   n_free;
    int   n_cribs;
} RkScratch;

// ===================================================================
//  Beam warm start
// ===================================================================

double rk_beam_warmstart(const int cipher[], int n, int family,
                         const float *ngram, int ngram_size, int beam_width, int key_out[]) {
    if (n <= 0) return 0.0;
    int ng = ngram_size; if (ng < 1) ng = 1; if (ng > RK_MAXNG) ng = RK_MAXNG;

    int B = beam_width;
    if (B < 26) B = 26;
    if (B > RK_BEAM_MAX) B = RK_BEAM_MAX;
    if ((long) B * n > RK_BEAM_BUDGET) { B = RK_BEAM_BUDGET / n; if (B < 26) B = 26; }

    int   *cur_key = (int *)    malloc((size_t) B * n * sizeof(int));
    int   *nxt_key = (int *)    malloc((size_t) B * n * sizeof(int));
    double *cur_sc = (double *) malloc((size_t) B * sizeof(double));
    double *nxt_sc = (double *) malloc((size_t) B * sizeof(double));
    struct { double score; int s, g; } *kept = malloc((size_t) B * sizeof(*kept));
    if (!cur_key || !nxt_key || !cur_sc || !nxt_sc || !kept) {
        free(cur_key); free(nxt_key); free(cur_sc); free(nxt_sc); free(kept);
        for (int i = 0; i < n; i++) key_out[i] = 0;         // degrade gracefully
        return -1e300;
    }

    int ncur = 1;                                            // one empty prefix
    cur_sc[0] = 0.0;

    for (int i = 0; i < n; i++) {
        int nkept = 0, worst = 0;
        for (int s = 0; s < ncur; s++) {
            for (int g = 0; g < g_alpha; g++) {
                double delta = 0.0;
                if (i >= ng - 1) {
                    int kwin[RK_MAXNG], pwin[RK_MAXNG];
                    for (int j = 0; j < ng; j++) {
                        int pos = i - ng + 1 + j;
                        int kk = (pos == i) ? g : cur_key[(size_t) s * n + pos];
                        kwin[j] = kk;
                        pwin[j] = rk_decode_char(cipher[pos], kk, family);
                    }
                    delta = ngram_weight_at(ngram, ngram_index_int(kwin, ng))
                          + ngram_weight_at(ngram, ngram_index_int(pwin, ng));
                }
                double sc = cur_sc[s] + delta;
                if (nkept < B) {
                    kept[nkept].score = sc; kept[nkept].s = s; kept[nkept].g = g;
                    if (nkept == 0 || sc < kept[worst].score) worst = nkept;
                    nkept++;
                } else if (sc > kept[worst].score) {
                    kept[worst].score = sc; kept[worst].s = s; kept[worst].g = g;
                    worst = 0;
                    for (int t = 1; t < nkept; t++) if (kept[t].score < kept[worst].score) worst = t;
                }
            }
        }
        for (int t = 0; t < nkept; t++) {
            int s = kept[t].s, g = kept[t].g;
            if (i > 0) memcpy(&nxt_key[(size_t) t * n], &cur_key[(size_t) s * n], (size_t) i * sizeof(int));
            nxt_key[(size_t) t * n + i] = g;
            nxt_sc[t] = kept[t].score;
        }
        int *tk = cur_key; cur_key = nxt_key; nxt_key = tk;
        double *ts = cur_sc; cur_sc = nxt_sc; nxt_sc = ts;
        ncur = nkept;
    }

    int best = 0;
    for (int t = 1; t < ncur; t++) if (cur_sc[t] > cur_sc[best]) best = t;
    double best_score = cur_sc[best];
    memcpy(key_out, &cur_key[(size_t) best * n], (size_t) n * sizeof(int));

    free(cur_key); free(nxt_key); free(cur_sc); free(nxt_sc); free(kept);
    return best_score;
}

// ===================================================================
//  Scoring helper: seam reward at the K|P join (self-keyed only)
// ===================================================================

// n-gram windows straddling the boundary of the contiguous 2N passage T = K ++ P,
// normalised onto ngram_score's per-window frame so it composes with the two half-stream
// scores. Only (ngram_size-1) windows, so this is cheap.
static double rk_seam_adjust(const int *key, const int *pt, int n, int ngram_size,
                             const float *nd) {
    if (n < ngram_size || ngram_size < 2) return 0.0;
    int win[RK_MAXNG];
    int ng = ngram_size > RK_MAXNG ? RK_MAXNG : ngram_size;
    double sum = 0.0;
    for (int start = n - ng + 1; start <= n - 1; start++) {
        for (int j = 0; j < ng; j++) {
            int pos = start + j;
            win[j] = (pos < n) ? key[pos] : pt[pos - n];
        }
        sum += ngram_weight_at(nd, ngram_index_int(win, ng));
    }
    double scale = g_ngram_logprob ? 1.0 : pow((double) g_alpha, ng);
    int denom = n - ng;
    return (denom > 0) ? scale * sum / denom : 0.0;
}

// ===================================================================
//  CipherModel hooks
// ===================================================================

static int rk_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const RkScratch *a = (const RkScratch *) ctx->model_scratch;
    int nf = a->n_families;
    if (nf > cap) nf = cap;
    for (int f = 0; f < nf; f++) {
        out[f].period = a->n;
        out[f].j = f;                        // config index -> warm-start row / family
        out[f].k = 0;
        out[f].aux[0] = a->families[f];      // family value
        out[f].aux[1] = 0;
    }
    return nf;
}

static int rk_key_len(const SolverCtx *ctx, const SolverConfig *cc) {
    (void) cc;
    const RkScratch *a = (const RkScratch *) ctx->model_scratch;
    return a->known_key ? 0 : a->n;          // 0 => SWEEP (known key); N => climb (blind)
}

static void rk_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const RkScratch *a = (const RkScratch *) ctx->model_scratch;
    int n = a->n, family = cc->aux[0];
    st->aux[0] = family;                      // family
    st->aux[1] = cc->j;                       // config index
    st->key_len = n;
    if (a->known_key) {
        for (int i = 0; i < n; i++) st->key[i] = a->known_key_stream[i];
    } else {
        const int *w = &a->warm[(size_t) cc->j * n];
        double pmax = RK_JITTER_MAX * frand();          // random jitter level this restart
        for (int i = 0; i < n; i++)
            st->key[i] = (frand() < pmax) ? rand_int(0, g_alpha) : w[i];
    }
    // Anchor crib positions: the known plaintext letter forces the key letter here. Skipped
    // in known-key mode (the supplied key wins). Only anchor when the crib is ACHIEVABLE for
    // this family -- Porta maps between alphabet halves, so a crib inconsistent with the
    // ciphertext half has no key; leave such a position to the search rather than pin it wrong.
    if (a->n_cribs > 0 && !a->known_key)
        for (int p = 0; p < n; p++)
            if (a->crib_pt[p] >= 0) {
                int fk = rk_key_from_pt(ctx->cipher[p], a->crib_pt[p], family);
                if (rk_decode_char(ctx->cipher[p], fk, family) == a->crib_pt[p])
                    st->key[p] = fk;
            }
}

static void rk_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                       bool *force_primary) {
    (void) cc; (void) force_primary;
    const RkScratch *a = (const RkScratch *) ctx->model_scratch;
    int n = st->key_len;
    if (n <= 0) return;

    // With cribs, only the free (unpinned) positions are searched -- pinned positions stay
    // at their crib-forced key value throughout, so the crib region is always exact.
    if (a->n_cribs > 0) {
        if (a->n_free <= 0) return;
        int i = a->free_pos[rand_int(0, a->n_free)];
        st->key[i] = rand_int(0, g_alpha);
        return;
    }

    double r = frand();
    if (r < 0.85) {                                     // single position -> random letter
        int i = rand_int(0, n);
        st->key[i] = rand_int(0, g_alpha);
    } else if (r < 0.95) {                              // short run of 2-3 positions
        int i = rand_int(0, n);
        int run = 2 + rand_int(0, 2);
        for (int t = 0; t < run && i + t < n; t++) st->key[i + t] = rand_int(0, g_alpha);
    } else {                                            // small +/- shift of one position
        int i = rand_int(0, n);
        int d = 1 + rand_int(0, 3);
        if (frand() < 0.5) d = -d;
        st->key[i] = (st->key[i] + d + g_alpha) % g_alpha;
    }
}

static void rk_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int n = src->key_len;
    memcpy(dst->key, src->key, (size_t) n * sizeof(int));
    dst->key_len = n;
    dst->aux[0] = src->aux[0];
    dst->aux[1] = src->aux[1];
}

static void rk_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                            int *out, double *score_adjust) {
    (void) cc;
    const RkScratch *a = (const RkScratch *) ctx->model_scratch;
    ColossusConfig *cfg = ctx->cfg;
    int family = st->aux[0], n = ctx->cipher_len;

    running_key_decrypt(out, ctx->cipher, n, st->key, family);          // out = plaintext half

    // Fold in the KEY stream's own English fitness (the joint dual-stream objective): the
    // engine adds this to state_score(out), which supplies the plaintext half. Both terms
    // must sit on the SAME scale -- with no crib, state_score returns the RAW mean-n-gram
    // of the plaintext (no weight_ngram factor), so the key term is the raw mean-n-gram of
    // the key too, giving the two streams EQUAL weight (== the mean n-gram of the whole
    // K||P passage). Weighting the key by weight_ngram here would swamp the plaintext ~12:1
    // and let the search maximise key fluency alone -> fluent-but-wrong decrypts. Self-keyed
    // adds the small K|P seam reward (the only term tying the two halves into one passage).
    double adj = ngram_score(st->key, n, ctx->ngram_data, cfg->ngram_size);
    if (a->self_keyed)
        adj += rk_seam_adjust(st->key, out, n, cfg->ngram_size, ctx->ngram_data);
    // (Cribs need no score term -- they are hard-anchored in seed/perturb, so the crib
    //  region is exact by construction; only the n-gram of the free positions is optimised.)
    *score_adjust = adj;
}

// ===================================================================
//  Reporting
// ===================================================================

static _Thread_local int g_rk_full[MAX_CIPHER_LENGTH];   // assembled reported plaintext (K||P or P)

static void rk_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) ctx; (void) cc; (void) decrypted;
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  family %s, score=%.4f  [%.1fs, %d restarts]\n",
        rk_family_name(st->aux[0]), score, elapsed, stats->n_restarts);
    fflush(stdout);
}

static void rk_report(const SolverCtx *ctx, const SolverConfig *cc,
                      const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const RkScratch *a = (const RkScratch *) ctx->model_scratch;
    int n = ctx->cipher_len, family = st->aux[0];

    // Assemble the reported plaintext: the full 2N passage (K||P) when self-keyed, else
    // the plaintext half only (the key is external / independent).
    int outlen;
    if (a->self_keyed) {
        for (int i = 0; i < n; i++) g_rk_full[i]     = st->key[i];
        for (int i = 0; i < n; i++) g_rk_full[n + i] = decrypted[i];
        outlen = 2 * n;
    } else {
        for (int i = 0; i < n; i++) g_rk_full[i] = decrypted[i];
        outlen = n;
    }

    int n_words_found = 0;
    static _Thread_local char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < outlen; i++) plaintext_string[i] = index_to_char(g_rk_full[i]);
    plaintext_string[outlen] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    const char *mode = a->known_key ? "known-key" : (a->self_keyed ? "self-keyed" : "independent");
    printf("\nResult Score: %.2f | Words: %d | family=%s | mode=%s\n",
        score, n_words_found, rk_family_name(family), mode);

    print_text(ctx->cipher, n);
    printf("\n");
    printf("Running key: ");
    print_text(st->key, n);
    printf("\n");
    print_text(g_rk_full, outlen);
    printf("\n");
    print_spaces_line(g_spaces_table, g_rk_full, outlen);
    printf("%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(g_rk_full, ctx->result->decrypted, outlen);
        ctx->result->decrypted_len = outlen;
    }

    // One-liner summary: >>> score, [words,] type, family=, mode=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, family=%s, mode=%s, ",
            score, n_words_found, cfg->cipher_type, rk_family_name(family), mode);
    else
        printf(">>> %.2f, %d, family=%s, mode=%s, ",
            score, cfg->cipher_type, rk_family_name(family), mode);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(ctx->cipher, n);
    printf(", ");
    print_text(g_rk_full, outlen);
    printf("\n");
}

static const CipherModel RUNNING_KEY_MODEL = {
    .name = "running-key", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = rk_enumerate, .key_len = rk_key_len,
    .seed = rk_seed, .perturb = rk_perturb, .copy_state = rk_copy,
    .decrypt = rk_decrypt_hook, .report = rk_report,
    .report_verbose = rk_report_verbose,
};

// ===================================================================
//  Key-file loader (known-key mode)
// ===================================================================

// Load a running-key TEXT from a file, keeping only letters (uppercased), into out[]
// as alphabet indices. Returns the number of letters read (<= cap).
static int rk_load_key(const char *filename, int out[], int cap) {
    FILE *fp = fopen(filename, "r");
    if (!fp) return -1;
    int len = 0, ch;
    while (len < cap && (ch = fgetc(fp)) != EOF) {
        int idx = g_char_to_idx[toupper((unsigned char) ch) & 127];
        if (idx >= 0 && idx < g_alpha) out[len++] = idx;
    }
    fclose(fp);
    return len;
}

// ===================================================================
//  Entry point
// ===================================================================

void solve_running_key(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Running Key needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < cfg->ngram_size + 1) {
        printf("\n\nERROR: ciphertext too short for a Running Key solve.\n\n");
        return;
    }
    if (2 * cipher_len > MAX_CIPHER_LENGTH) {
        printf("\n\nERROR: Running Key ciphertext too long (2N must be <= %d).\n\n", MAX_CIPHER_LENGTH);
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Running Key ciphertext has a non-alphabet symbol at "
                   "position %d; it must be solid letters.\n\n", i);
            return;
        }

    // Families to sweep (n-gram picks the winner), or one pinned by -variant / -beaufort.
    RkScratch scratch;
    scratch.n = cipher_len;
    scratch.n_families = 0;
    if (cfg->beaufort) {
        scratch.families[scratch.n_families++] = RK_BEAUFORT;
    } else if (cfg->variant) {
        scratch.families[scratch.n_families++] = RK_VARIANT;
    } else {
        scratch.families[scratch.n_families++] = RK_VIGENERE;
        scratch.families[scratch.n_families++] = RK_BEAUFORT;
        scratch.families[scratch.n_families++] = RK_VARIANT;
        scratch.families[scratch.n_families++] = RK_PORTA;
    }

    scratch.known_key = cfg->runningkey_present;
    scratch.self_keyed = !cfg->runningkey_present && !cfg->runningkey_independent;
    scratch.warm = NULL;
    scratch.known_key_stream = NULL;
    scratch.crib_pt = NULL;
    scratch.free_pos = NULL;
    scratch.n_free = cipher_len;
    scratch.n_cribs = n_cribs;

    // Build the crib anchor map: crib_pt[p] = crib plaintext letter (or -1), and the list of
    // free (searchable) positions. Cribs align with the ciphertext = the plaintext half.
    int *crib_pt = NULL, *free_pos = NULL;
    if (n_cribs > 0) {
        crib_pt  = (int *) malloc((size_t) cipher_len * sizeof(int));
        free_pos = (int *) malloc((size_t) cipher_len * sizeof(int));
        if (!crib_pt || !free_pos) { printf("\n\nERROR: out of memory.\n\n"); free(crib_pt); free(free_pos); return; }
        for (int p = 0; p < cipher_len; p++) crib_pt[p] = -1;
        for (int i = 0; i < n_cribs; i++)
            if (crib_positions[i] >= 0 && crib_positions[i] < cipher_len)
                crib_pt[crib_positions[i]] = crib_indices[i];
        int nf = 0;
        for (int p = 0; p < cipher_len; p++) if (crib_pt[p] < 0) free_pos[nf++] = p;
        scratch.crib_pt = crib_pt;
        scratch.free_pos = free_pos;
        scratch.n_free = nf;
    }

    int *known = NULL;
    if (cfg->runningkey_present) {
        known = (int *) malloc((size_t) MAX_CIPHER_LENGTH * sizeof(int));
        if (!known) { printf("\n\nERROR: out of memory.\n\n"); free(crib_pt); free(free_pos); return; }
        int klen = rk_load_key(cfg->runningkey_file, known, MAX_CIPHER_LENGTH);
        if (klen < 0) {
            printf("\n\nERROR: could not open running-key file '%s'.\n\n", cfg->runningkey_file);
            free(known); free(crib_pt); free(free_pos); return;
        }
        if (klen < cipher_len) {
            printf("\n\nERROR: running-key text has %d letters, shorter than the %d-letter "
                   "ciphertext.\n\n", klen, cipher_len);
            free(known); free(crib_pt); free(free_pos); return;
        }
        scratch.known_key_stream = known;
        if (cfg->verbose)
            printf("\nrunning-key: %d letters, KNOWN key '%s' (%d letters), family sweep\n",
                cipher_len, cfg->runningkey_file, klen);
    } else {
        // Blind: a deterministic beam warm start per family, then anneal K.
        scratch.warm = (int *) malloc((size_t) scratch.n_families * cipher_len * sizeof(int));
        if (!scratch.warm) { printf("\n\nERROR: out of memory.\n\n"); free(crib_pt); free(free_pos); return; }
        if (cfg->verbose)
            printf("\nrunning-key: %d letters, %s blind, %d-family sweep, beam warm start\n",
                cipher_len, scratch.self_keyed ? "self-keyed" : "independent", scratch.n_families);
        for (int f = 0; f < scratch.n_families; f++)
            rk_beam_warmstart(cipher_indices, cipher_len, scratch.families[f],
                shared->ngram_data, cfg->ngram_size, RK_BEAM_DEFAULT,
                &scratch.warm[(size_t) f * cipher_len]);
    }

    // n_cribs = 0 to the engine: cribs are hard-anchored in seed/perturb (not scored), so
    // state_score returns the plaintext half's raw n-gram, keeping the two streams balanced.
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    run_solver(&RUNNING_KEY_MODEL, &ctx);

    if (scratch.warm) free(scratch.warm);
    if (known) free(known);
    if (crib_pt) free(crib_pt);
    if (free_pos) free(free_pos);
}
