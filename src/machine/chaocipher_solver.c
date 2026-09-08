#include "chaocipher_solver.h"
#include "chaocipher.h"
#include "engine.h"
#include "scoring.h"
#include <string.h>

// =====================================================================
//  Chaocipher solver (TYPE chaocipher / chao)
// =====================================================================
//
// Chaocipher's key is the pair of STARTING alphabets (left = ciphertext disk, right =
// plaintext disk); everything downstream is deterministic (see chaocipher.c).
//
// KEY FINDING (characterised in tests/test_chaocipher_solver.c): Chaocipher is a NEEDLE for
// local search. A single starting-alphabet swap cascades the whole downstream decrypt, so
// the n-gram (and displacement-error, and forward-match) landscapes are flat/rugged
// everywhere except exactly at the solution -- blind AND full-known-plaintext annealing both
// stall at ~random even at 460 letters with a large budget. This is the same isolated-needle
// situation as Condi, and the same conclusion follows: don't hill-climb a needle, reconstruct
// it. So the two regimes are:
//
//  * BLIND ciphertext-only (default): a SHAPE_ANNEAL n-gram climb over the two starting
//    alphabets. Shipped for completeness (and so -method/-nthreads/PSO all apply), but it is
//    a DOCUMENTED LIMITATION -- it does not reliably recover the key at any practical length.
//
//  * KNOWN-PLAINTEXT (a contiguous known-plaintext prefix supplied via -crib): a DETERMINISTIC
//    BACKTRACKING RECONSTRUCTION -- the method used to solve the real Chaocipher exhibits, and
//    the natural fit for this repo's "deterministic/constructive for needles" rule. It walks
//    the known plaintext + ciphertext streams position by position, maintaining the two
//    evolving disks as partially-known permutations plus an ORIGIN-INDEX map back to the
//    STARTING frame; each step pins the shared position of pt[s]/ct[s], places letters into the
//    starting alphabets, and prunes on contradictions. The rotation gauge (rotating both
//    alphabets by the same amount is an equivalent key) fixes the first placement. It recovers
//    the key (up to that gauge) from a known-plaintext prefix and then decrypts the whole
//    message, including any un-cribbed tail. Typical solves are well under 0.1s; a few keys
//    take up to ~10s of early branching (bounded by a node cap).

// --- blind anneal reward weight (none: n-gram only) ---------------------------------------

// Constructive-KPA tuning.
#define CHAO_KPA_MIN_KNOWN    28            // minimum leading known-pt letters to attempt reconstruction
#define CHAO_KPA_PROCESS_CAP  2000          // recursion-depth safety bound (real cribs stop far sooner via
                                            //   the both-disks-fully-placed early exit; this only caps the
                                            //   rare case where some letter never appears, keeping the stack safe)
#define CHAO_KPA_MAX_NODES    250000000L    // backtracking node cap (safety valve; ~10s worst case)

// =====================================================================
//  Shared report helper (used by both the constructive path and the anneal report hook)
// =====================================================================

// Canonical representative: rotate BOTH alphabets so the LEFT disk begins with 'A' (index 0).
// Rotating both by the same offset is an equivalent key, so the decrypt is unchanged.
static void chao_canonicalize(const int left[CHAO_N], const int right[CHAO_N],
                              int cl[CHAO_N], int cr[CHAO_N]) {
    int k = 0; while (k < CHAO_N && left[k] != 0) k++;
    if (k >= CHAO_N) k = 0;
    for (int i = 0; i < CHAO_N; i++) { cl[i] = left[(i + k) % CHAO_N]; cr[i] = right[(i + k) % CHAO_N]; }
}

static void chao_alpha_string(const int a[CHAO_N], char out[]) {
    for (int i = 0; i < CHAO_N; i++) out[i] = index_to_char(a[i]);
    out[CHAO_N] = '\0';
}

static void chao_emit_report(ColossusConfig *cfg, SharedData *shared, int *cipher, int n,
        const int left[CHAO_N], const int right[CHAO_N], int *decrypted, char *cribtext,
        double score, SolveResult *result) {

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);

    int cl[CHAO_N], cr[CHAO_N];
    chao_canonicalize(left, right, cl, cr);
    char ls[CHAO_N + 1], rs[CHAO_N + 1];
    chao_alpha_string(cl, ls); chao_alpha_string(cr, rs);

    printf("\nResult Score: %.2f | Words: %d | left(ct)=%s | right(pt)=%s\n",
        score, n_words_found, ls, rs);
    print_text(cipher, n);
    printf("\n");
    print_text(decrypted, n);
    printf("\n");
    print_solution_check(decrypted, n);
    print_spaces_line(g_spaces_table, decrypted, n);
    if (cribtext) printf("%s\n", cribtext);

    if (result) {
        result->solved = true;
        result->cipher_type = cfg->cipher_type;
        result->score = score;
        result->n_words = n_words_found;
        result->cycleword_len = CHAO_N;
        for (int i = 0; i < CHAO_N; i++) {
            result->ciphertext_keyword[i] = cl[i];   // recovered LEFT (ciphertext) disk
            result->plaintext_keyword[i]  = cr[i];    // recovered RIGHT (plaintext) disk
        }
        vec_copy(decrypted, result->decrypted, n);
        result->decrypted_len = n;
    }

    // One-liner summary: >>> score, [words,] type, left=, right=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, left=%s, right=%s, ",
            score, n_words_found, cfg->cipher_type, ls, rs);
    else
        printf(">>> %.2f, %d, left=%s, right=%s, ",
            score, cfg->cipher_type, ls, rs);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(cipher, n);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

// =====================================================================
//  Constructive known-plaintext reconstruction (deterministic backtracking)
// =====================================================================

static const int *g_kpa_ct;                 // ciphertext stream (0..25)
static const int *g_kpa_pt;                 // known plaintext stream (0..25)
static int  g_kpa_len;                      // number of known (pt,ct) steps to process
static long g_kpa_nodes;
static int  g_kpa_abort;
static int  g_kpa_startL[CHAO_N], g_kpa_startR[CHAO_N];   // recovered STARTING alphabets (-1 = unknown)

static int chao_idxof(const int a[CHAO_N], int v) {
    for (int i = 0; i < CHAO_N; i++) if (a[i] == v) return i;
    return -1;
}

// Verify the remaining steps s..len-1 against fully-placed disks (no branching left).
static int chao_kpa_verify(const int curL[CHAO_N], const int curR[CHAO_N], int s) {
    int L[CHAO_N], R[CHAO_N];
    for (int k = 0; k < CHAO_N; k++) { L[k] = curL[k]; R[k] = curR[k]; }
    for (; s < g_kpa_len; s++) {
        int p = g_kpa_pt[s];
        int i = chao_idxof(R, p);
        if (i < 0 || L[i] != g_kpa_ct[s]) return 0;
        chaocipher_step_left(L, i);
        chaocipher_step_right(R, i);
    }
    return 1;
}

// DFS over the known stream. curL/curR are the evolving disks (-1 = unknown cell); oL/oR map
// each current cell back to its STARTING-frame index (so placing a letter records it into
// g_kpa_start{L,R}). Returns 1 on success (leaving g_kpa_start{L,R} holding the solution).
static int chao_kpa_dfs(int s, int curL[CHAO_N], int curR[CHAO_N], int oL[CHAO_N], int oR[CHAO_N],
                        int placed_first, int nL, int nR) {
    if (g_kpa_abort) return 0;
    if (++g_kpa_nodes > CHAO_KPA_MAX_NODES) { g_kpa_abort = 1; return 0; }

    if (nL == CHAO_N && nR == CHAO_N)                 // both disks fully placed: verify the rest
        return chao_kpa_verify(curL, curR, s) ? 1 : 0;
    if (s >= g_kpa_len || s >= CHAO_KPA_PROCESS_CAP)  // consumed the known stream (partial key ok)
        return 1;

    int p = g_kpa_pt[s], c = g_kpa_ct[s];
    int ip = chao_idxof(curR, p), ic = chao_idxof(curL, c);

    int cand[CHAO_N], nc = 0;
    if (ip >= 0 && ic >= 0) {                         // both letters already placed: positions must agree
        if (ip != ic) return 0;
        cand[nc++] = ip;
    } else if (ip >= 0) {                             // pt placed at ip -> ct must go there
        if (curL[ip] != -1 && curL[ip] != c) return 0;
        cand[nc++] = ip;
    } else if (ic >= 0) {                             // ct placed at ic -> pt must go there
        if (curR[ic] != -1 && curR[ic] != p) return 0;
        cand[nc++] = ic;
    } else {                                          // both new: place at a both-unknown position
        if (!placed_first) cand[nc++] = 0;            // rotation gauge fixes the first placement
        else for (int j = 0; j < CHAO_N; j++) if (curL[j] == -1 && curR[j] == -1) cand[nc++] = j;
    }

    for (int t = 0; t < nc; t++) {
        int i = cand[t];
        int cL[CHAO_N], cR[CHAO_N], dL[CHAO_N], dR[CHAO_N], sL[CHAO_N], sR[CHAO_N];
        memcpy(cL, curL, sizeof cL); memcpy(cR, curR, sizeof cR);
        memcpy(dL, oL, sizeof dL);   memcpy(dR, oR, sizeof dR);
        memcpy(sL, g_kpa_startL, sizeof sL); memcpy(sR, g_kpa_startR, sizeof sR);
        int npf = placed_first, nnL = nL, nnR = nR, ok = 1;

        if (curR[i] == -1) {
            if (g_kpa_startR[oR[i]] != -1) ok = 0;
            else { curR[i] = p; g_kpa_startR[oR[i]] = p; nnR++; npf = 1; }
        } else if (curR[i] != p) ok = 0;

        if (ok) {
            if (curL[i] == -1) {
                if (g_kpa_startL[oL[i]] != -1) ok = 0;
                else { curL[i] = c; g_kpa_startL[oL[i]] = c; nnL++; }
            } else if (curL[i] != c) ok = 0;
        }

        if (ok) {
            chaocipher_step_left(curL, i);  chaocipher_step_left(oL, i);
            chaocipher_step_right(curR, i); chaocipher_step_right(oR, i);
            if (chao_kpa_dfs(s + 1, curL, curR, oL, oR, npf, nnL, nnR)) return 1;
        }
        memcpy(curL, cL, sizeof cL); memcpy(curR, cR, sizeof cR);
        memcpy(oL, dL, sizeof dL);   memcpy(oR, dR, sizeof dR);
        memcpy(g_kpa_startL, sL, sizeof sL); memcpy(g_kpa_startR, sR, sizeof sR);
    }
    return 0;
}

// Fill any unknown (-1) cells of a partial alphabet with the leftover letters (ascending), so
// it is a full permutation. Cells never exercised by the crib are unrecoverable; filling them
// deterministically leaves every pinned (recovered) cell -- and thus the decrypt -- unchanged.
static void chao_complete_alphabet(int a[CHAO_N]) {
    int used[CHAO_N] = {0};
    for (int i = 0; i < CHAO_N; i++) if (a[i] >= 0) used[a[i]] = 1;
    int nxt = 0;
    for (int i = 0; i < CHAO_N; i++)
        if (a[i] < 0) { while (nxt < CHAO_N && used[nxt]) nxt++; a[i] = nxt; used[nxt] = 1; }
}

// Recover the two starting alphabets from the first K (pt,ct) pairs. Returns true and fills
// startL/startR (completed to full permutations) on success.
static bool chao_kpa_reconstruct(const int *cipher, const int *known_pt, int K,
                                 int startL[CHAO_N], int startR[CHAO_N]) {
    g_kpa_ct = cipher; g_kpa_pt = known_pt; g_kpa_len = K;
    g_kpa_nodes = 0; g_kpa_abort = 0;
    for (int i = 0; i < CHAO_N; i++) { g_kpa_startL[i] = -1; g_kpa_startR[i] = -1; }

    int curL[CHAO_N], curR[CHAO_N], oL[CHAO_N], oR[CHAO_N];
    for (int i = 0; i < CHAO_N; i++) { curL[i] = -1; curR[i] = -1; oL[i] = i; oR[i] = i; }

    if (!chao_kpa_dfs(0, curL, curR, oL, oR, 0, 0, 0)) return false;

    for (int i = 0; i < CHAO_N; i++) { startL[i] = g_kpa_startL[i]; startR[i] = g_kpa_startR[i]; }
    chao_complete_alphabet(startL);
    chao_complete_alphabet(startR);
    return true;
}

// =====================================================================
//  Blind n-gram anneal (CipherModel) -- documented limitation
// =====================================================================

void chaocipher_perturb_state(SolverState *st) {
    // Swap two positions in one alphabet (chosen at random); occasionally also swap in the
    // other, so the two disks mix rather than stall with one alphabet nearly right.
    int *a = (frand() < 0.5) ? st->ct_keyword : st->pt_keyword;
    int i = rand_int(0, CHAO_N), j = rand_int(0, CHAO_N);
    while (j == i) j = rand_int(0, CHAO_N);
    int t = a[i]; a[i] = a[j]; a[j] = t;

    if (frand() < 0.15) {
        int *b = (a == st->ct_keyword) ? st->pt_keyword : st->ct_keyword;
        int p = rand_int(0, CHAO_N), q = rand_int(0, CHAO_N);
        while (q == p) q = rand_int(0, CHAO_N);
        t = b[p]; b[p] = b[q]; b[q] = t;
    }
}

static int chao_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    (void) ctx;
    if (cap < 1) return 0;
    out[0].period = 0; out[0].j = 0; out[0].k = 0;
    out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;                                    // single config: no period / sweep
}

static void chao_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx; (void) cc;
    for (int i = 0; i < CHAO_N; i++) { st->ct_keyword[i] = i; st->pt_keyword[i] = i; }
    shuffle(st->ct_keyword, CHAO_N);
    shuffle(st->pt_keyword, CHAO_N);
    st->key_len = CHAO_N;
}

static void chao_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                         SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    chaocipher_perturb_state(st);
}

static void chao_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < CHAO_N; i++) {
        dst->ct_keyword[i] = src->ct_keyword[i];
        dst->pt_keyword[i] = src->pt_keyword[i];
    }
    dst->key_len = src->key_len;
}

static void chao_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                              SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    chaocipher_decrypt(ctx->cipher, ctx->cipher_len, st->ct_keyword, st->pt_keyword, out);
    *score_adjust = 0.0;                         // no decoupling reward (a needle; see header)
}

static void chao_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) ctx; (void) cc; (void) decrypted;
    int cl[CHAO_N], cr[CHAO_N];
    chao_canonicalize(st->ct_keyword, st->pt_keyword, cl, cr);
    char ls[CHAO_N + 1], rs[CHAO_N + 1];
    chao_alpha_string(cl, ls); chao_alpha_string(cr, rs);
    printf("\n  score=%.4f  [%.1fs, %d restarts]\n    left(ct)=%s\n    right(pt)=%s\n",
        score, engine_elapsed_sec(stats), stats->n_restarts, ls, rs);
    fflush(stdout);
}

static void chao_report(const SolverCtx *ctx, const SolverConfig *cc,
                        const SolverState *st, double score, int *decrypted) {
    (void) cc;
    chao_emit_report(ctx->cfg, ctx->shared, ctx->cipher, ctx->cipher_len,
        st->ct_keyword, st->pt_keyword, decrypted, ctx->cribtext, score, ctx->result);
}

static const CipherModel CHAOCIPHER_MODEL = {
    .name = "chaocipher", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = chao_enumerate, .key_len = NULL,
    .seed = chao_seed, .perturb = chao_perturb, .copy_state = chao_copy,
    .decrypt = chao_decrypt_hook, .report = chao_report,
    .report_verbose = chao_report_verbose,
};

// =====================================================================
//  Entry point
// =====================================================================

// Per-solve known-plaintext prefix (built from cribs; single-threaded setup).
static int g_chao_known_pt[MAX_CIPHER_LENGTH];

void solve_chaocipher(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (g_alpha != CHAO_N) {
        printf("\n\nERROR: Chaocipher needs the 26-letter A..Z alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 2) {
        printf("\n\nERROR: ciphertext too short for a Chaocipher solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Chaocipher ciphertext must be solid A..Z letters "
                   "(bad symbol at position %d).\n\n", i);
            return;
        }

    // Build the known-plaintext position map from the cribs, then measure the leading
    // contiguous known-plaintext run (the constructive reconstruction needs a prefix).
    for (int i = 0; i < cipher_len; i++) g_chao_known_pt[i] = -1;
    for (int i = 0; i < n_cribs; i++) {
        int pos = crib_positions[i];
        if (pos >= 0 && pos < cipher_len && crib_indices[i] >= 0 && crib_indices[i] < CHAO_N)
            g_chao_known_pt[pos] = crib_indices[i];
    }
    int K = 0;
    while (K < cipher_len && g_chao_known_pt[K] >= 0) K++;

    // --- Known-plaintext mode: constructive reconstruction from the crib prefix ---
    if (K >= CHAO_KPA_MIN_KNOWN) {
        if (cfg->verbose)
            printf("\nchaocipher: %d letters, known-plaintext reconstruction from a %d-letter "
                   "crib prefix\n", cipher_len, K);
        int sL[CHAO_N], sR[CHAO_N];
        if (chao_kpa_reconstruct(cipher_indices, g_chao_known_pt, K, sL, sR)) {
            static int decrypted[MAX_CIPHER_LENGTH];
            chaocipher_decrypt(cipher_indices, cipher_len, sL, sR, decrypted);
            double score = ngram_score(decrypted, cipher_len,
                (float *) shared->ngram_data, cfg->ngram_size);
            chao_emit_report(cfg, shared, cipher_indices, cipher_len, sL, sR,
                decrypted, cribtext_str, score, result);
            return;
        }
        printf("\nchaocipher: constructive reconstruction hit the node cap / a contradiction; "
               "falling back to the (weak) blind anneal.\n");
    }

    // --- Blind (or partial-crib) mode: n-gram anneal (documented limitation) ---
    if (cfg->verbose)
        printf("\nchaocipher: %d letters, blind n-gram anneal over the two starting alphabets "
               "(NOTE: Chaocipher is a needle for local search -- see the solver header)\n",
               cipher_len);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.result = result;
    run_solver(&CHAOCIPHER_MODEL, &ctx);
}
