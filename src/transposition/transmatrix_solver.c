#include "transmatrix_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// Transposition cipher solvers
//
// Pure transposition ciphers (transmatrix, transperoffset) are cracked by
// optimizing the small parameter vector of the transform itself, rather than
// enumerating every configuration. The search reuses the same shotgun /
// slippery hill-climbing engine and n-gram scoring as the polyalphabetic path.
//   - TRANSMATRIX    : params (p1,p2,p3) = (w1, w2, clockwise)
//   - TRANSPEROFFSET : params (p1,p2)    = (period d, offset n)

// Apply a candidate transform to a fresh copy of the ciphertext.
static void apply_transposition(ColossusConfig *cfg, int cipher_indices[],
    int cipher_len, int p1, int p2, int p3, int decrypted[]) {

    vec_copy(cipher_indices, decrypted, cipher_len);
    if (cfg->cipher_type == TRANSMATRIX) {
        transmatrix(decrypted, cipher_len, p1, p2, p3); // p1=w1, p2=w2, p3=clockwise
    } else { // TRANSPEROFFSET
        transperoffset(decrypted, cipher_len, p1, p2);  // p1=period(d), p2=offset(n)
    }
}

// Pick a random valid parameter triple for the current transposition type.
static void random_transposition_params(ColossusConfig *cfg, int cipher_len,
    int *p1, int *p2, int *p3) {

    if (cfg->cipher_type == TRANSMATRIX) {
        // Valid grid widths are [2, len-1]; matrix_rotate is identity outside.
        *p1 = rand_int(2, cipher_len);
        *p2 = rand_int(2, cipher_len);
        *p3 = rand_int(0, 2);            // clockwise in {0,1}
    } else { // TRANSPEROFFSET
        // The decimation step d must be coprime to len for a bijection.
        do {
            *p1 = rand_int(1, cipher_len);
        } while (gcd(*p1, cipher_len) != 1);
        *p2 = rand_int(0, cipher_len);   // offset n in [0, len-1]
        *p3 = 0;                         // unused
    }
}

// Perturb a single parameter (one local neighbour move).
static void perturbate_transposition_params(ColossusConfig *cfg, int cipher_len,
    int *p1, int *p2, int *p3) {

    if (cfg->cipher_type == TRANSMATRIX) {
        int which = rand_int(0, 3);
        if (which == 2) {
            *p3 = 1 - *p3;               // flip direction
        } else {
            int *w = (which == 0) ? p1 : p2;
            if (cipher_len > 4 && frand() < 0.5) {
                // Local +/-1 step, clamped to [2, len-1].
                *w += (frand() < 0.5) ? -1 : 1;
                if (*w < 2) *w = 2;
                if (*w > cipher_len - 1) *w = cipher_len - 1;
            } else {
                *w = rand_int(2, cipher_len);
            }
        }
    } else { // TRANSPEROFFSET
        if (frand() < 0.5) {
            // Local rotation step: aligns the wrap boundary / any crib.
            *p2 += (frand() < 0.5) ? -1 : 1;
            *p2 = ((*p2 % cipher_len) + cipher_len) % cipher_len;
        } else {
            // d is non-smooth, so re-randomise among coprime values.
            do {
                *p1 = rand_int(1, cipher_len);
            } while (gcd(*p1, cipher_len) != 1);
        }
    }
}

// ---- transmatrix / transperoffset model (cipher-agnostic engine) ----------
// A pure parameter-vector search: the candidate state is the transform's own
// (p1,p2,p3) triple carried in st->aux[0..2]. SHAPE_SHOTGUN reproduces the
// slip-probability acceptance and RNG draw order of the original dedicated
// climber, so the search is bit-identical at a fixed seed.

static int transmat_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    (void) ctx; (void) cap;
    out[0].period = 0; out[0].j = 0; out[0].k = 0;
    out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;     // a single config; the whole search is the inner climb
}

static void transmat_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    random_transposition_params(ctx->cfg, ctx->cipher_len, &st->aux[0], &st->aux[1], &st->aux[2]);
}

static void transmat_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                             bool *force_primary) {
    (void) cc; (void) force_primary;
    perturbate_transposition_params(ctx->cfg, ctx->cipher_len, &st->aux[0], &st->aux[1], &st->aux[2]);
}

static void transmat_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    dst->aux[0] = src->aux[0]; dst->aux[1] = src->aux[1]; dst->aux[2] = src->aux[2];
}

static void transmat_decrypt(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                             int *out, double *score_adjust) {
    (void) cc; (void) score_adjust;
    apply_transposition(ctx->cfg, ctx->cipher, ctx->cipher_len,
        st->aux[0], st->aux[1], st->aux[2], out);
}

static void transmat_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    int buf[MAX_CIPHER_LENGTH];
    apply_transposition(cfg, ctx->cipher, ctx->cipher_len, st->aux[0], st->aux[1], st->aux[2], buf);

    double elapsed = engine_elapsed_sec(stats);
    double n_iter_per_sec = (elapsed > 0.) ? ((double) stats->n_iterations)/elapsed : 0.;

    printf("\n%.2f\t[sec]\n", elapsed);
    printf("%.0fK\t[it/sec/thread]\n", 1.e-3*n_iter_per_sec);
    printf("%d\t[restarts]\n", stats->n_restarts);
    printf("%d\t[backtracks]\n", stats->n_backtracks);
    printf("%d\t[slips]\n", stats->n_slips);
    printf("%.4f\t[entropy]\n", entropy(buf, ctx->cipher_len));
    printf("%.2f\t[score]\n", score);
    if (cfg->cipher_type == TRANSMATRIX) {
        printf("w1 = %d, w2 = %d, direction = %s\t[params]\n",
            st->aux[0], st->aux[1], st->aux[2] ? "cw" : "ccw");
    } else { // TRANSPEROFFSET
        printf("period = %d, offset = %d\t[params]\n", st->aux[0], st->aux[1]);
    }
    printf("\n");
    print_text(buf, ctx->cipher_len); printf("\n");
    fflush(stdout);
}

static void transmat_report(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    SharedData *shared = ctx->shared;
    int cipher_len = ctx->cipher_len;
    int *cipher_indices = ctx->cipher;
    char *cribtext_str = ctx->cribtext;
    int n_cribs = ctx->n_cribs;
    int best_p1 = st->aux[0], best_p2 = st->aux[1], best_p3 = st->aux[2];
    int n_words_found = 0;

    // Parameter report line.
    if (cfg->cipher_type == TRANSMATRIX) {
        printf("\ntransmatrix: w1 = %d, w2 = %d, direction = %s\n",
            best_p1, best_p2, best_p3 ? "cw" : "ccw");
    } else {
        printf("\ntransperiodoffset: period = %d, offset = %d\n", best_p1, best_p2);
    }

    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < cipher_len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[cipher_len] = '\0';

    if (cfg->dictionary_present && shared->dict != NULL) {
        n_words_found = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);
    }

    // Results output.
    printf("\nResult Score: %.2f | Words: %d\n", score, n_words_found);

    print_text(cipher_indices, cipher_len);
    printf("\n");
    print_text(decrypted, cipher_len);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, cipher_len);
    printf("%s\n", cribtext_str);

    if (PARTIAL_CRIB_MATCH && n_cribs > 0) {
        // Indexed by cipher position via cribtext_str (same scheme as solve_cipher).
        for (int i = 0; i < cipher_len; i++) {
            if (cribtext_str[i] == '_') {
                printf("_");
            } else {
                int diff = abs(decrypted[i] - (g_char_to_idx[toupper((unsigned char)cribtext_str[i]) & 127]));
                if (diff < 10) printf("%d", diff); else printf("*");
            }
        }
    }
    printf("\n\n");

    // One-liner summary. Field order matches the existing transmatrix /
    // transperoffset summaries, minus the keyword/cycleword fields (a pure
    // transposition has none).
    if (cfg->cipher_type == TRANSMATRIX) {
        if (cfg->dictionary_present) {
            printf(">>> %.2f, %d, %d, %d, %d, %d, %s, ", score, n_words_found,
                cfg->cipher_type, best_p1, best_p2, best_p3,
                cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        } else {
            printf(">>> %.2f, %d, %d, %d, %d, %s, ", score,
                cfg->cipher_type, best_p1, best_p2, best_p3,
                cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        }
    } else {
        if (cfg->dictionary_present) {
            printf(">>> %.2f, %d, %d, %d, %d, %s, ", score, n_words_found,
                cfg->cipher_type, best_p1, best_p2,
                cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        } else {
            printf(">>> %.2f, %d, %d, %d, %s, ", score,
                cfg->cipher_type, best_p1, best_p2,
                cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        }
    }

    print_text(cipher_indices, cipher_len);
    printf(", ");
    print_text(decrypted, cipher_len);
    printf("\n");
}

static const CipherModel TRANSMATRIX_MODEL = {
    .name = "transmatrix",
    .shape = SHAPE_SHOTGUN,
    .needs_hist = false,
    .enumerate_configs = transmat_enumerate,
    .key_len = NULL,
    .seed = transmat_seed,
    .perturb = transmat_perturb,
    .copy_state = transmat_copy,
    .decrypt = transmat_decrypt,
    .report = transmat_report,
    .report_verbose = transmat_report_verbose,
};

void solve_transposition(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str; // ciphertext is carried as cipher_indices.

    if (cipher_len < 3) {
        printf("\n\nERROR: ciphertext too short for a transposition solve.\n\n");
        return ;
    }

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    run_solver(&TRANSMATRIX_MODEL, &ctx);
}


