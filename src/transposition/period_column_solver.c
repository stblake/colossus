#include "period_column_solver.h"
#include "scoring.h"
#include "trans_common.h"
#include <time.h>

// =====================================================================
//  Period column order transposition solver -- see period_column_solver.h
// =====================================================================

// Upper bound on the enumerated single-stage list. Bounds the depth-2 cost at
// O(K^2) decrypt+score; K is the number of (width, period, utp) triples. For a
// typical few-hundred-char cipher K is a few hundred to a couple thousand, so
// K^2 is a few million O(len) evaluations -- sub-second. A pathologically
// composite length is clamped here (documented; the dropped tail is the
// largest-width stages, which are the least likely genuine grid shapes).
#define PC_MAX_STAGE_LIST 20000

// Enumerate every candidate single stage: each complete-grid width dx | len in
// [3, len], each period p in [2, dx-1], each direction utp in {0, 1}.
static int pc_build_stage_list(int len, PCStage *list, int cap) {
    int n = 0;
    for (int dx = 3; dx <= len; dx++) {
        if (len % dx != 0) continue;                 // complete grid only
        for (int p = 2; p < dx; p++) {
            if (n + 2 > cap) return n;
            list[n].dx = dx; list[n].period = p; list[n].utp = 0; n++;
            list[n].dx = dx; list[n].period = p; list[n].utp = 1; n++;
        }
    }
    return n;
}

double period_column_search(const int *cipher, int len, int max_depth,
    float *ngram_data, int ngram_size,
    int crib_indices[], int crib_positions[], int n_cribs,
    float weight_ngram, float weight_crib,
    int best_decrypted[], PCStage stages[], int *n_stages) {

    static PCStage list[PC_MAX_STAGE_LIST];
    static int buf1[MAX_CIPHER_LENGTH], buf2[MAX_CIPHER_LENGTH];

    int K = pc_build_stage_list(len, list, PC_MAX_STAGE_LIST);
    double best = 0.0;
    int have = 0;
    *n_stages = 0;
    if (K <= 0) {                                    // len too small / prime-ish: no stage
        for (int i = 0; i < len; i++) best_decrypted[i] = cipher[i];
        return 0.0;
    }

    // Depth 1: every single stage applied to the ciphertext.
    for (int a = 0; a < K; a++) {
        period_column_transform(cipher, buf1, len, list[a].dx, list[a].period, list[a].utp);
        double s = state_score(buf1, len, crib_indices, crib_positions, n_cribs,
                               ngram_data, ngram_size, weight_ngram, weight_crib, 0.f, 0.f);
        if (!have || s > best) {
            best = s; have = 1;
            for (int i = 0; i < len; i++) best_decrypted[i] = buf1[i];
            stages[0] = list[a]; *n_stages = 1;
        }
    }

    // Depth 2: every ordered pair. buf1 = stage-a of the cipher (hoisted out of the
    // inner loop), buf2 = stage-b of buf1. The global best across depth 1 and 2 wins,
    // so a genuinely single-stage cipher still reports one stage.
    if (max_depth >= 2) {
        for (int a = 0; a < K; a++) {
            period_column_transform(cipher, buf1, len, list[a].dx, list[a].period, list[a].utp);
            for (int b = 0; b < K; b++) {
                period_column_transform(buf1, buf2, len, list[b].dx, list[b].period, list[b].utp);
                double s = state_score(buf2, len, crib_indices, crib_positions, n_cribs,
                                       ngram_data, ngram_size, weight_ngram, weight_crib, 0.f, 0.f);
                if (s > best) {
                    best = s;
                    for (int i = 0; i < len; i++) best_decrypted[i] = buf2[i];
                    stages[0] = list[a]; stages[1] = list[b]; *n_stages = 2;
                }
            }
        }
    }
    return best;
}

// Format "[4x42,TP,P:3][56x3,UTP,P:2]" (with an optional "dir=bt " prefix when the
// cipher was read in reverse) for the report / >>> CSV param field.
static void pc_param_summary(int len, const PCStage *stages, int n_stages,
                             int read_dir, char *buf, size_t n) {
    size_t off = 0;
    if (read_dir == COL_READ_BT) off += snprintf(buf + off, n - off, "dir=bt ");
    off += snprintf(buf + off, n - off, "depth=%d ", n_stages);
    for (int s = 0; s < n_stages && off < n; s++) {
        int dx = stages[s].dx, dy = (len + dx - 1) / dx;
        off += snprintf(buf + off, n - off, "[%dx%d,%s,P:%d]",
                        dx, dy, stages[s].utp == 0 ? "TP" : "UTP", stages[s].period);
    }
}

// Reverse cipher[0..len-1] into out (read the stream bottom-to-top / backwards).
static void pc_reverse(const int *cipher, int len, int *out) {
    for (int i = 0; i < len; i++) out[i] = cipher[len - 1 - i];
}

void solve_period_column(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str;
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a period-column solve.\n\n");
        return;
    }

    int max_depth = cfg->trans_depth;
    if (max_depth < 1) max_depth = 2;
    if (max_depth > 2) {
        printf("\nNote: period-column exhaustive search supports depth <= 2; clamping -depth %d to 2.\n",
               max_depth);
        max_depth = 2;
    }

    // Read direction: forwards (tb, canonical) and/or the reversed stream (bt).
    // Default is COL_READ_TB, so a normal run is one forward search -- bit-identical
    // to before. -readdir bt reads the cipher backwards; -readdir both tries each.
    int dirs[2], ndirs = 0;
    if (cfg->read_direction == COL_READ_BOTH) { dirs[ndirs++] = COL_READ_TB; dirs[ndirs++] = COL_READ_BT; }
    else dirs[ndirs++] = cfg->read_direction;

    // -verbose: describe the (deterministic, exhaustive) search space up front so
    // the run isn't a silent block until the final report.
    if (cfg->verbose) {
        int nwidths = 0, nstage = 0;
        printf("\nperiod-column: exhaustive depth-%d search; complete-grid widths dx | %d in [3,%d]:",
               max_depth, cipher_len, cipher_len);
        for (int dx = 3; dx <= cipher_len; dx++)
            if (cipher_len % dx == 0) { printf(" %d", dx); nwidths++; nstage += 2 * (dx - 2); }
        printf("\n  %d width(s), %d candidate stage(s) => depth-1 = %d, depth-2 = %d decrypt+score(s) per direction\n",
               nwidths, nstage, nstage, max_depth >= 2 ? nstage * nstage : 0);
        printf("  read direction(s):");
        for (int d = 0; d < ndirs; d++) printf(" %s", dirs[d] == COL_READ_BT ? "bt" : "tb");
        printf("\n");
    }

    static int best_decrypted[MAX_CIPHER_LENGTH], oriented[MAX_CIPHER_LENGTH];
    static int best_cipher[MAX_CIPHER_LENGTH], dec[MAX_CIPHER_LENGTH];
    PCStage stages[2], best_stages[2];
    int n_stages = 0, best_ns = 0, best_dir = COL_READ_TB;
    double score = -1e30;
    clock_t t0 = clock();

    for (int d = 0; d < ndirs; d++) {
        if (dirs[d] == COL_READ_BT) pc_reverse(cipher_indices, cipher_len, oriented);
        else for (int i = 0; i < cipher_len; i++) oriented[i] = cipher_indices[i];

        int ns = 0;
        double s = period_column_search(oriented, cipher_len, max_depth,
            shared->ngram_data, cfg->ngram_size,
            crib_indices, crib_positions, n_cribs,
            cfg->weight_ngram, cfg->weight_crib,
            dec, stages, &ns);

        if (cfg->verbose) {
            // Engine-style running block: timing, score, params, then the current
            // best plaintext for this read direction.
            char p[128];
            pc_param_summary(cipher_len, stages, ns, dirs[d], p, sizeof p);
            double elapsed = ((double) clock() - t0) / CLOCKS_PER_SEC;
            printf("\n%.2f\t[sec]\n", elapsed);
            printf("%.4f\t[entropy]\n", entropy(dec, cipher_len));
            printf("%.2f\t[score]\n", s);
            printf("%s\t[params]\n\n", p);
            print_text(dec, cipher_len); printf("\n");
            print_solution_check(dec, cipher_len);
            fflush(stdout);
        }

        if (d == 0 || s > score) {
            score = s; best_ns = ns; best_dir = dirs[d];
            for (int i = 0; i < cipher_len; i++) { best_decrypted[i] = dec[i]; best_cipher[i] = oriented[i]; }
            best_stages[0] = stages[0]; best_stages[1] = stages[1];
        }
    }
    n_stages = best_ns;

    char params[128];
    pc_param_summary(cipher_len, best_stages, n_stages, best_dir, params, sizeof params);

    printf("\nperiod-column: %d composed stage(s), read %s (widths are complete-grid divisors of %d)\n",
           n_stages, best_dir == COL_READ_BT ? "bottom-to-top" : "top-to-bottom", cipher_len);
    for (int s = 0; s < n_stages; s++) {
        int dx = best_stages[s].dx, dy = (cipher_len + dx - 1) / dx;
        printf("  stage %d: %d x %d grid, period %d, %s\n", s + 1, dx, dy,
               best_stages[s].period, best_stages[s].utp == 0 ? "transpose" : "untranspose");
    }

    report_transposition(cfg, shared, best_cipher, cipher_len, best_decrypted,
        score, cribtext_str, n_cribs, params);
}
