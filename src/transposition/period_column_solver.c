#include "period_column_solver.h"
#include "scoring.h"
#include "trans_common.h"

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

// Format "[4x42,TP,P:3][56x3,UTP,P:2]" for the report / >>> CSV param field.
static void pc_param_summary(int len, const PCStage *stages, int n_stages,
                             char *buf, size_t n) {
    size_t off = 0;
    off += snprintf(buf + off, n - off, "depth=%d ", n_stages);
    for (int s = 0; s < n_stages && off < n; s++) {
        int dx = stages[s].dx, dy = (len + dx - 1) / dx;
        off += snprintf(buf + off, n - off, "[%dx%d,%s,P:%d]",
                        dx, dy, stages[s].utp == 0 ? "TP" : "UTP", stages[s].period);
    }
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

    static int best_decrypted[MAX_CIPHER_LENGTH];
    PCStage stages[2];
    int n_stages = 0;
    double score = period_column_search(cipher_indices, cipher_len, max_depth,
        shared->ngram_data, cfg->ngram_size,
        crib_indices, crib_positions, n_cribs,
        cfg->weight_ngram, cfg->weight_crib,
        best_decrypted, stages, &n_stages);

    char params[128];
    pc_param_summary(cipher_len, stages, n_stages, params, sizeof params);

    printf("\nperiod-column: %d composed stage(s) (widths are complete-grid divisors of %d)\n",
           n_stages, cipher_len);
    for (int s = 0; s < n_stages; s++) {
        int dx = stages[s].dx, dy = (cipher_len + dx - 1) / dx;
        printf("  stage %d: %d x %d grid, period %d, %s\n", s + 1, dx, dy,
               stages[s].period, stages[s].utp == 0 ? "transpose" : "untranspose");
    }

    report_transposition(cfg, shared, cipher_indices, cipher_len, best_decrypted,
        score, cribtext_str, n_cribs, params);
}
