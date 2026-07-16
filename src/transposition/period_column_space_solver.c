#include "period_column_space_solver.h"
#include "trans_common.h"           // report_transposition, trans_word_set, WordSet
#include "scoring.h"
#include <math.h>
#include <string.h>
#include <time.h>

// =====================================================================
//  Period column order (space-robust) solver -- see the header.
// =====================================================================
//
//  Wraps the deterministic period-column search in a FROM-SCRATCH search over
//  small EDITS to the observed cipher, keeping its spaces (which are structurally
//  significant) intact. Two symmetric edit kinds repair a poorly enciphered text:
//    * INSERT a blank "gap" cell -- a synthesised space that rides the
//      transposition like a real space, transparent to the letters-only n-gram
//      fitness but re-aligning every following letter (repairs a DROPPED char);
//    * DELETE an observed cell -- drop an extraneous cell from the stream
//      (repairs a spuriously ADDED char).
//  Either edit changes the grid layout (and hence the decrypt) and lets the
//  edited length hit a better complete-grid factorisation.

// A gap is encoded as a space sentinel, identical in kind to a real space cell
// from ord()/decode_cipher (space ' ' -> -(32)-1), so it prints back as a space
// and is dropped by the letters-only n-gram projection.
#define PC_GAP_SENT (-(int)' ' - 1)   // = -33

// Insertion-sort a short slot array ascending (count is tiny: <= max_ins).
static void gp_sort(int *g, int n) {
    for (int i = 1; i < n; i++) {
        int v = g[i], j = i - 1;
        while (j >= 0 && g[j] > v) { g[j + 1] = g[j]; j--; }
        g[j + 1] = v;
    }
}

// Materialise the edited cipher: delete the observed cells flagged in deleted[],
// and insert a gap cell at each slot in ins[] (0..N, "before observed cell i"; N
// appends). ins must be sorted ascending; repeated slots stack consecutive gaps.
// Returns the edited length M = N + n_ins - (number of deleted cells).
static int build_edited(const int *cipher, int N, const int *ins, int n_ins,
                        const unsigned char *deleted, int *out) {
    int ii = 0, a = 0;
    for (int i = 0; i <= N; i++) {
        while (ii < n_ins && ins[ii] == i) { out[a++] = PC_GAP_SENT; ii++; }
        if (i < N && !deleted[i]) out[a++] = cipher[i];
    }
    return a;
}

// One fitness evaluation of an edited stream: the exhaustive period-column search
// at inner_depth (its returned n-gram score), plus an optional word-coverage bonus
// on the n-gram-best decrypt. Writes the decrypt/stages to scratch.
static double pc_fitness(const int *edited, int M, int inner_depth,
                         float *ngram_data, int ngram_size, float weight_ngram,
                         const WordSet *ws, double weight_word,
                         int *scratch_dec, PCStage *scratch_stg, int *scratch_ns) {
    double s = period_column_search(edited, M, inner_depth, ngram_data, ngram_size,
                                    NULL, NULL, 0, weight_ngram, 0.0f,
                                    scratch_dec, scratch_stg, scratch_ns);
    if (ws && weight_word > 0.0) s += weight_word * word_coverage(scratch_dec, M, ws);
    return s;
}

// Format the composed stages "[4x42,TP,P:3][56x3,UTP,P:2]" into buf.
static void pcs_stage_str(int len, const PCStage *stages, int n_stages, char *buf, size_t n) {
    size_t off = 0;
    buf[0] = '\0';
    for (int s = 0; s < n_stages && off < n; s++) {
        int dx = stages[s].dx, dy = (len + dx - 1) / dx;
        off += snprintf(buf + off, n - off, "[%dx%d,%s,P:%d]",
                        dx, dy, stages[s].utp == 0 ? "TP" : "UTP", stages[s].period);
    }
}

// Live (-verbose) progress block for the (non-engine) edit search, in the SAME
// shape colossus prints while running an engine solve (see
// report_transposition_verbose / columnar_model_report_verbose): timing + search
// counters, a param line, then the CURRENT BEST plaintext.
static void pcs_verbose_block(clock_t t0, long n_iter, int n_restarts,
                              double score, int ni, int nd,
                              const int *decrypt, int len,
                              const PCStage *stages, int n_stages) {
    double elapsed = ((double) clock() - t0) / CLOCKS_PER_SEC;
    double ips = (elapsed > 0.) ? ((double) n_iter) / elapsed : 0.;
    char stg[128];
    pcs_stage_str(len, stages, n_stages, stg, sizeof stg);
    printf("\n%.2f\t[sec]\n", elapsed);
    printf("%.0fK\t[it/sec]\n", 1.e-3 * ips);
    printf("%d\t[restarts]\n", n_restarts);
    printf("%.4f\t[entropy]\n", entropy((int *) decrypt, len));
    printf("%.2f\t[score]\n", score);
    printf("gaps=%d drops=%d %s\t[params]\n\n", ni, nd, stg);
    print_text((int *) decrypt, len); printf("\n");
    print_solution_check((int *) decrypt, len);
    fflush(stdout);
}

double period_column_space_search(const int *cipher, int len,
    int max_ins, int max_dels, int inner_depth, int final_depth,
    int n_restarts, int n_hillclimbs, double init_temp, double min_temp,
    float *ngram_data, int ngram_size,
    float weight_ngram, float weight_word,
    SharedData *shared, bool use_words,
    int out_edited[], int *out_len,
    int out_decrypt[], PCStage stages[], int *n_stages,
    int gap_positions[], int *n_gaps,
    int del_positions[], int *n_dels, bool verbose) {

    if (max_ins < 0) max_ins = 0;
    if (max_dels < 0) max_dels = 0;
    if (max_dels > len - 4) max_dels = (len > 4) ? len - 4 : 0;   // keep >= 4 cells
    if (inner_depth < 1) inner_depth = 1;
    if (final_depth < 1) final_depth = 2;
    if (n_restarts < 1) n_restarts = 1;
    if (n_hillclimbs < 1) n_hillclimbs = 1;

    // Bound the insertion budget so the padded grid never overruns the scratch arrays.
    if (len + max_ins > MAX_CIPHER_LENGTH) max_ins = MAX_CIPHER_LENGTH - len;
    if (max_ins < 0) max_ins = 0;

    const WordSet *ws = (use_words && shared) ? trans_word_set(shared) : NULL;

    static int cur_ins[MAX_CIPHER_LENGTH], cand_ins[MAX_CIPHER_LENGTH], best_ins[MAX_CIPHER_LENGTH];
    static int cur_del[MAX_CIPHER_LENGTH], cand_del[MAX_CIPHER_LENGTH], best_del[MAX_CIPHER_LENGTH];
    static unsigned char deleted[MAX_CIPHER_LENGTH];
    static int edited[MAX_CIPHER_LENGTH], scratch_dec[MAX_CIPHER_LENGTH];
    PCStage scratch_stg[2];
    int scratch_ns = 0;

    double best_overall = -1e30;
    int best_nins = 0, best_ndel = 0;

    clock_t t0 = clock();       // -verbose timing / it-per-sec, engine-style
    long n_iter = 0;            // fitness evaluations
    int total_restarts = 0;

    if (verbose)
        printf("\nperiod-column-space: edit budgets ins 0..%d x del 0..%d, %d restart(s) x %d climb(s) each "
               "(inner depth %d, final depth %d)\n", max_ins, max_dels, n_restarts, n_hillclimbs,
               inner_depth, final_depth);

    for (int ni = 0; ni <= max_ins; ni++) {
      for (int nd = 0; nd <= max_dels; nd++) {

        if (ni == 0 && nd == 0) {
            // Baseline: the plain period-column search on the untouched cipher.
            // (0,0) makes this solver a strict superset of the plain one.
            memset(deleted, 0, (size_t) len);
            double s = pc_fitness(cipher, len, inner_depth, ngram_data, ngram_size,
                                  weight_ngram, ws, weight_word,
                                  scratch_dec, scratch_stg, &scratch_ns);
            n_iter++;
            if (s > best_overall) {
                best_overall = s; best_nins = 0; best_ndel = 0;
                if (verbose) pcs_verbose_block(t0, n_iter, total_restarts, best_overall,
                                               0, 0, scratch_dec, len, scratch_stg, scratch_ns);
            }
            continue;
        }

        double budget_best = -1e30;
        for (int rs = 0; rs < n_restarts; rs++) {
            total_restarts++;
            // From-scratch seed: ni random insertion slots + nd random deletions.
            for (int k = 0; k < ni; k++) cur_ins[k] = rand_int(0, len + 1);   // [0, len]
            gp_sort(cur_ins, ni);
            for (int k = 0; k < nd; k++) cur_del[k] = rand_int(0, len);       // [0, len-1]

            memset(deleted, 0, (size_t) len);
            for (int k = 0; k < nd; k++) deleted[cur_del[k]] = 1;
            int M = build_edited(cipher, len, cur_ins, ni, deleted, edited);
            double cur = pc_fitness(edited, M, inner_depth, ngram_data, ngram_size,
                                    weight_ngram, ws, weight_word,
                                    scratch_dec, scratch_stg, &scratch_ns);
            n_iter++;
            if (cur > budget_best) budget_best = cur;
            if (cur > best_overall) {
                best_overall = cur; best_nins = ni; best_ndel = nd;
                memcpy(best_ins, cur_ins, (size_t) ni * sizeof(int));
                memcpy(best_del, cur_del, (size_t) nd * sizeof(int));
                if (verbose) pcs_verbose_block(t0, n_iter, total_restarts, best_overall,
                                               ni, nd, scratch_dec, M, scratch_stg, scratch_ns);
            }

            for (int it = 0; it < n_hillclimbs; it++) {
                double temp = (n_hillclimbs > 1)
                    ? init_temp * pow(min_temp / init_temp, (double) it / (n_hillclimbs - 1))
                    : min_temp;

                memcpy(cand_ins, cur_ins, (size_t) ni * sizeof(int));
                memcpy(cand_del, cur_del, (size_t) nd * sizeof(int));

                // Perturb one edit: a deletion if we have any (and, when we also have
                // insertions, half the time), else an insertion slot.
                int do_del = (nd > 0) && (ni == 0 || frand() < 0.5);
                if (do_del) {
                    int k = rand_int(0, nd);                 // [0, nd-1]
                    if (frand() < 0.5) {
                        cand_del[k] = rand_int(0, len);      // global reseed
                    } else {
                        cand_del[k] += (frand() < 0.5) ? -1 : 1;
                        if (cand_del[k] < 0) cand_del[k] = 0;
                        if (cand_del[k] > len - 1) cand_del[k] = len - 1;
                    }
                } else {
                    int k = rand_int(0, ni);                 // [0, ni-1]
                    if (frand() < 0.5) {
                        cand_ins[k] = rand_int(0, len + 1);  // global reseed
                    } else {
                        cand_ins[k] += (frand() < 0.5) ? -1 : 1;
                        if (cand_ins[k] < 0) cand_ins[k] = 0;
                        if (cand_ins[k] > len) cand_ins[k] = len;
                    }
                    gp_sort(cand_ins, ni);
                }

                memset(deleted, 0, (size_t) len);
                for (int k = 0; k < nd; k++) deleted[cand_del[k]] = 1;
                int Mc = build_edited(cipher, len, cand_ins, ni, deleted, edited);
                double cand = pc_fitness(edited, Mc, inner_depth, ngram_data, ngram_size,
                                         weight_ngram, ws, weight_word,
                                         scratch_dec, scratch_stg, &scratch_ns);
                n_iter++;

                double delta = cand - cur;
                if (delta >= 0.0 || frand() < exp(delta / temp)) {
                    memcpy(cur_ins, cand_ins, (size_t) ni * sizeof(int));
                    memcpy(cur_del, cand_del, (size_t) nd * sizeof(int));
                    cur = cand;
                }
                if (cur > budget_best) budget_best = cur;
                if (cur > best_overall) {
                    best_overall = cur; best_nins = ni; best_ndel = nd;
                    memcpy(best_ins, cur_ins, (size_t) ni * sizeof(int));
                    memcpy(best_del, cur_del, (size_t) nd * sizeof(int));
                    // cur == cand here (a new best is only reached by accepting cand),
                    // so scratch_dec/Mc hold the just-evaluated best decrypt.
                    if (verbose) pcs_verbose_block(t0, n_iter, total_restarts, best_overall,
                                                   ni, nd, scratch_dec, Mc, scratch_stg, scratch_ns);
                }
            }
        }
        if (verbose) printf("  edits(ins=%d,del=%d): best-at-budget %.4f (global best %.4f)\n",
                            ni, nd, budget_best, best_overall);
      }
    }
    if (verbose) printf("  finishing best (ins=%d del=%d) at depth %d\n",
                        best_nins, best_ndel, final_depth);

    // Rebuild the best edited cipher and finish it at final_depth.
    memset(deleted, 0, (size_t) len);
    for (int k = 0; k < best_ndel; k++) deleted[best_del[k]] = 1;
    int M = build_edited(cipher, len, best_ins, best_nins, deleted, out_edited);
    *out_len = M;

    double final_score = period_column_search(out_edited, M, final_depth,
        ngram_data, ngram_size, NULL, NULL, 0, weight_ngram, 0.0f,
        out_decrypt, stages, n_stages);
    if (ws && weight_word > 0.0)
        final_score += weight_word * word_coverage(out_decrypt, M, ws);

    *n_gaps = best_nins;
    for (int k = 0; k < best_nins; k++) gap_positions[k] = best_ins[k];
    // Report the distinct deleted observed indices, ascending.
    int nd_eff = 0;
    for (int i = 0; i < len; i++) if (deleted[i]) del_positions[nd_eff++] = i;
    *n_dels = nd_eff;

    return final_score;
}

// Reverse cipher[0..len-1] into out (read the stream bottom-to-top / backwards).
static void pcs_reverse(const int *cipher, int len, int *out) {
    for (int i = 0; i < len; i++) out[i] = cipher[len - 1 - i];
}

// Format "dir=bt depth=2 gaps=1 drops=0 [4x42,TP,P:3][56x3,UTP,P:2]" for the report /
// CSV field (the "dir=bt " prefix only when the cipher was read in reverse).
static void pcs_param_summary(int len, const PCStage *stages, int n_stages,
                              int n_gaps, int n_dels, int read_dir, char *buf, size_t n) {
    size_t off = 0;
    if (read_dir == COL_READ_BT) off += snprintf(buf + off, n - off, "dir=bt ");
    off += snprintf(buf + off, n - off, "depth=%d gaps=%d drops=%d ", n_stages, n_gaps, n_dels);
    for (int s = 0; s < n_stages && off < n; s++) {
        int dx = stages[s].dx, dy = (len + dx - 1) / dx;
        off += snprintf(buf + off, n - off, "[%dx%d,%s,P:%d]",
                        dx, dy, stages[s].utp == 0 ? "TP" : "UTP", stages[s].period);
    }
}

void solve_period_column_space(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str;
    (void) cribtext_str;
    (void) crib_indices;
    (void) crib_positions;

    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a period-column-space solve.\n\n");
        return;
    }
    if (n_cribs > 0)
        printf("\nNote: period-column-space ignores cribs (inserted/deleted cells shift plaintext "
               "positions; use -type period-column for a known-clean cipher with cribs).\n");

    int final_depth = cfg->trans_depth;
    if (final_depth < 1) final_depth = 2;
    if (final_depth > 2) {
        printf("\nNote: period-column exhaustive search supports depth <= 2; clamping -depth %d to 2.\n",
               final_depth);
        final_depth = 2;
    }
    int max_ins = cfg->max_gaps;
    int max_dels = cfg->max_dels;
    if (max_ins < 0) max_ins = 0;
    if (max_dels < 0) max_dels = 0;

    static int out_edited[MAX_CIPHER_LENGTH], out_decrypt[MAX_CIPHER_LENGTH];
    static int gap_positions[MAX_CIPHER_LENGTH], del_positions[MAX_CIPHER_LENGTH];
    PCStage stages[2];
    int n_stages = 0, n_gaps = 0, n_dels = 0, out_len = cipher_len;

    // Inner-fitness depth is kept at 1 (fast) so the edit hill climb can afford many
    // evaluations; the single best edited cipher is then finished at final_depth.
    bool use_words = cfg->weight_word > 0.0f;

    // Read direction: forwards (tb, canonical) and/or the reversed stream (bt).
    // Default is COL_READ_TB, so a normal run searches the forward stream only.
    int dirs[2], ndirs = 0;
    if (cfg->read_direction == COL_READ_BOTH) { dirs[ndirs++] = COL_READ_TB; dirs[ndirs++] = COL_READ_BT; }
    else dirs[ndirs++] = cfg->read_direction;

    static int oriented[MAX_CIPHER_LENGTH];
    static int b_edited[MAX_CIPHER_LENGTH], b_dec[MAX_CIPHER_LENGTH];
    static int b_gaps[MAX_CIPHER_LENGTH], b_dels[MAX_CIPHER_LENGTH];
    PCStage b_stages[2];
    int b_ns = 0, b_ng = 0, b_nd = 0, b_len = cipher_len, best_dir = COL_READ_TB;
    double score = -1e30;

    for (int d = 0; d < ndirs; d++) {
        if (dirs[d] == COL_READ_BT) pcs_reverse(cipher_indices, cipher_len, oriented);
        else for (int i = 0; i < cipher_len; i++) oriented[i] = cipher_indices[i];

        if (cfg->verbose && ndirs > 1)
            printf("\n=== reading %s ===\n", dirs[d] == COL_READ_BT ? "bottom-to-top (bt)" : "top-to-bottom (tb)");

        double s = period_column_space_search(oriented, cipher_len,
            max_ins, max_dels, /*inner_depth=*/1, final_depth,
            cfg->n_restarts, cfg->n_hill_climbs, cfg->init_temp, cfg->min_temp,
            shared->ngram_data, cfg->ngram_size,
            cfg->weight_ngram, cfg->weight_word,
            shared, use_words,
            out_edited, &out_len, out_decrypt, stages, &n_stages,
            gap_positions, &n_gaps, del_positions, &n_dels, cfg->verbose);

        if (d == 0 || s > score) {
            score = s; best_dir = dirs[d];
            b_ns = n_stages; b_ng = n_gaps; b_nd = n_dels; b_len = out_len;
            b_stages[0] = stages[0]; b_stages[1] = stages[1];
            for (int i = 0; i < out_len; i++) { b_edited[i] = out_edited[i]; b_dec[i] = out_decrypt[i]; }
            for (int i = 0; i < n_gaps; i++) b_gaps[i] = gap_positions[i];
            for (int i = 0; i < n_dels; i++) b_dels[i] = del_positions[i];
        }
    }
    // Adopt the winning direction's results.
    n_stages = b_ns; n_gaps = b_ng; n_dels = b_nd; out_len = b_len;
    stages[0] = b_stages[0]; stages[1] = b_stages[1];
    for (int i = 0; i < out_len; i++) { out_edited[i] = b_edited[i]; out_decrypt[i] = b_dec[i]; }
    for (int i = 0; i < n_gaps; i++) gap_positions[i] = b_gaps[i];
    for (int i = 0; i < n_dels; i++) del_positions[i] = b_dels[i];

    char params[192];
    pcs_param_summary(out_len, stages, n_stages, n_gaps, n_dels, best_dir, params, sizeof params);

    printf("\nperiod-column-space: %d composed stage(s), read %s; %d inserted gap(s), %d deleted cell(s) "
           "(observed length %d -> edited %d; widths are complete-grid divisors of %d)\n",
           n_stages, best_dir == COL_READ_BT ? "bottom-to-top" : "top-to-bottom",
           n_gaps, n_dels, cipher_len, out_len, out_len);
    for (int s = 0; s < n_stages; s++) {
        int dx = stages[s].dx, dy = (out_len + dx - 1) / dx;
        printf("  stage %d: %d x %d grid, period %d, %s\n", s + 1, dx, dy,
               stages[s].period, stages[s].utp == 0 ? "transpose" : "untranspose");
    }
    if (n_gaps > 0) {
        printf("  gaps inserted before observed cipher position(s):");
        for (int g = 0; g < n_gaps; g++) printf(" %d", gap_positions[g]);
        printf("\n");
    }
    if (n_dels > 0) {
        printf("  observed cipher cell(s) deleted at position(s):");
        for (int d = 0; d < n_dels; d++) printf(" %d", del_positions[d]);
        printf("\n");
    }

    // Report against the EDITED cipher (length out_len) so the echo and decrypt line
    // up cell-for-cell and the inserted gaps show as spaces. Cribs are not used, so a
    // fresh all-'_' crib row of the edited length keeps report_transposition happy.
    static char cribM[MAX_CIPHER_LENGTH];
    for (int i = 0; i < out_len; i++) cribM[i] = '_';
    cribM[out_len] = '\0';

    report_transposition(cfg, shared, out_edited, out_len, out_decrypt,
        score, cribM, 0, params);
}
