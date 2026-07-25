#include "indep_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Independent periodic substitution (TYPE indep_periodic)
// =====================================================================
//
// A period-P cipher in which each residue class i % P is enciphered with its OWN,
// INDEPENDENT mixed substitution alphabet (unlike Quagmire, whose columns are all
// shifts of one keyed alphabet). decrypted[i] = map[i % P][cipher[i]], where each
// map[] is a permutation of the runtime alphabet sending a cipher letter to its
// plaintext letter.
//
// The columns CANNOT be solved one at a time: column j read alone is every P-th
// letter of the message and carries no n-gram signal. The only constraint tying
// the P alphabets down is that CONSECUTIVE plaintext letters -- which come from
// DIFFERENT columns -- must form English n-grams. So all P alphabets are climbed
// JOINTLY against the reassembled text. With P alphabets x (alpha-1) d.o.f. this is
// a hard search, so we use simulated annealing + iterated local search (restarts
// that perturb the global best, not just fresh seeds) and frequency-seeded starts.

// Build the plaintext (as indices) from the P per-column maps. Sentinels (<0,
// carried-through spaces/punctuation) pass through unchanged.
static void indep_decrypt(int cipher_indices[], int cipher_len, int period,
    int maps[][ALPHABET_SIZE], int decrypted[]) {
    for (int i = 0; i < cipher_len; i++) {
        int c = cipher_indices[i];
        decrypted[i] = (c >= 0) ? maps[i % period][c] : c;
    }
}

// Frequency-seed each column's map: most-frequent cipher letter in the column ->
// most-frequent English letter, etc. Gives the annealer a sensible starting point.
static void indep_seed(int cipher_indices[], int cipher_len, int period,
    int maps[][ALPHABET_SIZE]) {

    // English letters (runtime alphabet indices) sorted by descending frequency.
    int eng_rank[ALPHABET_SIZE];
    for (int i = 0; i < g_alpha; i++) eng_rank[i] = i;
    for (int a = 0; a < g_alpha; a++)
        for (int b = a + 1; b < g_alpha; b++)
            if (g_monograms[eng_rank[b]] > g_monograms[eng_rank[a]]) {
                int t = eng_rank[a]; eng_rank[a] = eng_rank[b]; eng_rank[b] = t;
            }

    for (int j = 0; j < period; j++) {
        int hist[ALPHABET_SIZE];
        for (int c = 0; c < g_alpha; c++) hist[c] = 0;
        for (int i = j; i < cipher_len; i += period)
            if (cipher_indices[i] >= 0) hist[cipher_indices[i]]++;
        // cipher letters sorted by descending column frequency (stable on ties)
        int crank[ALPHABET_SIZE];
        for (int c = 0; c < g_alpha; c++) crank[c] = c;
        for (int a = 0; a < g_alpha; a++)
            for (int b = a + 1; b < g_alpha; b++)
                if (hist[crank[b]] > hist[crank[a]]) {
                    int t = crank[a]; crank[a] = crank[b]; crank[b] = t;
                }
        for (int r = 0; r < g_alpha; r++) maps[j][crank[r]] = eng_rank[r];
    }
}

// Coordinate ascent: fully local-optimize ONE column's alphabet (best-improving
// pairwise swaps to convergence) holding the other columns fixed. This is the move
// that breaks the coordination barrier -- single random swaps barely move the score
// while neighbouring columns are still wrong, but greedily perfecting one column
// against the (partially correct) others gives a real gradient. Returns the score.
static double indep_column_opt(ColossusConfig *cfg,
    int cipher_indices[], int cipher_len, int period, int j,
    int maps[][ALPHABET_SIZE], int crib_indices[], int crib_positions[], int n_cribs,
    float *ngram_data, int decrypted[], double cur_score) {

    bool improved = true;
    while (improved) {
        improved = false;
        for (int a = 0; a < g_alpha; a++) {
            for (int b = a + 1; b < g_alpha; b++) {
                int t = maps[j][a]; maps[j][a] = maps[j][b]; maps[j][b] = t;
                indep_decrypt(cipher_indices, cipher_len, period, maps, decrypted);
                double sc = state_score(decrypted, cipher_len,
                    crib_indices, crib_positions, n_cribs, ngram_data, cfg->ngram_size,
                    cfg->weight_ngram, cfg->weight_crib, cfg->weight_ioc, cfg->weight_entropy);
                if (sc > cur_score) { cur_score = sc; improved = true; }
                else { t = maps[j][a]; maps[j][a] = maps[j][b]; maps[j][b] = t; }
            }
        }
    }
    return cur_score;
}

static double shotgun_indep_climber(ColossusConfig *cfg,
    int cipher_indices[], int cipher_len, int period,
    int crib_indices[], int crib_positions[], int n_cribs,
    float *ngram_data, int best_decrypted[], int best_maps[][ALPHABET_SIZE]) {

    int cur[MAX_COLS][ALPHABET_SIZE], loc[MAX_COLS][ALPHABET_SIZE];
    int seed[MAX_COLS][ALPHABET_SIZE];
    int decrypted[MAX_CIPHER_LENGTH];
    double best_score = 0., current_score = 0.;
    bool have_best = false;
    size_t state_bytes = (size_t) period * ALPHABET_SIZE * sizeof(int);
    (void) loc;

    indep_seed(cipher_indices, cipher_len, period, seed);

    clock_t start_time = clock();
    long n_iterations = 0, n_slips = 0, n_backtracks = 0;
    double elapsed, n_iter_per_sec, entropy_score;

    for (long rs = 0; rs < cfg->n_restarts; rs++) {

        if (have_best && frand() < cfg->backtracking_probability) {
            // Iterated local search: perturb the global best (basin hopping).
            memcpy(cur, best_maps, state_bytes);
            int kicks = rand_int(3, 9);
            for (int k = 0; k < kicks; k++) {
                int j = rand_int(0, period), a = rand_int(0, g_alpha), b = rand_int(0, g_alpha);
                int t = cur[j][a]; cur[j][a] = cur[j][b]; cur[j][b] = t;
            }
            n_backtracks += 1;
        } else {
            // Fresh frequency-seeded start, lightly shuffled for diversity.
            memcpy(cur, seed, state_bytes);
            for (int k = 0; k < (int)(rs % 25); k++) {
                int j = rand_int(0, period), a = rand_int(0, g_alpha), b = rand_int(0, g_alpha);
                int t = cur[j][a]; cur[j][a] = cur[j][b]; cur[j][b] = t;
            }
        }
        (void) n_slips;
        indep_decrypt(cipher_indices, cipher_len, period, cur, decrypted);
        current_score = state_score(decrypted, cipher_len,
            crib_indices, crib_positions, n_cribs, ngram_data, cfg->ngram_size,
            cfg->weight_ngram, cfg->weight_crib, cfg->weight_ioc, cfg->weight_entropy);

        // Local search = coordinate-ascent sweeps over the columns: optimize each
        // column fully given the others, repeating until a whole sweep makes no
        // progress. Columns are visited in random order each sweep. n_hill_climbs
        // caps the number of sweeps.
        int order[MAX_COLS];
        for (int j = 0; j < period; j++) order[j] = j;
        int max_sweeps = cfg->n_hill_climbs > 0 ? cfg->n_hill_climbs : 1;
        if (max_sweeps > 200) max_sweeps = 200;   // coordinate ascent converges fast
        for (int sweep = 0; sweep < max_sweeps; sweep++) {
            for (int x = period - 1; x > 0; x--) {   // shuffle column order
                int y = rand_int(0, x + 1); int t = order[x]; order[x] = order[y]; order[y] = t;
            }
            double before = current_score;
            for (int jj = 0; jj < period; jj++) {
                n_iterations += 1;
                current_score = indep_column_opt(cfg, cipher_indices, cipher_len, period,
                    order[jj], cur, crib_indices, crib_positions, n_cribs,
                    ngram_data, decrypted, current_score);
            }
            if (current_score <= before + 1.e-9) break;   // converged
        }

        if (!have_best || current_score > best_score) {
            best_score = current_score;
            memcpy(best_maps, cur, state_bytes);
            have_best = true;
            if (cfg->verbose) {
                indep_decrypt(cipher_indices, cipher_len, period, best_maps, decrypted);
                entropy_score = entropy(decrypted, cipher_len);
                elapsed = ((double) clock() - start_time)/CLOCKS_PER_SEC;
                n_iter_per_sec = (elapsed > 0.) ? ((double) n_iterations)/elapsed : 0.;
                printf("\n%.2f\t[sec]\n%.0fK\t[col-opts/sec]\n%ld\t[restarts]\n%ld\t[backtracks]\n"
                       "%.4f\t[entropy]\nperiod %d\t[params]\n%.2f\t[score]\n",
                    elapsed, 1.e-3*n_iter_per_sec, rs, n_backtracks,
                    entropy_score, period, best_score);
                // The N independent alphabets (each row maps cipher A.. -> plaintext),
                // mirroring how the Quagmire climber prints its keyed alphabets.
                for (int j = 0; j < period; j++) {
                    printf("alphabet %d: ", j);
                    for (int c = 0; c < g_alpha; c++) printf("%c", index_to_char(best_maps[j][c]));
                    printf("\n");
                }
                printf("\n");
                print_text(decrypted, cipher_len); printf("\n");
                print_solution_check(decrypted, cipher_len);
                fflush(stdout);
            }
        }
    }

    indep_decrypt(cipher_indices, cipher_len, period, best_maps, best_decrypted);
    return best_score;
}

void solve_indep_periodic(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str;

    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for an independent-periodic solve.\n\n");
        return ;
    }

    // Candidate periods. With -cyclewordlen, use it; otherwise estimate them by
    // columnar-IoC the same way the periodic polyalphabetic ciphers do (each
    // residue class mod P is monoalphabetic, so the IoC peaks at the true P).
    int periods[MAX_CYCLEWORD_LEN];
    int n_periods = 0;
    if (cfg->cycleword_len_present) {
        periods[n_periods++] = cfg->cycleword_len;
    } else {
        estimate_cycleword_lengths(cipher_indices, cipher_len,
            cfg->max_cycleword_len, cfg->n_sigma_threshold, cfg->ioc_threshold,
            &n_periods, periods, cfg->verbose);
        if (n_periods == 0) {
            printf("\nNo periodicities found above threshold; nothing to attack.\n");
            return ;
        }
    }

    int best_decrypted[MAX_CIPHER_LENGTH], best_maps[MAX_COLS][ALPHABET_SIZE];
    int try_decrypted[MAX_CIPHER_LENGTH], try_maps[MAX_COLS][ALPHABET_SIZE];
    double best_score = -1.e18;
    int best_period = periods[0];

    for (int pi = 0; pi < n_periods; pi++) {
        int p = periods[pi];
        if (p < 1 || p > MAX_COLS || p > cipher_len / 2) continue;
        double sc = shotgun_indep_climber(cfg, cipher_indices, cipher_len, p,
            crib_indices, crib_positions, n_cribs,
            shared->ngram_data, try_decrypted, try_maps);
        if (cfg->verbose)
            printf("\nperiod %d: score %.2f\n", p, sc);
        if (sc > best_score) {
            best_score = sc;
            best_period = p;
            memcpy(best_decrypted, try_decrypted, (size_t) cipher_len * sizeof(int));
            memcpy(best_maps, try_maps, (size_t) p * ALPHABET_SIZE * sizeof(int));
        }
    }

    char param_summary[64];
    snprintf(param_summary, sizeof(param_summary), "period=%d", best_period);
    report_transposition(cfg, shared, cipher_indices, cipher_len, best_decrypted,
        best_score, cribtext_str, n_cribs, param_summary);

    // Recovered per-column alphabets (cipher A..Z -> plaintext), for reproduction.
    for (int j = 0; j < best_period; j++) {
        printf("col %d:", j);
        for (int c = 0; c < g_alpha; c++) printf("%c", index_to_char(best_maps[j][c]));
        printf("\n");
    }
}


