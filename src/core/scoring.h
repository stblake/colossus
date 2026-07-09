#ifndef SCORING_H
#define SCORING_H
#include "colossus.h"

// n-gram / crib fitness and the keyword/cycleword randomization helpers shared by
// the engine and the polyalphabetic model.
double state_score(int decrypted[], int cipher_len,
            int crib_indices[], int crib_positions[], int n_cribs,
            float *ngram_data, int ngram_size,
            float weight_ngram, float weight_crib,
            float weight_ioc, float weight_entropy);
double ngram_score(int decrypted[], int cipher_len, float *ngram_data, int ngram_size);

// Raw (un-normalized) within-word n-gram sum: the same window walk as ngram_score
// -- windows that straddle a negative sentinel (space/punct) are skipped -- but
// returns the bare sum of ngram_data over the surviving windows, with NO division
// by length and NO g_alpha^n scale. Because it is a plain additive sum, the score
// of two concatenated texts equals the sum of their individual scores plus the
// windows spanning the join, which is what makes an exact seam decomposition (and
// hence Held-Karp best row ordering) possible. Not used by any existing solve.
double ngram_sum_raw(const int *text, int len, const float *ngram_data, int ngram_size);
double crib_score(int text[], int len, int crib_indices[], int crib_positions[], int n_cribs);
// Crib dragging: mean over words of each word's best-offset partial match (see -cribdrag).
double cribdrag_score(int text[], int len, const CribDrag *cd);

float* load_ngrams(char *ngram_file, int ngram_size, bool verbose);
int ngram_index_int(int *ngram, int ngram_size);
int ngram_index_str(char *ngram, int ngram_size);

void perturbate_keyword(int state[], int len, int keyword_len);
void random_keyword(int keyword[], int len, int keyword_len);
void random_cycleword(int cycleword[], int max, int keyword_len);
void perturbate_cycleword(int state[], int max, int len);
int rand_int_frequency_weighted(int state[], int min_index, int max_index);
#endif
