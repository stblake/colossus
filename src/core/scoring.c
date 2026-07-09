#include "scoring.h"

// Crib-dragging scoring globals (see -cribdrag). Default off => state_score is
// bit-identical to the pre-feature behaviour. Set from main() once the words are
// parsed, mirroring the g_ngram_logprob / g_ngram_reverse scoring-toggle pattern.
const CribDrag *g_cribdrag = NULL;
float g_cribdrag_weight = 0.0f;

double state_score(int decrypted[], int cipher_len,
            int crib_indices[], int crib_positions[], int n_cribs,
            float *ngram_data, int ngram_size,
            float weight_ngram, float weight_crib,
            float weight_ioc, float weight_entropy) {

    double score, decrypted_ngram_score = 0., decrypted_crib_score = 0., drag_score = 0.;

    int have_drag = (g_cribdrag && g_cribdrag_weight > 1.e-4 && g_cribdrag->nwords > 0);

    if (weight_crib > 1.e-4) {
        decrypted_crib_score = crib_score(decrypted, cipher_len, crib_indices, crib_positions, n_cribs);
    }

    if (weight_ngram > 1.e-4) {
        decrypted_ngram_score = ngram_score(decrypted, cipher_len, ngram_data, ngram_size);
    }

    if (have_drag) {
        drag_score = cribdrag_score(decrypted, cipher_len, g_cribdrag);
    }

    if (n_cribs > 0 || have_drag) {
        double num = weight_ngram * decrypted_ngram_score;
        double den = weight_ngram;
        if (n_cribs > 0)  { num += weight_crib * decrypted_crib_score; den += weight_crib; }
        if (have_drag)    { num += g_cribdrag_weight * drag_score;      den += g_cribdrag_weight; }
        score = num / den;
    } else {
        score = decrypted_ngram_score;
    }

    return score;
}

// Crib dragging: each supplied word is slid across the decrypt and its BEST-matching
// offset scored (max over offsets). The reward is the mean of the per-word bests
// (== AND: every word is expected to appear somewhere). Per-letter matching mirrors
// crib_score's PARTIAL_CRIB_MATCH convention so the hill-climber gets a gradient.
double cribdrag_score(int text[], int len, const CribDrag *cd) {
    if (!cd || cd->nwords <= 0) return 0.;
    double total = 0.;
    for (int w = 0; w < cd->nwords; w++) {
        int L = cd->wordlen[w];
        if (L <= 0 || L > len) continue;   // word can't fit -> contributes 0
        const int *word = cd->words[w];
        double best = 0.;
        for (int off = 0; off + L <= len; off++) {
            double m = 0.;
            for (int j = 0; j < L; j++) {
                int diff = abs(text[off + j] - word[j]);
#if PARTIAL_CRIB_MATCH
                m += (diff == 0) ? 1. : 1./(1. + diff * diff);
#else
                m += (diff == 0) ? 1. : 0.;
#endif
            }
            m /= (double) L;
            if (m > best) {
                best = m;
                if (best >= 1.) break;     // exact hit: no better offset possible
            }
        }
        total += best;
    }
    return total / ((double) cd->nwords);
}



double crib_score(int text[], int len, int crib_indices[], int crib_positions[], int n_cribs) {
    if (n_cribs == 0) return 0.;
#if PARTIAL_CRIB_MATCH
    int diff;
    double score = 0.;
    for (int i = 0; i < n_cribs; i++) {
        diff = abs(text[crib_positions[i]] - crib_indices[i]);
        if (diff == 0) {
            score += 1.;
        } else {
            score += 1./(1. + diff * diff);
        }
    }
    return score / ((double) n_cribs);
#else
    int n_matches = 0;
    for (int i = 0; i < n_cribs; i++) {
        if (text[crib_positions[i]] == crib_indices[i]) {
            n_matches += 1;
        }
    }
    return ((double) n_matches)/((double) n_cribs);
#endif
}

double ngram_score(int decrypted[], int cipher_len, float *ngram_data, int ngram_size) {
    int index, base;
    double score = 0.;

    // pow(g_alpha, ngram_size) is a positive constant for the whole run
    // (ngram_size never changes), yet was previously recomputed via a libm pow()
    // on EVERY score -- i.e. every hill-climber iteration. Memoize it. pow()
    // returns the identical double for identical args, so the cached value equals
    // the recomputed one bit-for-bit; the score is unchanged.
    static int cached_ngram_size = -1;
    static double scale = 0.;
    if (ngram_size != cached_ngram_size) {
        // Legacy table entries are ~1/n_ngrams, so the historical g_alpha^ngram_size
        // factor brings the mean back to O(1). The log-prob table already holds O(1)
        // log10 values, so it needs no rescaling (scale = 1) -- the score is then a
        // mean log-probability, the AZDecrypt fitness.
        scale = g_ngram_logprob ? 1.0 : pow(g_alpha, ngram_size);
        cached_ngram_size = ngram_size;
    }

    // Compact the decrypted text to its LETTERS ONLY, dropping negative sentinels
    // (spaces / punctuation carried through from the ciphertext). n-grams are then
    // formed from CONSECUTIVE LETTERS, so a sentinel is transparent to the window
    // ("THE MOST" -> THEMO, HEMOS, ...) rather than voiding every window that spans
    // it. This is the fix for space-bearing transpositions: with a space every few
    // characters the old "skip any window containing a sentinel" rule discarded
    // almost every window, flattening the fitness landscape so the climber had no
    // gradient (AZDecrypt-style solvers score the space as a symbol instead). Word-
    // boundary (space-placement) signal is supplied separately by the dictionary
    // word-coverage term (-weightword). The sentinel COUNT is invariant under a
    // transposition (the same multiset is permuted), so `m` is constant across every
    // candidate decrypt of a given cipher -- the normalisation below never changes
    // which solution wins. When the text is all letters the compaction is the
    // identity, `m == cipher_len`, and every step is bit-for-bit identical to the
    // historical scorer (same indices, same additions, same order, same divisor).
    static int letters[MAX_CIPHER_LENGTH];
    int m = 0;
    for (int i = 0; i < cipher_len; i++)
        if (decrypted[i] >= 0) letters[m++] = decrypted[i];

    // Rolling base-26 index over the compacted stream. The packed window index is
    // little-endian idx_i = sum_{j} letters[i+j] * 26^j, so advancing one position
    // is exact integer arithmetic:
    //   idx_{i+1} = (idx_i - letters[i]) / 26 + letters[i+n] * 26^(n-1).
    // (idx_i - letters[i]) is divisible by 26, so the integer division is exact.
    int n_windows = m - ngram_size + 1;
    if (n_windows > 0) {
        int top = 1;                    // 26^(ngram_size-1)
        for (int j = 0; j < ngram_size - 1; j++) top *= g_alpha;

        index = 0;
        base = 1;
        for (int j = 0; j < ngram_size; j++) {
            index += letters[j]*base;
            base *= g_alpha;
        }
        score += ngram_data[index];

        for (int i = 1; i < n_windows; i++) {
            index = (index - letters[i - 1]) / g_alpha + letters[i + ngram_size - 1] * top;
            score += ngram_data[index];
        }
    }
    int denom = m - ngram_size;
    score = (denom > 0) ? scale*score/denom : 0.0;
    return score;
}

double ngram_sum_raw(const int *text, int len, const float *ngram_data, int ngram_size) {
    double score = 0.;
    int n_windows = len - ngram_size + 1;
    if (n_windows <= 0) return 0.;

    int top = 1;                          // g_alpha^(ngram_size-1)
    for (int j = 0; j < ngram_size - 1; j++) top *= g_alpha;

    int index = 0, base = 1, bad = 0;
    for (int j = 0; j < ngram_size; j++) {
        int v = text[j];
        if (v < 0) { bad++; v = 0; }       // sentinel: contributes 0 to the packed index
        index += v * base;
        base *= g_alpha;
    }
    if (bad == 0) score += ngram_data[index];

    for (int i = 1; i < n_windows; i++) {
        int out_v = text[i - 1];
        int in_v  = text[i + ngram_size - 1];
        if (out_v < 0) { bad--; out_v = 0; }
        int in_iv = in_v;
        if (in_v < 0) { bad++; in_iv = 0; }
        index = (index - out_v) / g_alpha + in_iv * top;
        if (bad == 0) score += ngram_data[index];
    }
    return score;
}

void perturbate_cycleword(int state[], int max, int len) {
    int i = rand_int(0, len);
    state[i] = rand_int(0, max);
}

void perturbate_keyword(int state[], int len, int keyword_len) {
    int i, j, k, l, temp;

    if (frand() < 0.2) { 
        // Swap two letters of the key.
        i = rand_int(0, keyword_len);
        j = rand_int(0, keyword_len);
        temp = state[i];
        state[i] = state[j];
        state[j] = temp;
    } else {
        // Swap a letter from the key with an alphabet letter. 
#if FREQUENCY_WEIGHTED_SELECTION
        i = rand_int_frequency_weighted(state, 0, keyword_len);
        j = rand_int_frequency_weighted(state, keyword_len, len);
#else
        i = rand_int(0, keyword_len);
        j = rand_int(keyword_len, len);
#endif
        temp = state[i];
        state[i] = state[j];
        for (k = j + 1; k < len; k++) state[k - 1] = state[k];
        for (k = keyword_len; k < len; k++) {
            if (state[k] > temp || k == len - 1) {
                for (l = len - 1; l > k; l--) state[l] = state[l - 1];
                state[k] = temp;
                break ;
            }
        }
    }
}

void random_keyword(int keyword[], int len, int keyword_len) {
    int i, j, candidate, indx, n_chars;
    bool distinct, present;
    n_chars = 0;
    while (n_chars < keyword_len) {
        distinct = true;
        candidate = rand_int(0, g_alpha);
        for (i = 0; i < n_chars; i++) {
            if (keyword[i] == candidate) {
                distinct = false;
                break ;
            }
        }
        if (distinct) keyword[n_chars++] = candidate;
    }
    indx = keyword_len;
    for (i = 0; i < g_alpha; i++) {
        present = false;
        for (j = 0; j < keyword_len; j++) {
            if (keyword[j] == i) {
                present = true; 
                break ;
            }
        }
        if (! present) keyword[indx++] = i;
    }
}

void random_cycleword(int cycleword[], int max, int keyword_len) {
    for (int i = 0; i < keyword_len; i++) {
        cycleword[i] = rand_int(0, max);
    }
}

int rand_int_frequency_weighted(int state[], int min_index, int max_index) {
    double total = 0.0;
    double cumsum = 0.0;

    for (int i = min_index; i < max_index; i++) {
        total += english_monograms[state[i]];
    }

    if (total == 0.0) {
        return rand_int(min_index, max_index - 1); 
    }

    // Multiply the random float [0.0, 1.0) by the total weight.
    double target = frand() * total; 

    // Accumulate raw weights.
    for (int i = min_index; i < max_index; i++) {
        cumsum += english_monograms[state[i]];
        if (cumsum >= target) {
            return i;
        }
    }

    return max_index - 1;
}

float* load_ngrams(char *ngram_file, int ngram_size, bool verbose) {
    FILE *fp;
    int i, n_ngrams, freq, indx;
    char ngram[MAX_NGRAM_SIZE];
    float *ngram_data, total;

    if (verbose) printf("\nLoading ngrams...");
    n_ngrams = int_pow(g_alpha, ngram_size);
    ngram_data = malloc(n_ngrams*sizeof(float));
    for (i = 0; i < n_ngrams; i++) ngram_data[i] = 0.;

    fp = fopen(ngram_file, "r");
    // Loop on the parse succeeding (both fields read), not on feof: !feof is
    // still false after the last good line, so feof-looping re-reads the final
    // line and would mis-assign on any trailing/malformed line.
    while (fscanf(fp, "%s\t%d", ngram, &freq) == 2) {
        indx = ngram_index_str(ngram, ngram_size);
        if (indx < 0) continue;   // n-gram uses a letter not in the runtime alphabet
        ngram_data[indx] = freq;
    }
    fclose(fp);

    if (g_ngram_logprob) {
        // AZDecrypt / Practical-Cryptography fitness: each cell holds log10 P(n-gram),
        // and every UNSEEN n-gram is set to a floor probability so implausible n-grams
        // are penalised (the legacy table leaves them at 0, i.e. merely unrewarded).
        // The per-window sum of these log-probs is the standard n-gram fitness; ngram_score
        // keeps the scale at 1 in this mode so the result is a mean log-probability.
        double count_total = 0.;
        for (i = 0; i < n_ngrams; i++) count_total += ngram_data[i];   // raw counts
        if (count_total <= 0.) count_total = 1.;
        double floor = log10(0.01 / count_total);   // ~ a rare-but-not-impossible n-gram
        for (i = 0; i < n_ngrams; i++)
            ngram_data[i] = (ngram_data[i] > 0.) ? (float) log10(ngram_data[i] / count_total)
                                                 : (float) floor;
    } else {
        // Legacy reward-only scheme: normalized log(1 + count); unseen -> 0.
        total = 0.;
        for (i = 0; i < n_ngrams; i++) {
            ngram_data[i] = log(1. + ngram_data[i]);
            total += ngram_data[i];
        }
        for (i = 0; i < n_ngrams; i++) ngram_data[i] /= total;
    }

    if (g_ngram_reverse) {
        // Reversal-invariant scoring: give each n-gram and its digit-reversed twin the
        // same (max) weight, so a plaintext written with any words/segments reversed
        // scores like clean English (a reversed word's n-grams are the reverses of the
        // forward word's). Applied AFTER normalization so weights are final; max makes
        // it symmetric and order-independent. The packed index is
        // sum digit[j]*g_alpha^j (position 0 = least significant, per ngram_score), so
        // the reversed n-gram maps digit at j to weight g_alpha^(ngram_size-1-j).
        float *sym = malloc(n_ngrams*sizeof(float));
        int topw = int_pow(g_alpha, ngram_size - 1);   // g_alpha^(ngram_size-1)
        for (i = 0; i < n_ngrams; i++) {
            int idx = i, rev = 0, hi = topw, j;
            for (j = 0; j < ngram_size; j++) {
                rev += (idx % g_alpha) * hi;
                idx /= g_alpha;
                hi  /= g_alpha;
            }
            sym[i] = (ngram_data[i] > ngram_data[rev]) ? ngram_data[i] : ngram_data[rev];
        }
        free(ngram_data);
        ngram_data = sym;
    }

    if (verbose) printf("...finished.\n\n");
    return ngram_data;
}

int ngram_index_str(char *ngram, int ngram_size) {
    int c, index = 0, base = 1;
    for (int i = 0; i < ngram_size; i++) {
        c = g_char_to_idx[toupper((unsigned char) ngram[i]) & 127];
        // An n-gram containing a letter outside the runtime alphabet (e.g. 'P'
        // under -excludeletter P) cannot occur in the plaintext, so it has no
        // slot; signal the caller to skip it.
        if (c < 0) return -1;
        index += c*base;
        base *= g_alpha;
    }
    return index;
}

int ngram_index_int(int *ngram, int ngram_size) {
    int index = 0, base = 1;
    for (int i = 0; i < ngram_size; i++) {
        index += ngram[i]*base;
        base *= g_alpha;
    }
    return index;
}
