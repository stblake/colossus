#include "scoring.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <limits.h>

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
        // AZDecrypt-style multiplicative entropy coupling (opt-in, -weightentropy > 0):
        //   fitness_ngram = (ngram_score - floor) * H^weight_entropy
        // where H is the Shannon entropy of the decrypted letter distribution. Pure
        // n-gram fitness barely separates the true plaintext from the homophonic
        // "collapse" fixed point (fold many symbols onto E/T/A to tile common n-grams);
        // that collapse has a low-entropy letter distribution, so scaling the n-gram
        // magnitude by H^w suppresses it PROPORTIONALLY to its own n-gram fit (scale-
        // invariant, unlike an additive monogram penalty). Shifting by g_ngram_floor
        // makes the log-prob magnitude non-negative so the multiply is correctly
        // oriented (higher entropy -> larger score). Default weight_entropy == 0 leaves
        // decrypted_ngram_score untouched -> every existing solve is bit-identical.
        if (weight_entropy > 1.e-4) {
            double H = entropy(decrypted, cipher_len);
            decrypted_ngram_score = (decrypted_ngram_score - g_ngram_floor)
                                    * pow(H, weight_entropy);
        }
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

// Rolling big-endian window sum over `src[0..m-1]`. Kept as a static inline so ngram_score
// can call it with COMPILE-TIME `ngram_size`/`alpha`/`top` for the common alphabets (25 for
// the J-merged squares, 26 for the plain types), letting the compiler strength-reduce the
// base multiplies and fully unroll the fixed-size init loop, and with runtime values for the
// 27/36 alphabets. The arithmetic (and hence every packed index, gather, and addition order)
// is identical either way -- bit-for-bit. (A software-prefetch variant of the loop-carried
// gather was tried here and REJECTED: on the only path where it fires -- the memory-bound
// quintgram table -- it measured ~5% SLOWER, because the recurrence already runs address
// generation far enough ahead of the result-independent gathers to saturate memory-level
// parallelism, so the extra look-ahead index arithmetic is pure overhead.)
static inline double ngram_walk(const int * restrict src, int m,
                                const float * restrict nd,
                                int ngram_size, int alpha, int top) {
    double score = 0.;
    int n_windows = m - ngram_size + 1;
    if (n_windows > 0) {
        int index = 0;
        for (int j = 0; j < ngram_size; j++)
            index = index * alpha + src[j];
        score += nd[index];

        for (int i = 1; i < n_windows; i++) {
            index = (index - src[i - 1] * top) * alpha + src[i + ngram_size - 1];
            score += nd[index];
        }
    }
    return score;
}

// Byte-table twin of ngram_walk (see the -reversengrams/.ngbin path): the packed
// index is computed IDENTICALLY (same big-endian roll), but the per-window gather
// dereferences the dense 8-bit table `nd8` and dequantizes through the 256-entry
// `lut` (lut[b] == w_floor + b*w_scale). Same specialization scheme as ngram_walk so
// the base multiplies strength-reduce for the hot (alpha, size) pairs. Only reached
// when g_ngram_u8 != NULL, so the float path above stays bit-for-bit unchanged.
static inline double ngram_walk_u8(const int * restrict src, int m,
                                   const unsigned char * restrict nd8,
                                   const float * restrict lut,
                                   int ngram_size, int alpha, int top) {
    double score = 0.;
    int n_windows = m - ngram_size + 1;
    if (n_windows > 0) {
        int index = 0;
        for (int j = 0; j < ngram_size; j++)
            index = index * alpha + src[j];
        score += lut[nd8[index]];

        for (int i = 1; i < n_windows; i++) {
            index = (index - src[i - 1] * top) * alpha + src[i + ngram_size - 1];
            score += lut[nd8[index]];
        }
    }
    return score;
}

double ngram_score(int decrypted[], int cipher_len, float *ngram_data, int ngram_size) {
    double score = 0.;

    // pow(g_alpha, ngram_size) and top = g_alpha^(ngram_size-1) are positive constants
    // for the whole run (g_alpha and ngram_size never change), yet were recomputed on
    // EVERY score -- i.e. every hill-climber iteration (pow() via libm, top via a small
    // loop). Memoize both, keyed on ngram_size (g_alpha is fixed per run, as the cached
    // scale already assumes). The cached values equal the recomputed ones bit-for-bit,
    // so the score is unchanged.
    static _Thread_local int cached_ngram_size = -1;
    static _Thread_local double scale = 0.;
    static _Thread_local int cached_top = 1;
    if (ngram_size != cached_ngram_size) {
        // Legacy table entries are ~1/n_ngrams, so the historical g_alpha^ngram_size
        // factor brings the mean back to O(1). The log-prob table already holds O(1)
        // log10 values, so it needs no rescaling (scale = 1) -- the score is then a
        // mean log-probability, the AZDecrypt fitness.
        scale = g_ngram_logprob ? 1.0 : pow(g_alpha, ngram_size);
        int t = 1;                      // g_alpha^(ngram_size-1)
        for (int j = 0; j < ngram_size - 1; j++) t *= g_alpha;
        cached_top = t;
        cached_ngram_size = ngram_size;
    }
    const int top = cached_top;

    // Source stream for the window walk. n-grams are formed from CONSECUTIVE LETTERS, so
    // a negative sentinel (space / punctuation carried through from the ciphertext) must
    // be transparent to the window ("THE MOST" -> THEMO, HEMOS, ...) rather than voiding
    // every window that spans it -- the fix for space-bearing transpositions, where the
    // old "skip any window containing a sentinel" rule discarded almost every window and
    // flattened the gradient (AZDecrypt-style solvers score the space as a symbol; word-
    // boundary signal is supplied separately by -weightword). When the cipher is all
    // letters (g_score_no_sentinel, set once after decode) the compaction is the identity
    // -- `m == cipher_len` and letters[i] == decrypted[i] -- so we score decrypted[] in
    // place, saving an O(n) copy+branch pass per iteration; every index, addition, order
    // and divisor is then bit-for-bit identical to compacting first. Otherwise compact to
    // letters only. The sentinel COUNT is invariant under a transposition (the same
    // multiset is permuted), so `m` is constant across every candidate decrypt of a given
    // cipher and the normalisation below never changes which solution wins.
    const int * restrict src;
    int m;
    if (g_score_no_sentinel) {
        src = decrypted;
        m = cipher_len;
    } else {
        static _Thread_local int letters[MAX_CIPHER_LENGTH];
        int k = 0;
        for (int i = 0; i < cipher_len; i++)
            if (decrypted[i] >= 0) letters[k++] = decrypted[i];
        src = letters;
        m = k;
    }

    // Rolling base-26 index over the stream. The packed window index is BIG-endian
    // idx_i = sum_{j} src[i+j] * 26^(n-1-j) (position 0 = most significant), so advancing
    // one position is a MULTIPLY, not a divide:
    //   idx_{i+1} = (idx_i - src[i]*26^(n-1)) * 26 + src[i+n].
    // Dropping the divide matters because idx is loop-carried (each window's index depends
    // on the previous), so an integer division by the runtime variable g_alpha would
    // serialize a ~20-cycle divide per window; the multiply is ~3. Bit-identical to the
    // old little-endian packing: load_ngrams builds the table with the SAME convention
    // (ngram_index_str), so every n-gram lands the SAME float at a different slot, and the
    // per-window sum runs in the SAME order.
    // The window walk is factored into ngram_walk (below) so the hot (g_alpha, ngram_size)
    // pairs can be dispatched to a copy with the alphabet base and top digit as COMPILE-TIME
    // constants: the compiler then strength-reduces the base multiplies into shift+add and
    // fully unrolls the fixed-size init loop, shortening the loop-carried integer chain. The
    // specialised pairs are the ones the scoring-bound families actually use -- the J-merged
    // squares (Playfair/Bifid/ADFGX/... all g_alpha == 25), the plain 26-letter types, the
    // 27-symbol cube/grid types (Trifid/Digrafid) and the 36-symbol ADFGVX -- at quadgram
    // and quintgram sizes. The integer index is identical whether the multiply is an imul or
    // a lea chain, so it is bit-for-bit identical to the generic path, which still handles any
    // other alphabet size or n-gram order. (g_alpha==36 quintgram is intentionally NOT
    // specialised: a base-36 quintgram table is 60M entries / 240MB and the codebase avoids
    // it.) NOTE: once the compact-shadow table is armed (large tables only, see below), the
    // 27/36 pairs whose float table exceeds the shadow threshold take the shadow path instead
    // -- these float specialisations then serve the small-table (27^4) case and the
    // shadow-disabled fallback.
    // Compressed 8-bit table (see load_ngrams/.ngbin): score through g_ngram_lut with
    // the SAME (alpha, size) specialization scheme as the float path. The branch is
    // hoisted out of the per-window loop (checked once per ngram_score call, and always
    // NULL/perfectly-predicted for existing solves), so the else branch runs the float
    // dispatch with byte-for-byte the historical arithmetic -- every existing solve is
    // numerically unchanged and the regression suite stays bit-identical. The dense byte
    // table is practical through order 6 (309 MB at alpha=26), so 26/6 & 25/6 are
    // specialized here even though the float path (which never holds a 6-gram table by
    // policy) does not.
    if (g_ngram_u8) {
        const unsigned char * restrict nd8 = g_ngram_u8;
        const float * restrict lut = g_ngram_lut;
        if      (g_alpha == 25 && ngram_size == 4) score = ngram_walk_u8(src, m, nd8, lut, 4, 25,     15625);
        else if (g_alpha == 26 && ngram_size == 4) score = ngram_walk_u8(src, m, nd8, lut, 4, 26,     17576);
        else if (g_alpha == 25 && ngram_size == 5) score = ngram_walk_u8(src, m, nd8, lut, 5, 25,    390625);
        else if (g_alpha == 26 && ngram_size == 5) score = ngram_walk_u8(src, m, nd8, lut, 5, 26,    456976);
        else if (g_alpha == 25 && ngram_size == 6) score = ngram_walk_u8(src, m, nd8, lut, 6, 25,   9765625);
        else if (g_alpha == 26 && ngram_size == 6) score = ngram_walk_u8(src, m, nd8, lut, 6, 26,  11881376);
        else if (g_alpha == 27 && ngram_size == 4) score = ngram_walk_u8(src, m, nd8, lut, 4, 27,     19683);
        else if (g_alpha == 27 && ngram_size == 5) score = ngram_walk_u8(src, m, nd8, lut, 5, 27,    531441);
        else                                       score = ngram_walk_u8(src, m, nd8, lut, ngram_size, g_alpha, top);
    } else {
        const float * restrict nd = ngram_data;
        if      (g_alpha == 25 && ngram_size == 4) score = ngram_walk(src, m, nd, 4, 25,  15625);
        else if (g_alpha == 26 && ngram_size == 4) score = ngram_walk(src, m, nd, 4, 26,  17576);
        else if (g_alpha == 25 && ngram_size == 5) score = ngram_walk(src, m, nd, 5, 25, 390625);
        else if (g_alpha == 26 && ngram_size == 5) score = ngram_walk(src, m, nd, 5, 26, 456976);
        else if (g_alpha == 27 && ngram_size == 4) score = ngram_walk(src, m, nd, 4, 27,  19683);
        else if (g_alpha == 27 && ngram_size == 5) score = ngram_walk(src, m, nd, 5, 27, 531441);
        else if (g_alpha == 36 && ngram_size == 4) score = ngram_walk(src, m, nd, 4, 36,  46656);
        else                                       score = ngram_walk(src, m, nd, ngram_size, g_alpha, top);
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

    // When a compressed 8-bit table is active, gather through the dequantizing LUT
    // (the float pointer is then NULL, load_ngrams having returned NULL for a .ngbin).
    // This is off the hot solve path (seam decomposition only), so a per-gather select
    // is fine; nd8 == NULL leaves the original float gather byte-for-byte unchanged.
    const unsigned char *nd8 = g_ngram_u8;

    // Big-endian packing, same convention (and same table) as ngram_score: the roll
    // is a multiply, and the outgoing letter sits at the top digit (weight `top`). A
    // sentinel contributed 0 to the index, so subtracting out_v*top (out_v forced to
    // 0 for a sentinel) removes exactly what was added -- bit-identical to before.
    int index = 0, bad = 0;
    for (int j = 0; j < ngram_size; j++) {
        int v = text[j];
        if (v < 0) { bad++; v = 0; }       // sentinel: contributes 0 to the packed index
        index = index * g_alpha + v;
    }
    if (bad == 0) score += nd8 ? (double) g_ngram_lut[nd8[index]] : ngram_data[index];

    for (int i = 1; i < n_windows; i++) {
        int out_v = text[i - 1];
        int in_v  = text[i + ngram_size - 1];
        if (out_v < 0) { bad--; out_v = 0; }
        int in_iv = in_v;
        if (in_v < 0) { bad++; in_iv = 0; }
        index = (index - out_v * top) * g_alpha + in_iv;
        if (bad == 0) score += nd8 ? (double) g_ngram_lut[nd8[index]] : ngram_data[index];
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

// Load a dense 8-bit .ngbin table (see write_ngram_bin / scoring.h). mmaps the payload
// into g_ngram_u8, fills the byte->weight LUT, and forces log-prob semantics. Returns
// NULL -- the caller stores it as the (now unused) float table and ngram_score /
// ngram_sum_raw read g_ngram_u8 instead. Hard-errors on any header/runtime mismatch.
static float* load_ngrams_bin(const char *path, int ngram_size, bool verbose) {
    if (g_ngram_reverse) {
        fprintf(stderr, "Error: -reversengrams is not supported with a compressed .ngbin table "
                        "(read-only mmap).\n       Use the text table or pre-bake a reversed .ngbin.\n");
        exit(1);
    }

    int fd = open(path, O_RDONLY);
    if (fd < 0) { fprintf(stderr, "Error: cannot open ngram file '%s'\n", path); exit(1); }

    struct stat st;
    if (fstat(fd, &st) != 0 || (size_t) st.st_size < (size_t) NGBIN_HEADER_SIZE) {
        fprintf(stderr, "Error: '%s' is too small to be a .ngbin table\n", path);
        exit(1);
    }
    size_t flen = (size_t) st.st_size;

    NgramBinHeader h;
    if (pread(fd, &h, sizeof h, 0) != (ssize_t) sizeof h) {
        fprintf(stderr, "Error: cannot read .ngbin header from '%s'\n", path);
        exit(1);
    }

    uint64_t expect = 1;
    for (int i = 0; i < ngram_size; i++) expect *= (uint64_t) g_alpha;

    if (memcmp(h.magic, NGBIN_MAGIC, NGBIN_MAGIC_LEN) != 0 || h.version != NGBIN_VERSION) {
        fprintf(stderr, "Error: '%s' is not a v%d .ngbin table\n", path, NGBIN_VERSION);
        exit(1);
    }
    if (h.mode != NGBIN_MODE_LOGPROB) {
        fprintf(stderr, "Error: .ngbin '%s' has unsupported weighting mode %u\n", path, h.mode);
        exit(1);
    }
    if (h.order != ngram_size) {
        fprintf(stderr, "Error: .ngbin '%s' is order %u but -ngramsize is %d\n",
                path, h.order, ngram_size);
        exit(1);
    }
    if (h.alphabet_size != g_alpha) {
        fprintf(stderr, "Error: .ngbin '%s' was built for a %u-letter alphabet but the runtime "
                        "alphabet is %d.\n       A compressed table is alphabet-specific -- rebuild it "
                        "for this cipher type (see -writengrambin).\n", path, h.alphabet_size, g_alpha);
        exit(1);
    }
    if (h.n_entries != expect || expect > (uint64_t) INT_MAX) {
        fprintf(stderr, "Error: .ngbin '%s' entry count %llu != expected %llu (or exceeds the int "
                        "index range)\n", path, (unsigned long long) h.n_entries,
                (unsigned long long) expect);
        exit(1);
    }
    if (flen < (size_t) NGBIN_HEADER_SIZE + h.n_entries) {
        fprintf(stderr, "Error: .ngbin '%s' payload is truncated\n", path);
        exit(1);
    }

    // mmap the whole file; the payload begins right after the 64-byte header.
    void *base = mmap(NULL, flen, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    const unsigned char *payload;
    if (base == MAP_FAILED) {
        // Fall back to a plain read (still 1 byte/entry). Leaked at exit, like the
        // never-freed float table -- the process is short-lived.
        unsigned char *buf = malloc(flen);
        int fd2 = open(path, O_RDONLY);
        if (!buf || fd2 < 0 || pread(fd2, buf, flen, 0) != (ssize_t) flen) {
            fprintf(stderr, "Error: cannot read .ngbin payload from '%s'\n", path);
            exit(1);
        }
        close(fd2);
        payload = buf + NGBIN_HEADER_SIZE;
        g_ngram_mmap_len = 0;
    } else {
        payload = (const unsigned char *) base + NGBIN_HEADER_SIZE;
        g_ngram_mmap_len = flen;
    }

    for (int b = 0; b < 256; b++)
        g_ngram_lut[b] = (float) (h.w_floor + (double) b * h.w_scale);
    g_ngram_u8 = payload;
    g_ngram_floor = h.w_floor;      // entropy-term base (matches the float logprob path)
    g_ngram_logprob = true;         // the format stores log10 probabilities

    if (verbose)
        printf("\nLoaded compressed ngrams: order %d, alphabet %d, %llu entries (%.1f MB, %s).\n\n",
               ngram_size, g_alpha, (unsigned long long) h.n_entries,
               h.n_entries / 1048576.0, (base == MAP_FAILED) ? "read" : "mmap");
    return NULL;
}

float* load_ngrams(char *ngram_file, int ngram_size, bool verbose) {
    FILE *fp;
    int i, n_ngrams, freq, indx;
    char ngram[MAX_NGRAM_SIZE];
    float *ngram_data, total;

    // A dense 8-bit .ngbin (magic "COLNGBIN") is mmap'd and scored via g_ngram_lut; any
    // other file is the historical space-separated text table parsed below.
    {
        FILE *pf = fopen(ngram_file, "rb");
        if (pf) {
            char magic[NGBIN_MAGIC_LEN];
            size_t got = fread(magic, 1, NGBIN_MAGIC_LEN, pf);
            fclose(pf);
            if (got == (size_t) NGBIN_MAGIC_LEN && memcmp(magic, NGBIN_MAGIC, NGBIN_MAGIC_LEN) == 0)
                return load_ngrams_bin(ngram_file, ngram_size, verbose);
        }
    }

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
        g_ngram_floor = floor;   // min per-window value: shift base for the entropy term
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
        // it symmetric and order-independent. Independent of endianness: this
        // extracts each digit low-to-high (% g_alpha) and re-places it high-to-low
        // (weights topw..1), which maps a slot to its digit-reversed twin's slot
        // under either packing -- so it is unchanged by the big-endian switch.
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

int write_ngram_bin(const char *path, const float *ngram_data, int ngram_size) {
    if (!ngram_data || g_ngram_u8) {
        fprintf(stderr, "Error: -writengrambin needs a TEXT ngram table as input (the loaded "
                        "table is already compressed).\n");
        return 1;
    }
    if (!g_ngram_logprob) {
        fprintf(stderr, "Error: -writengrambin requires -logprob (the .ngbin format stores "
                        "log10 probabilities).\n");
        return 1;
    }
    uint64_t n = 1;
    for (int i = 0; i < ngram_size; i++) n *= (uint64_t) g_alpha;
    if (n > (uint64_t) INT_MAX) {
        fprintf(stderr, "Error: a %d-gram table over a %d-letter alphabet has %llu entries, too "
                        "large for the dense .ngbin format (int index limit).\n",
                ngram_size, g_alpha, (unsigned long long) n);
        return 1;
    }

    // The in-memory log-prob table sets unseen cells to g_ngram_floor (the minimum), so
    // quantize [floor, max] into 0..255: unseen -> 0, dequantizing back to floor exactly.
    double w_floor = g_ngram_floor;
    double w_max = ngram_data[0];
    for (uint64_t i = 1; i < n; i++) if (ngram_data[i] > w_max) w_max = ngram_data[i];
    if (w_max <= w_floor) w_max = w_floor + 1.0;          // degenerate guard (empty table)
    double w_scale = (w_max - w_floor) / 255.0;

    NgramBinHeader h;
    memset(&h, 0, sizeof h);
    memcpy(h.magic, NGBIN_MAGIC, NGBIN_MAGIC_LEN);
    h.version = NGBIN_VERSION;
    h.order = (uint8_t) ngram_size;
    h.alphabet_size = (uint8_t) g_alpha;
    h.mode = NGBIN_MODE_LOGPROB;
    h.n_entries = n;
    h.w_floor = w_floor;
    h.w_scale = w_scale;
    h.total_count = 0.01 * pow(10.0, -w_floor);          // reference: floor = log10(0.01/total)

    FILE *fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "Error: cannot open '%s' for writing\n", path); return 1; }
    if (fwrite(&h, 1, sizeof h, fp) != sizeof h) {
        fprintf(stderr, "Error: .ngbin header write failed\n"); fclose(fp); return 1;
    }

    unsigned char *buf = malloc((size_t) n);
    if (!buf) { fprintf(stderr, "Error: out of memory writing .ngbin\n"); fclose(fp); return 1; }
    double inv = 1.0 / w_scale, max_abs_err = 0.0;
    for (uint64_t i = 0; i < n; i++) {
        long b = lround((ngram_data[i] - w_floor) * inv);
        if (b < 0) b = 0; else if (b > 255) b = 255;
        buf[i] = (unsigned char) b;
        double err = fabs((w_floor + (double) b * w_scale) - ngram_data[i]);
        if (err > max_abs_err) max_abs_err = err;
    }
    size_t wrote = fwrite(buf, 1, (size_t) n, fp);
    free(buf);
    fclose(fp);
    if (wrote != (size_t) n) { fprintf(stderr, "Error: .ngbin payload write failed\n"); return 1; }

    printf("Wrote %s: order %d, alphabet %d, %llu entries (%.2f MB), floor=%.4f scale=%.6g, "
           "max quant err=%.4f.\n", path, ngram_size, g_alpha, (unsigned long long) n,
           (NGBIN_HEADER_SIZE + n) / 1048576.0, w_floor, w_scale, max_abs_err);
    return 0;
}

int ngram_index_str(char *ngram, int ngram_size) {
    // Big-endian packing (idx = sum c*g_alpha^(n-1-j)), matching ngram_score's
    // window walk exactly -- the two conventions must agree or the table is misread.
    int c, index = 0;
    for (int i = 0; i < ngram_size; i++) {
        c = g_char_to_idx[toupper((unsigned char) ngram[i]) & 127];
        // An n-gram containing a letter outside the runtime alphabet (e.g. 'P'
        // under -excludeletter P) cannot occur in the plaintext, so it has no
        // slot; signal the caller to skip it.
        if (c < 0) return -1;
        index = index * g_alpha + c;
    }
    return index;
}

int ngram_index_int(int *ngram, int ngram_size) {
    int index = 0;                         // big-endian, matching ngram_index_str
    for (int i = 0; i < ngram_size; i++)
        index = index * g_alpha + ngram[i];
    return index;
}
