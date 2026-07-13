// -spaces support: see spaces.h. The n-gram file is a plain CHARACTER n-gram table over an
// EXTENDED 27-symbol alphabet {A..Z, ' '} -- not a word-vocabulary model. Each line is an
// `order`-character window (letters and/or embedded literal spaces) followed by its corpus
// frequency, e.g. "AAAAA 19157" (order 5, no space in the window) or "AAA A 1264" (order 5,
// the window is "AAA" + a space + "A"). Because space is just another symbol scored by the
// SAME rolling n-gram model, finding the best place(s) to insert spaces into an unbroken
// letter run is an exact Viterbi decode: maximize the total log-probability of the emitted
// {A..Z,' '} sequence under this Markov model, chain rule over one symbol at a time. This is
// the reason a dense flat array (matching load_ngrams()'s convention, base 27 instead of
// g_alpha) is used rather than a hash map -- the alphabet is tiny and fixed, so nothing here
// is sparse the way a word-vocabulary model would be.

#include "colossus.h"
#include "spaces.h"

SpacesNgramTable *g_spaces_table = NULL;

#define SPACES_ALPHA 27
#define SPACES_SPACE_SYM 26
#define SPACES_MAX_ORDER 8              // matches MAX_NGRAM_SIZE; higher orders are impractical
                                          // to allocate densely (27^9 alone is ~7.6e12 cells)
#define SPACES_MAX_LINE 128

struct SpacesNgramTable {
    float *data;   // log10 P(window), size SPACES_ALPHA^order, indexed base-27 positional
    int order;
};

static long long sp_pow(int base, int exp) {
    long long r = 1;
    for (int i = 0; i < exp; i++) r *= base;
    return r;
}

static int sp_char_to_sym(char c) {
    if (c == ' ') return SPACES_SPACE_SYM;
    unsigned char u = (unsigned char) toupper((unsigned char) c);
    if (u >= 'A' && u <= 'Z') return u - 'A';
    return -1;
}

SpacesNgramTable *load_spaces_ngrams(const char *filename, int order, bool verbose) {
    if (order < 1 || order > SPACES_MAX_ORDER) {
        printf("ERROR: -spacesngramsize %d out of range (1..%d)\n", order, SPACES_MAX_ORDER);
        return NULL;
    }
    long long n_entries = sp_pow(SPACES_ALPHA, order);
    if (n_entries > (long long) 600000000) {   // ~2.4GB of float -- order 7+ rejected
        printf("ERROR: -spacesngramsize %d needs a %lld-cell table, too large to allocate\n",
               order, n_entries);
        return NULL;
    }
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("ERROR: could not open -spacesngramfile %s\n", filename);
        return NULL;
    }

    double *counts = calloc((size_t) n_entries, sizeof(double));
    double total = 0.0;
    char line[SPACES_MAX_LINE];

    while (fgets(line, sizeof(line), fp)) {
        int len = (int) strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;
        // The rightmost space separates the frequency (always plain digits) from the
        // n-gram window, which may itself contain embedded single spaces.
        char *last_space = strrchr(line, ' ');
        if (!last_space || last_space[1] == '\0') continue;
        long freq = atol(last_space + 1);
        if (freq <= 0) continue;
        int ngram_len = (int) (last_space - line);
        if (ngram_len != order) continue;
        long long idx = 0;
        bool ok = true;
        for (int k = 0; k < order; k++) {
            int sym = sp_char_to_sym(line[k]);
            if (sym < 0) { ok = false; break; }
            idx = idx * SPACES_ALPHA + sym;
        }
        if (!ok) continue;
        counts[idx] += (double) freq;
        total += (double) freq;
    }
    fclose(fp);

    if (total <= 0.0) {
        printf("ERROR: -spacesngramfile %s had no usable order-%d entries\n", filename, order);
        free(counts);
        return NULL;
    }

    SpacesNgramTable *tbl = malloc(sizeof(SpacesNgramTable));
    tbl->order = order;
    tbl->data = malloc((size_t) n_entries * sizeof(float));
    // Unseen windows must be penalised, not left at 0: most candidate segmentations probe
    // windows the corpus never saw, and a flat 0 gives the Viterbi decode no way to prefer a
    // real segmentation over a nonsense one. Same floor convention as load_ngrams()'s -logprob
    // branch (scoring.c): a fixed small pseudo-count relative to the corpus total.
    double floor_val = log10(0.01 / total);
    for (long long i = 0; i < n_entries; i++) {
        tbl->data[i] = (counts[i] > 0.0) ? (float) log10(counts[i] / total) : (float) floor_val;
    }
    free(counts);
    if (verbose) {
        printf("Loaded -spacesngramfile %s: order %d, %.0f total count.\n", filename, order, total);
    }
    return tbl;
}

void free_spaces_ngrams(SpacesNgramTable *tbl) {
    if (!tbl) return;
    free(tbl->data);
    free(tbl);
}

// One Viterbi state: the trailing (order-1) EMITTED symbols (letters or SPACES_SPACE_SYM),
// only the first ctx_len of which are valid (fewer at the very start of a run -- those
// windows are simply not scored yet, matching ngram_score()'s "only full windows count").
typedef struct {
    int ctx[SPACES_MAX_ORDER - 1];
    int ctx_len;
    double score;
    int prev;             // pool index of the predecessor state, -1 for the root
    int emitted_space;    // 1 if the edge into this state inserted a space, 0 if it consumed
                           // the next input letter
} SpState;

static int sp_ctx_eq(const SpState *a, const SpState *b) {
    if (a->ctx_len != b->ctx_len) return 0;
    for (int k = 0; k < a->ctx_len; k++) if (a->ctx[k] != b->ctx[k]) return 0;
    return 1;
}

static void sp_append(int ctx[], int *ctx_len, int cap, int sym) {
    if (cap == 0) return;                  // order == 1: no context is ever needed
    if (*ctx_len < cap) {
        ctx[(*ctx_len)++] = sym;
    } else {
        for (int k = 1; k < cap; k++) ctx[k-1] = ctx[k];
        ctx[cap-1] = sym;
    }
}

static double sp_window_score(const SpacesNgramTable *tbl, const int ctx[], int ctx_len, int sym) {
    int need = tbl->order - 1;
    if (ctx_len < need) return 0.0;        // not enough history yet -- free, uncounted
    long long idx = 0;
    for (int k = 0; k < need; k++) idx = idx * SPACES_ALPHA + ctx[k];
    idx = idx * SPACES_ALPHA + sym;
    return (double) tbl->data[idx];
}

// Exact Viterbi segmentation of a single pure-letter run (run_len >= 0 alphabet-index
// symbols already validated as plain A..Z). Reachable contexts per position are bounded by
// 2^(order-1) (each of the last order-1 letter boundaries either got a space or didn't --
// consecutive spaces are disallowed, so nothing larger is reachable), so this stays fast
// regardless of run length.
static char *spaces_segment_run(const SpacesNgramTable *tbl, const int sym[], int run_len) {
    if (run_len == 0) { char *e = malloc(1); e[0] = '\0'; return e; }

    int pool_cap = 1024, pool_n = 0;
    SpState *pool = malloc((size_t) pool_cap * sizeof(SpState));
#define SP_GROW_POOL() do { if (pool_n == pool_cap) { pool_cap *= 2; pool = realloc(pool, (size_t) pool_cap * sizeof(SpState)); } } while (0)

    int level_cap = 256;
    int *level = malloc((size_t) level_cap * sizeof(int));
    int *next_level = malloc((size_t) level_cap * sizeof(int));
    int level_n = 0, next_n = 0;

    SpState root = { .ctx_len = 0, .score = 0.0, .prev = -1, .emitted_space = 0 };
    SP_GROW_POOL(); pool[pool_n] = root; level[level_n++] = pool_n++;

    for (int i = 0; i < run_len; i++) {
        // Phase 1: from every state currently at position i, optionally insert a space
        // (skip if the trailing symbol is already a space -- no double spaces -- or if this
        // would be a leading space). Space-insertion states are appended into the SAME level
        // (still at position i); a space state can never itself take another Phase-1 step
        // (its trailing symbol is now SPACES_SPACE_SYM), so this closes in one pass.
        int raw_n = level_n;
        for (int s = 0; s < raw_n; s++) {
            SpState *st = &pool[level[s]];
            int last = (st->ctx_len > 0) ? st->ctx[st->ctx_len - 1] : -1;
            if (last == SPACES_SPACE_SYM || i == 0) continue;
            SpState ns = *st;
            ns.score += sp_window_score(tbl, st->ctx, st->ctx_len, SPACES_SPACE_SYM);
            sp_append(ns.ctx, &ns.ctx_len, tbl->order - 1, SPACES_SPACE_SYM);
            ns.prev = level[s];
            ns.emitted_space = 1;
            int dup = -1;
            for (int t = 0; t < level_n; t++) if (sp_ctx_eq(&pool[level[t]], &ns)) { dup = t; break; }
            if (dup >= 0) {
                if (ns.score > pool[level[dup]].score) pool[level[dup]] = ns;
            } else {
                if (level_n == level_cap) {
                    level_cap *= 2;
                    level = realloc(level, (size_t) level_cap * sizeof(int));
                    next_level = realloc(next_level, (size_t) level_cap * sizeof(int));
                }
                SP_GROW_POOL(); pool[pool_n] = ns; level[level_n++] = pool_n++;
            }
        }

        // Phase 2: consume input letter i from every state at position i, producing i+1.
        next_n = 0;
        for (int s = 0; s < level_n; s++) {
            SpState *st = &pool[level[s]];
            SpState ns = *st;
            ns.score += sp_window_score(tbl, st->ctx, st->ctx_len, sym[i]);
            sp_append(ns.ctx, &ns.ctx_len, tbl->order - 1, sym[i]);
            ns.prev = level[s];
            ns.emitted_space = 0;
            int dup = -1;
            for (int t = 0; t < next_n; t++) if (sp_ctx_eq(&pool[next_level[t]], &ns)) { dup = t; break; }
            if (dup >= 0) {
                if (ns.score > pool[next_level[dup]].score) pool[next_level[dup]] = ns;
            } else {
                if (next_n == level_cap) {
                    level_cap *= 2;
                    level = realloc(level, (size_t) level_cap * sizeof(int));
                    next_level = realloc(next_level, (size_t) level_cap * sizeof(int));
                }
                SP_GROW_POOL(); pool[pool_n] = ns; next_level[next_n++] = pool_n++;
            }
        }

        int *tmp = level; level = next_level; next_level = tmp;
        level_n = next_n;
    }

    int best = level[0];
    for (int s = 1; s < level_n; s++) if (pool[level[s]].score > pool[best].score) best = level[s];

    // Reconstruct by walking backpointers from `best` to the root. Every non-root state
    // emitted exactly one symbol: a space, or the next input letter in REVERSE consumption
    // order (letter-consuming edges walk sym[] backward from run_len-1 down to 0).
    char *rev = malloc((size_t) (2 * run_len + 1));
    int rn = 0, letter_pos = run_len - 1;
    int cur = best;
    while (pool[cur].prev != -1) {
        if (pool[cur].emitted_space) {
            rev[rn++] = ' ';
        } else {
            rev[rn++] = (char) ('A' + sym[letter_pos]);
            letter_pos--;
        }
        cur = pool[cur].prev;
    }
    char *out = malloc((size_t) rn + 1);
    for (int k = 0; k < rn; k++) out[k] = rev[rn - 1 - k];
    out[rn] = '\0';
    free(rev);

    free(pool);
    free(level);
    free(next_level);
    return out;
}

char *spaces_insert(const SpacesNgramTable *tbl, int indices[], int len) {
    if (!tbl) return NULL;

    char *out_buf = malloc((size_t) (2 * len + 2));
    int out_len = 0;
    int run_start = -1;

    for (int i = 0; i <= len; i++) {
        int is_live = (i < len) && (indices[i] >= 0);
        if (is_live && run_start < 0) run_start = i;
        if ((!is_live || i == len) && run_start >= 0) {
            int run_len = i - run_start;
            // Convert this run to A..Z symbols; if it contains anything other than a plain
            // letter (e.g. ADFGVX's digits, Trifid's '+'), segmenting it isn't meaningful --
            // reproduce it verbatim instead of guessing.
            int *sym = malloc((size_t) run_len * sizeof(int));
            bool plain = true;
            for (int k = 0; k < run_len; k++) {
                char c = (char) index_to_char(indices[run_start + k]);
                int s = sp_char_to_sym(c);
                if (s < 0 || s == SPACES_SPACE_SYM) { plain = false; break; }
                sym[k] = s;
            }
            if (plain) {
                char *seg = spaces_segment_run(tbl, sym, run_len);
                int seg_len = (int) strlen(seg);
                memcpy(out_buf + out_len, seg, (size_t) seg_len);
                out_len += seg_len;
                free(seg);
            } else {
                for (int k = 0; k < run_len; k++) out_buf[out_len++] = (char) index_to_char(indices[run_start + k]);
            }
            free(sym);
            run_start = -1;
        }
        if (i < len && !is_live) out_buf[out_len++] = (char) index_to_char(indices[i]);
    }
    out_buf[out_len] = '\0';
    return out_buf;
}

void print_spaces_line(const SpacesNgramTable *tbl, int indices[], int len) {
    if (!tbl) return;
    char *spaced = spaces_insert(tbl, indices, len);
    if (!spaced) return;
    printf("with spaces: %s\n", spaced);
    free(spaced);
}
