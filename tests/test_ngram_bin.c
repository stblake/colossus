// Unit tests for the dense 8-bit .ngbin compressed n-gram table (see
// src/core/scoring.c: write_ngram_bin / load_ngrams / load_ngrams_bin, and the
// g_ngram_u8 / g_ngram_lut hot path in ngram_score). We:
//   1. build a small text table, load it as the float -logprob table;
//   2. dump it to a .ngbin with write_ngram_bin;
//   3. reload the .ngbin (which mmaps g_ngram_u8 + fills g_ngram_lut) and assert every
//      dequantized cell is within one quantization step of the original float weight;
//   4. assert ngram_score via the byte table tracks ngram_score via the float table
//      (within one step), and that the g_ngram_u8 == NULL branch still runs the float
//      table unchanged.
// The header-mismatch error paths (wrong alphabet / order / magic) exit() the process,
// so they are covered by the CLI smoke test, not here.
//
// Build (see makefile `test` target):
//   gcc -I... tests/test_ngram_bin.c src/core/scoring.c src/core/utils.c -o tests/test_ngram_bin

#include "colossus.h"
#include "scoring.h"
#include "spaces.h"

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { printf("  FAIL: %s\n", msg); failures++; } } while (0)

#define TXT_PATH "/tmp/colossus_test_ngbin.txt"
#define BIN_PATH "/tmp/colossus_test_ngbin.ngbin"
#define SBIN_PATH "/tmp/colossus_test_spaces.ngbin"

// A small trigram table with varied counts; the rest of the 26^3 grid stays unseen.
static void write_text_table(void) {
    FILE *fp = fopen(TXT_PATH, "w");
    if (!fp) { printf("  FAIL: cannot open %s\n", TXT_PATH); failures++; return; }
    const char *ng[] = {"THE","AND","ING","HER","ERE","ENT","THA","NTH","WAS","EDT",
                        "TIS","OFT","STH","MEN","XQZ","JKV"};
    const int   ct[] = { 90000, 55000, 41000, 33000, 25000, 21000, 18000, 15000,
                          9000,  7000,  5000,  3000,  1200,   400,     3,     1};
    for (int i = 0; i < 16; i++) fprintf(fp, "%s %d\n", ng[i], ct[i]);
    fclose(fp);
}

static void test_ngbin_roundtrip(void) {
    printf("ngbin round-trip fidelity + score parity...\n");
    CHECK(sizeof(NgramBinHeader) == NGBIN_HEADER_SIZE, "NgramBinHeader is not 64 bytes");

    const int ng_size = 3;
    init_alphabet(NULL);          // populate g_char_to_idx (the real binary does this at startup)
    write_text_table();

    // 1. Float -logprob table.
    g_ngram_logprob = true;
    g_ngram_reverse = false;
    g_ngram_u8 = NULL;
    char txt[] = TXT_PATH;
    float *ft = load_ngrams(txt, ng_size, false);
    CHECK(ft != NULL, "float table load returned NULL");
    CHECK(g_ngram_u8 == NULL, "text load must not set the u8 table");

    // 2. Dump to .ngbin.
    int rc = write_ngram_bin(BIN_PATH, ft, ng_size);
    CHECK(rc == 0, "write_ngram_bin failed");

    // 3. Reload the .ngbin: sets g_ngram_u8 + g_ngram_lut, returns NULL.
    char bin[] = BIN_PATH;
    float *bt = load_ngrams(bin, ng_size, false);
    CHECK(bt == NULL, "compressed load should return NULL (u8 table is used instead)");
    CHECK(g_ngram_u8 != NULL, "compressed load must set g_ngram_u8");

    // Dequantization step recovered from the LUT: lut[b] = floor + b*scale.
    double lut_floor = g_ngram_lut[0];
    double scale = (double) g_ngram_lut[1] - (double) g_ngram_lut[0];
    CHECK(scale > 0.0, "LUT scale must be positive");
    CHECK(fabs(lut_floor - g_ngram_floor) < 1e-4, "g_ngram_floor should match the LUT floor");
    double cell_tol = 0.51 * scale + 1e-3;

    // Every cell within one quantization step of the float weight.
    int n_entries = int_pow(g_alpha, ng_size);
    double max_cell_err = 0.0;
    for (int i = 0; i < n_entries; i++) {
        double deq = g_ngram_lut[g_ngram_u8[i]];
        double err = fabs(deq - (double) ft[i]);
        if (err > max_cell_err) max_cell_err = err;
    }
    printf("  scale=%.5f  max cell err=%.5f  (tol=%.5f)\n", scale, max_cell_err, cell_tol);
    CHECK(max_cell_err <= cell_tol, "a dequantized cell exceeds one quantization step");

    // 4. Score parity on a sample stream: byte-table score tracks float-table score.
    seed_rand(20240607);
    int len = 120;
    int pt[120];
    for (int i = 0; i < len; i++) pt[i] = rand_int(0, g_alpha);
    g_score_no_sentinel = true;                 // all-letters stream

    double score_u8 = ngram_score(pt, len, NULL, ng_size);   // uses g_ngram_u8
    const unsigned char *saved = g_ngram_u8;
    g_ngram_u8 = NULL;                                        // force the float path
    double score_ft = ngram_score(pt, len, ft, ng_size);
    g_ngram_u8 = saved;

    printf("  score(byte)=%.5f  score(float)=%.5f\n", score_u8, score_ft);
    CHECK(fabs(score_u8 - score_ft) <= 0.6 * scale + 1e-3,
          "byte-table score diverges from float-table score by more than one step");
    // Sanity: the float path with g_ngram_u8==NULL is the genuine historical scorer.
    CHECK(isfinite(score_ft), "float-table score is not finite");

    remove(TXT_PATH);
    remove(BIN_PATH);
    printf("  ok\n");
}

// Write a minimal .ngbin (every cell at the floor byte 0) with an explicit alphabet_size,
// for exercising the spaces reader in src/core/spaces.c.
static void write_ngbin_all_floor(const char *path, int order, int alpha,
                                  double w_floor, double w_scale) {
    long long n = 1;
    for (int i = 0; i < order; i++) n *= alpha;
    NgramBinHeader h;
    memset(&h, 0, sizeof h);
    memcpy(h.magic, NGBIN_MAGIC, NGBIN_MAGIC_LEN);
    h.version = NGBIN_VERSION;
    h.order = (uint8_t) order;
    h.alphabet_size = (uint8_t) alpha;
    h.mode = NGBIN_MODE_LOGPROB;
    h.n_entries = (uint64_t) n;
    h.w_floor = w_floor;
    h.w_scale = w_scale;
    FILE *fp = fopen(path, "wb");
    fwrite(&h, 1, sizeof h, fp);
    unsigned char *buf = calloc((size_t) n, 1);     // all 0 => floor => "unseen"
    fwrite(buf, 1, (size_t) n, fp);
    free(buf);
    fclose(fp);
}

static void test_spaces_ngbin(void) {
    printf("spaces .ngbin reader (27-symbol) + mismatch rejection...\n");

    // Valid order-2, 27-symbol spaces table (all cells at the floor => each inserted space
    // only costs the floor penalty, so the Viterbi inserts none).
    write_ngbin_all_floor(SBIN_PATH, 2, 27, -6.0, 6.0 / 255.0);
    SpacesNgramTable *tbl = load_spaces_ngrams(SBIN_PATH, 2, false);
    CHECK(tbl != NULL, "load_spaces_ngrams should accept a valid 27-symbol .ngbin");
    if (tbl) {
        int idx[3] = {19, 7, 4};                    // "THE"
        char *s = spaces_insert(tbl, idx, 3);
        CHECK(s != NULL, "spaces_insert returned NULL");
        if (s) {
            CHECK(strcmp(s, "THE") == 0, "an all-floor table should segment THE with no spaces");
            free(s);
        }
        free_spaces_ngrams(tbl);
    }

    // Rejection: a 26-symbol (letter) header must be refused by the spaces reader
    // (returns NULL rather than exit()ing -- the expected ERROR line below is intentional).
    write_ngbin_all_floor(SBIN_PATH, 2, 26, -6.0, 6.0 / 255.0);
    SpacesNgramTable *bad = load_spaces_ngrams(SBIN_PATH, 2, false);
    CHECK(bad == NULL, "spaces reader must reject a 26-symbol .ngbin (alphabet mismatch)");

    remove(SBIN_PATH);
    printf("  ok\n");
}

int main(void) {
    printf("\n=== test_ngram_bin ===\n");
    init_alphabet(NULL);          // populate g_char_to_idx / g_idx_to_char_arr for both tests
    test_ngbin_roundtrip();
    test_spaces_ngbin();
    if (failures) { printf("\n%d CHECK(s) FAILED\n", failures); return 1; }
    printf("\nAll .ngbin tests passed.\n");
    return 0;
}
