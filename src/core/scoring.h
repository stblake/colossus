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

// ---- Compressed dense 8-bit n-gram table (.ngbin) --------------------------
// A dense byte-per-n-gram image of the in-memory table: the byte at big-endian
// index i (AAA=0, AAB=1, ... , exactly ngram_index_str's packing) is the 8-bit
// quantization of that n-gram's log10 weight, w = w_floor + byte*w_scale. The
// reader (load_ngrams, dispatched by the magic) mmaps the payload into g_ngram_u8
// and ngram_score dequantizes via g_ngram_lut. Only the log-prob weighting is
// stored (mode 0) and it is ALPHABET-SPECIFIC: a table built for g_alpha=26 is a
// different image from the 25-letter (J->I) one a Bifid solve needs, so the reader
// validates alphabet_size == g_alpha. Dense => practical only through order 6
// (26^6 = 309 MB; 26^7 = 8 GB and would overflow the int walk index). Written by
// write_ngram_bin (the `-writengrambin` dump reuses the exact per-type table).
#define NGBIN_MAGIC        "COLNGBIN"   // 8 bytes, compared without a NUL
#define NGBIN_MAGIC_LEN    8
#define NGBIN_VERSION      1
#define NGBIN_MODE_LOGPROB 0
#define NGBIN_HEADER_SIZE  64

#pragma pack(push, 1)
typedef struct {
    char     magic[8];        // "COLNGBIN"
    uint16_t version;         // NGBIN_VERSION
    uint8_t  order;           // n-gram order n (must equal -ngramsize)
    uint8_t  alphabet_size;   // g_alpha the table was built for
    uint8_t  mode;            // NGBIN_MODE_LOGPROB
    uint8_t  reserved0[3];
    uint64_t n_entries;       // alphabet_size ^ order
    double   w_floor;         // weight of byte 0 (unseen) = log10(0.01/total)
    double   w_scale;         // weight = w_floor + byte*w_scale
    double   total_count;     // reference only (derived from the floor)
    uint8_t  reserved1[16];
} NgramBinHeader;             // exactly 64 bytes, little-endian on disk
#pragma pack(pop)
_Static_assert(sizeof(NgramBinHeader) == NGBIN_HEADER_SIZE, "NgramBinHeader must be 64 bytes");

// Quantize an in-memory log-prob table (g_ngram_logprob must be set; unseen cells
// sit at g_ngram_floor) to a dense 8-bit .ngbin at `path`. Returns 0 on success.
int write_ngram_bin(const char *path, const float *ngram_data, int ngram_size);

void perturbate_keyword(int state[], int len, int keyword_len);
void random_keyword(int keyword[], int len, int keyword_len);
void random_cycleword(int cycleword[], int max, int keyword_len);
void perturbate_cycleword(int state[], int max, int len);
int rand_int_frequency_weighted(int state[], int min_index, int max_index);
#endif
