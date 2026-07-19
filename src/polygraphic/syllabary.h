#ifndef SYLLABARY_H
#define SYLLABARY_H
#include "colossus.h"

// =====================================================================
//  Syllabary cipher primitive (syllabary.c)
// =====================================================================
//
// An ACA "Syllabary" (110-154 ciphertext pairs). A 10x10 square holds a FIXED standard
// alphabet of 100 syllabary elements -- the 26 single letters, the digits 1-9, a null
// element (here written as an empty token), and 64 common syllables (AR, ING, THE, RED,
// ...). Row and column labels are permutations of 0-9. Each plaintext ELEMENT (a syllable
// or single letter) is enciphered by the 2-digit (row-label, column-label) code of its
// cell. Because a plaintext run can be split into elements in several ways (o+r+d+e+r vs
// or+de+r), the SAME plaintext has many spellings (isologs), which flattens letter and
// word frequencies.
//
//   Square (ACA "Unknown coordinates, Unknown keysquare" worked example):
//          6  7  1  9  4  3  2  5  0  8       <- column labels
//     8  | C  3  H  8  AR M  ING P  RI N
//     5  | CE A  1  AL AN AND ARE AS AT ATE
//     0  | ATI B  2  BE CA CO COM D  4  DA
//     2  | DE E  5  EA ED EN ENT ER ERE ERS
//     3  | ES EST F  6  G  7  HAS HE I  9
//     4  | IN ION IS IT IVE J  .  K  L  LA     ( . = the null element )
//     1  | LE ME ND NE NT O  OF ON OR OU
//     6  | Q  R  RA RE RED RES RO S  SE SH
//     7  | ST STO T  TE TED TER TH THE THI THR
//     9  | TI TO U  V  VE W  WE X  Y  Z
//   ^ row labels
//
//   pt "orders received" -> 13 67 05 27 67 65 67 27 86 27 30 99 27 05 (the all-single-letter
//   spelling; or+de+re+ceived etc. give shorter isologous ciphertexts).
//
// DECODING IS A UNIQUE BIJECTION (code -> the token in that cell); ENCODING is one-to-many
// (variant spellings). So cryptanalytically Syllabary is a SUBSTITUTION over 100 numeric
// codes -> 100 KNOWN multi-character tokens (a length-changing decode). The row/column
// label order folds into the square (like the Checkerboard / Nihilist-Sub), so the solver
// searches the composite code->token map directly (see syllabary_solver.c) -- this subsumes
// the ACA "Known/Unknown Coordinates x Known/Unknown Keysquare" variants uniformly.
//
// The canonical UNMIXED token order below is an approximation of the ACA appendix (not in
// the type PDF); it drives only the generator's keyed fill and is irrelevant to the solver,
// which searches the composite map. The KAT builds the exact PDF square from token strings.

#define SYL_NTOKENS   100
#define SYL_SIDE      10
#define SYL_MAXTOK    3          // longest token is 3 chars (STO, THE, THR, ...)
#define SYL_NULL_IDX  99         // canonical index of the null element (decodes to no letters)

// The fixed 100-element syllabary alphabet, in canonical (unmixed) order. Index 0..25 are the
// single letters A..Z, 26..34 the digits 1..9, 35..98 the 64 syllables, 99 the null ("").
extern const char *const syllabary_tokens[SYL_NTOKENS];

// Token index of a string (case-sensitive uppercase; "" -> the null index), or -1 if absent.
int syllabary_token_index(const char *s);

// A built square: the token arrangement, the label permutations, and the derived code maps.
typedef struct {
    int token_of_cell[SYL_NTOKENS];  // physical cell r*10+c -> token index 0..99 (a permutation)
    int row_label[SYL_SIDE];         // physical row -> its digit label (a permutation of 0..9)
    int col_label[SYL_SIDE];         // physical col -> its digit label (a permutation of 0..9)
    int token_of_code[100];          // code rlab*10+clab -> token index (decode map)
    int code_of_token[SYL_NTOKENS];  // token index -> its 2-digit code (encode map)
} SyllabarySquare;

// Tokenization modes for encrypt.
enum { SYL_TOK_SINGLE = 0,   // every plaintext letter as its single-letter token (unambiguous)
       SYL_TOK_GREEDY = 1,   // greedy longest matching token at each position
       SYL_TOK_RANDOM = 2 }; // random matching token at each position (isolog frequency flattening)

// Build the square from the physical-order token indices `cell_tokens[0..99]` (a permutation of
// 0..99, cell r*10+c) and the two label permutations. Returns 0 on success, -1 on a non-permutation.
int syllabary_build(SyllabarySquare *sq, const int cell_tokens[SYL_NTOKENS],
                    const int row_label[SYL_SIDE], const int col_label[SYL_SIDE]);

// Convenience builder from token STRINGS in physical order (row-major, row 0 first). Looks each
// up in the canonical table. Returns 0 on success, -1 if a string is unknown or a label bad.
int syllabary_build_from_strings(SyllabarySquare *sq, const char *const cell_str[SYL_NTOKENS],
                                 const int row_label[SYL_SIDE], const int col_label[SYL_SIDE]);

// A reproducible random square (random token arrangement + random label permutations) from
// `seed` -- the general "Unknown Coordinates, Unknown Keysquare" case. Uses rand_int.
void syllabary_build_random(SyllabarySquare *sq, uint32_t seed);

// Encipher plaintext letters `plain[0..n-1]` (indices 0..25) into the 2-digit code stream
// out_codes[] using tokenization `mode`. Returns the number of codes written.
int syllabary_encrypt(const int plain[], int n, const SyllabarySquare *sq, int mode, int out_codes[]);

// Decipher the 2-digit code stream `codes[0..ncodes-1]` into out_pt[] (letters 0..25): each code
// -> its cell's token -> that token's letters, concatenated (the null token contributes none).
// Returns the number of plaintext letters written.
int syllabary_decrypt(const int codes[], int ncodes, const SyllabarySquare *sq, int out_pt[]);

#endif
