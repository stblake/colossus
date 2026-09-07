#include "syllabary.h"

// =====================================================================
//  Syllabary cipher primitive -- see syllabary.h
// =====================================================================

// The fixed 100-element syllabary alphabet in canonical order: 26 single letters, digits 1-9,
// 64 syllables (alphabetical), and the null element (empty token). Order is arbitrary but fixed
// -- only the generator's keyed fill depends on it; the solver searches the composite map.
const char *const syllabary_tokens[SYL_NTOKENS] = {
    "A","B","C","D","E","F","G","H","I","J","K","L","M",
    "N","O","P","Q","R","S","T","U","V","W","X","Y","Z",
    "1","2","3","4","5","6","7","8","9",
    "AL","AN","AND","AR","ARE","AS","AT","ATE","ATI","BE","CA","CE","CO","COM",
    "DA","DE","EA","ED","EN","ENT","ER","ERE","ERS","ES","EST","HAS","HE",
    "IN","ING","ION","IS","IT","IVE","LA","LE","ME","ND","NE","NT",
    "OF","ON","OR","OU","RA","RE","RED","RES","RI","RO","SE","SH","ST","STO",
    "TE","TED","TER","TH","THE","THI","THR","TI","TO","VE","WE",
    ""   // SYL_NULL_IDX (99): the null element, decodes to no letters
};

int syllabary_token_index(const char *s) {
    for (int i = 0; i < SYL_NTOKENS; i++)
        if (strcmp(syllabary_tokens[i], s) == 0) return i;
    return -1;
}

// True if int array a[0..n-1] is a permutation of 0..n-1.
static int is_perm(const int a[], int n) {
    int seen[SYL_NTOKENS];
    for (int i = 0; i < n; i++) seen[i] = 0;
    for (int i = 0; i < n; i++) {
        if (a[i] < 0 || a[i] >= n || seen[a[i]]) return 0;
        seen[a[i]] = 1;
    }
    return 1;
}

int syllabary_build(SyllabarySquare *sq, const int cell_tokens[SYL_NTOKENS],
                    const int row_label[SYL_SIDE], const int col_label[SYL_SIDE]) {
    if (!is_perm(cell_tokens, SYL_NTOKENS) || !is_perm(row_label, SYL_SIDE) || !is_perm(col_label, SYL_SIDE))
        return -1;
    for (int i = 0; i < SYL_NTOKENS; i++) sq->token_of_cell[i] = cell_tokens[i];
    for (int i = 0; i < SYL_SIDE; i++) { sq->row_label[i] = row_label[i]; sq->col_label[i] = col_label[i]; }
    for (int r = 0; r < SYL_SIDE; r++)
        for (int c = 0; c < SYL_SIDE; c++) {
            int tok  = cell_tokens[r * SYL_SIDE + c];
            int code = row_label[r] * 10 + col_label[c];
            sq->token_of_code[code] = tok;
            sq->code_of_token[tok]  = code;
        }
    return 0;
}

int syllabary_build_from_strings(SyllabarySquare *sq, const char *const cell_str[SYL_NTOKENS],
                                 const int row_label[SYL_SIDE], const int col_label[SYL_SIDE]) {
    int cell_tokens[SYL_NTOKENS];
    for (int i = 0; i < SYL_NTOKENS; i++) {
        int t = syllabary_token_index(cell_str[i]);
        if (t < 0) return -1;
        cell_tokens[i] = t;
    }
    return syllabary_build(sq, cell_tokens, row_label, col_label);
}

void syllabary_build_random(SyllabarySquare *sq, uint32_t seed) {
    seed_rand(seed);
    int cell_tokens[SYL_NTOKENS], rl[SYL_SIDE], cl[SYL_SIDE];
    for (int i = 0; i < SYL_NTOKENS; i++) cell_tokens[i] = i;
    for (int i = SYL_NTOKENS - 1; i > 0; i--) { int j = rand_int(0, i + 1); int t = cell_tokens[i]; cell_tokens[i] = cell_tokens[j]; cell_tokens[j] = t; }
    for (int i = 0; i < SYL_SIDE; i++) { rl[i] = i; cl[i] = i; }
    for (int i = SYL_SIDE - 1; i > 0; i--) { int j = rand_int(0, i + 1); int t = rl[i]; rl[i] = rl[j]; rl[j] = t; }
    for (int i = SYL_SIDE - 1; i > 0; i--) { int j = rand_int(0, i + 1); int t = cl[i]; cl[i] = cl[j]; cl[j] = t; }
    syllabary_build(sq, cell_tokens, rl, cl);
}

// Does token `t` (a non-empty string) match the plaintext letters at plain[i..], within [i,n)?
static int token_matches(int t, const int plain[], int i, int n) {
    const char *s = syllabary_tokens[t];
    int L = (int) strlen(s);
    if (L == 0 || i + L > n) return 0;
    for (int k = 0; k < L; k++)
        if (s[k] < 'A' || s[k] > 'Z' || plain[i + k] != s[k] - 'A') return 0;
    return L;
}

int syllabary_encrypt(const int plain[], int n, const SyllabarySquare *sq, int mode, int out_codes[]) {
    int m = 0, i = 0;
    while (i < n) {
        if (plain[i] < 0 || plain[i] >= ALPHABET_SIZE) { i++; continue; }   // skip non-letters
        int chosen = plain[i], clen = 1;                                    // single-letter fallback
        if (mode == SYL_TOK_GREEDY) {
            for (int len = SYL_MAXTOK; len >= 2; len--) {
                int found = -1;
                for (int t = 0; t < SYL_NTOKENS && found < 0; t++)
                    if ((int) strlen(syllabary_tokens[t]) == len && token_matches(t, plain, i, n)) found = t;
                if (found >= 0) { chosen = found; clen = len; break; }
            }
        } else if (mode == SYL_TOK_RANDOM) {
            int cand[SYL_NTOKENS], clen_of[SYL_NTOKENS], nc = 0;
            for (int t = 0; t < SYL_NTOKENS; t++) {
                int L = token_matches(t, plain, i, n);
                if (L > 0) { cand[nc] = t; clen_of[nc] = L; nc++; }
            }
            if (nc > 0) { int pick = rand_int(0, nc); chosen = cand[pick]; clen = clen_of[pick]; }
        }
        out_codes[m++] = sq->code_of_token[chosen];
        i += clen;
    }
    return m;
}

int syllabary_decrypt(const int codes[], int ncodes, const SyllabarySquare *sq, int out_pt[]) {
    int m = 0;
    for (int i = 0; i < ncodes; i++) {
        int code = codes[i];
        if (code < 0 || code > 99) continue;
        const char *s = syllabary_tokens[sq->token_of_code[code]];
        for (int k = 0; s[k]; k++)
            if (s[k] >= 'A' && s[k] <= 'Z') out_pt[m++] = s[k] - 'A';
    }
    return m;
}
