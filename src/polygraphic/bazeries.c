//
//  Bazeries cipher primitives -- "simple substitution plus transposition" (ACA).
//
//  See bazeries.h for the full convention. The square build reuses bifid.c
//  (bifid_grid_from_keyword / bifid_build_inverse), so Bazeries adds almost no new cipher
//  math: it is a fixed monoalphabetic square-to-square substitution composed with a
//  digit-grouped reversal transposition, both keyed by one number N < 1,000,000.
//

#include "colossus.h"
#include "bazeries.h"

// Single-threaded scratch for the transposed/intermediate stream, kept off the stack so the
// per-iteration solver decrypt hook does not carry a MAX_CIPHER_LENGTH frame (cf. bifid.c).
static _Thread_local int g_baz_scratch[MAX_CIPHER_LENGTH];

// English number-word letters (uppercase; mapped through g_char_to_idx, which folds the
// active J->I alphabet -- no number word contains J). The transposition uses N's digits,
// not these words, so only the LETTERS matter for the keyed square.
static const char *const BAZ_ONES[20] = {
    "", "ONE", "TWO", "THREE", "FOUR", "FIVE", "SIX", "SEVEN", "EIGHT", "NINE",
    "TEN", "ELEVEN", "TWELVE", "THIRTEEN", "FOURTEEN", "FIFTEEN", "SIXTEEN",
    "SEVENTEEN", "EIGHTEEN", "NINETEEN"
};
static const char *const BAZ_TENS[10] = {
    "", "", "TWENTY", "THIRTY", "FORTY", "FIFTY", "SIXTY", "SEVENTY", "EIGHTY", "NINETY"
};

// Append the letters of word w to kw[] (mapping each through g_char_to_idx); advance *m.
static void baz_append_word(const char *w, int kw[], int *m) {
    for (const char *p = w; *p; p++) {
        int idx = g_char_to_idx[(unsigned char) *p];
        if (idx >= 0 && *m < BAZERIES_MAX_SPELL) kw[(*m)++] = idx;
    }
}

// Spell a value < 100.
static void baz_spell_below_100(int n, int kw[], int *m) {
    if (n < 20) { baz_append_word(BAZ_ONES[n], kw, m); return; }
    baz_append_word(BAZ_TENS[n / 10], kw, m);
    if (n % 10) baz_append_word(BAZ_ONES[n % 10], kw, m);
}

// Spell a value < 1000 (no "and"; "and" letters A,N,D already appear in THOUSAND, so the
// keyed square is identical either way).
static void baz_spell_below_1000(int n, int kw[], int *m) {
    if (n >= 100) {
        baz_append_word(BAZ_ONES[n / 100], kw, m);
        baz_append_word("HUNDRED", kw, m);
        n %= 100;
    }
    if (n) baz_spell_below_100(n, kw, m);
}

int bazeries_spell(long n, int kw[]) {
    int m = 0;
    if (n >= 1000) {
        baz_spell_below_1000((int) (n / 1000), kw, &m);
        baz_append_word("THOUSAND", kw, &m);
        n %= 1000;
        if (n) baz_spell_below_1000((int) n, kw, &m);
    } else {
        baz_spell_below_1000((int) n, kw, &m);
    }
    return m;
}

void bazeries_build_square(long key, int square[]) {
    int kw[BAZERIES_MAX_SPELL];
    int kwlen = bazeries_spell(key, kw);
    bifid_grid_from_keyword(kw, kwlen, square, BAZERIES_GRID);
}

void bazeries_build_sub(const int square[], int fsub[]) {
    // Plaintext square is column-major: letter L sits at (r=L%5, c=L/5) -> ct cell r*side+c.
    for (int L = 0; L < BAZERIES_GRID; L++) {
        int r = L % BAZERIES_SIDE, c = L / BAZERIES_SIDE;
        fsub[L] = square[r * BAZERIES_SIDE + c];
    }
}

void bazeries_build_invsub(const int square[], int invsub[]) {
    int pos[BAZERIES_GRID];
    bifid_build_inverse(square, pos, BAZERIES_GRID);          // ct letter -> ct cell
    for (int X = 0; X < BAZERIES_GRID; X++) {
        int p = pos[X], r = p / BAZERIES_SIDE, c = p % BAZERIES_SIDE;
        invsub[X] = c * BAZERIES_SIDE + r;                    // pt square's column-major letter
    }
}

void bazeries_digits(long key, int digits[], int *ndigits) {
    int tmp[BAZERIES_MAX_DIGITS + 4], t = 0;
    if (key <= 0) { digits[0] = 0; *ndigits = 1; return; }
    while (key > 0 && t < (int) (sizeof tmp / sizeof tmp[0])) { tmp[t++] = (int) (key % 10); key /= 10; }
    for (int i = 0; i < t; i++) digits[i] = tmp[t - 1 - i];   // reverse to MSB-first
    *ndigits = t;
}

void bazeries_transpose(const int in[], int n, const int digits[], int ndigits, int out[]) {
    int pos = 0, di = 0;
    while (pos < n) {
        int g = digits[di % ndigits];
        di++;
        if (g <= 0) continue;                                // zero-length group: skip
        if (g > n - pos) g = n - pos;                        // ragged final group
        for (int i = 0; i < g; i++) out[pos + i] = in[pos + g - 1 - i];   // reverse in place
        pos += g;
    }
}

void bazeries_encrypt(const int plain[], int n, long key, int out[]) {
    int square[BAZERIES_GRID], fsub[BAZERIES_GRID];
    int digits[BAZERIES_MAX_DIGITS], nd;
    bazeries_build_square(key, square);
    bazeries_build_sub(square, fsub);
    bazeries_digits(key, digits, &nd);

    int *rv = g_baz_scratch;
    bazeries_transpose(plain, n, digits, nd, rv);            // transpose (reverse digit groups)
    for (int i = 0; i < n; i++) out[i] = fsub[rv[i]];        // then substitute pt->ct
}

void bazeries_decrypt(const int cipher[], int n, long key, int out[]) {
    int square[BAZERIES_GRID], invsub[BAZERIES_GRID];
    int digits[BAZERIES_MAX_DIGITS], nd;
    bazeries_build_square(key, square);
    bazeries_build_invsub(square, invsub);
    bazeries_digits(key, digits, &nd);

    int *rv = g_baz_scratch;
    for (int i = 0; i < n; i++) rv[i] = invsub[cipher[i]];   // inverse-substitute ct->pt (RV)
    bazeries_transpose(rv, n, digits, nd, out);              // un-transpose (involution)
}
