#include "keyphrase.h"
#include <ctype.h>

// =====================================================================
//  Key Phrase cipher primitives -- see keyphrase.h
// =====================================================================

void keyphrase_encrypt(const int in[], int len, const int key[KEYPHRASE_LEN], int out[]) {
    for (int i = 0; i < len; i++)
        out[i] = (in[i] >= 0 && in[i] < KEYPHRASE_LEN) ? key[in[i]] : in[i];
}

int keyphrase_build_key(const char *phrase, int key[KEYPHRASE_LEN]) {
    int n = 0;
    for (int i = 0; phrase[i] && n < KEYPHRASE_LEN; i++) {
        unsigned char c = (unsigned char) toupper((unsigned char) phrase[i]);
        if (c >= 'A' && c <= 'Z') key[n++] = g_char_to_idx[c];
    }
    return n;
}
