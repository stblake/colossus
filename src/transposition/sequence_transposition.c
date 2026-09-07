//
// Sequence Transposition cipher primitives (ACA)
//

#include "colossus.h"

/*
   Sequence Transposition (ACA)

   A COLUMNAR-style transposition whose column assignment is driven by a Gromark-style
   chain-addition digit sequence rather than by a grid of fixed width.

   1. A P-digit primer (the ACA standard is P = 5) seeds a running digit sequence the length
      of the plaintext by chain addition -- the SAME recurrence Gromark uses:

          SS[i] = primer[i]                        for i < P
          SS[i] = (SS[i-P] + SS[i-P+1]) mod 10     for i >= P

      (add successive pairs of digits, dropping tens). One digit 0..9 per plaintext letter.

   2. There are 10 columns, one per digit value 0..9. Plaintext letter i belongs to column
      SS[i]. Within a column the letters keep their plaintext order.

   3. A 10-letter keyword is numbered by stable alphabetical rank (1..10, ties left to right;
      rank 10 is written as digit 0). Keyword position k heads the column of that digit, and
      the columns are drawn off in keyword-position order 0,1,...,9. So the ciphertext is the
      concatenation of the 10 columns in the order given by the read-order permutation
      `pi` (pi[k] = the digit column emitted at read step k, a permutation of 0..9).

   ACA worked example: plaintext THE EARLY BIRD GETS THE WORM, keyword GUMMYBEARS, primer
   69315 -> SS = 6931552460760673663092 9, pi = [4,9,5,6,0,2,3,1,7,8], and the columns drawn
   off in that order give the ciphertext YHOMARTBDETHIGWLRESEERT. (The transmitted form also
   prepends the primer and appends a 1-digit checksum, which are transmission artifacts, not
   part of the transposition -- the solver works on the bare-letter body.)

   The transposition is a pure permutation of the symbols, so these primitives are
   index-agnostic: they move whatever integer symbols are handed to them (0..25 alphabet
   indices in normal use), and encrypt/decrypt round-trip exactly.
*/

// Read-order permutation from a keyword: pi[k] = the digit column emitted at read step k.
// The first SEQ_TRANS_BUCKETS letters of the keyword are ranked 1..10 by stable alphabetical
// order (ties broken left to right); rank r maps to digit (r mod 10) so rank 10 -> digit 0.
// A proper 10-letter keyword yields a permutation of 0..9; any unfilled slots keep the
// identity order (a malformed short keyword is the caller's responsibility).
void sequence_transposition_pi_from_keyword(const char *keyword, int pi[]) {
    char letters[SEQ_TRANS_BUCKETS];
    int m = 0;
    for (int i = 0; keyword[i] && m < SEQ_TRANS_BUCKETS; i++) {
        unsigned char c = (unsigned char) keyword[i];
        if (isalpha(c)) letters[m++] = (char) toupper(c);
    }
    for (int k = 0; k < SEQ_TRANS_BUCKETS; k++) pi[k] = k;   // identity default
    for (int k = 0; k < m; k++) {
        int rank = 1;
        for (int j = 0; j < m; j++)
            if (letters[j] < letters[k] || (letters[j] == letters[k] && j < k)) rank++;
        pi[k] = rank % SEQ_TRANS_BUCKETS;                    // 1..9 stay; 10 -> 0
    }
}

// Encrypt: gather the plaintext column by column (all SS==pi[0], then all SS==pi[1], ...),
// each column in plaintext order. out[] holds `len` symbols.
void sequence_transposition_encrypt(const int plain[], int len, const int primer[],
                     int primer_len, const int pi[], int out[]) {
    static _Thread_local int ss[MAX_CIPHER_LENGTH];
    gromark_chain_key(primer, primer_len, len, ss);
    int pos = 0;
    for (int k = 0; k < SEQ_TRANS_BUCKETS; k++) {
        int d = pi[k];
        for (int i = 0; i < len; i++)
            if (ss[i] == d) out[pos++] = plain[i];
    }
}

// Decrypt: recompute SS[], count each column, lay the columns out in pi read order to find
// each column's start offset in the ciphertext, then redistribute the ciphertext symbols back
// to their plaintext positions. O(len).
void sequence_transposition_decrypt(const int cipher[], int len, const int primer[],
                     int primer_len, const int pi[], int out[]) {
    static _Thread_local int ss[MAX_CIPHER_LENGTH];
    gromark_chain_key(primer, primer_len, len, ss);
    int count[SEQ_TRANS_BUCKETS] = {0};
    for (int i = 0; i < len; i++) count[ss[i]]++;
    int next[SEQ_TRANS_BUCKETS], run = 0;
    for (int k = 0; k < SEQ_TRANS_BUCKETS; k++) { int d = pi[k]; next[d] = run; run += count[d]; }
    for (int i = 0; i < len; i++) out[i] = cipher[next[ss[i]]++];
}
