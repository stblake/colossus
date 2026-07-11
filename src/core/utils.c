
// Cryptanalytic and general software utilities 


#include "colossus.h"

_Thread_local uint32_t rng_state = 123456789;

// --- Runtime alphabet (see colossus.h). Defaults to the full A..Z so the
// historical 26-letter behaviour is bit-identical until -excludeletter/-alphabet
// is given. ---
int  g_alpha = ALPHABET_SIZE;
int  g_char_to_idx[128];
char g_idx_to_char_arr[MAX_ALPHABET_SIZE + 1] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
double g_monograms[MAX_ALPHABET_SIZE];

// n-gram scoring mode. false (default) keeps the historical reward-only normalized
// log(1+count) table, so every existing solve is bit-identical. true selects the
// AZDecrypt / Practical-Cryptography fitness: log10 probabilities with a floor that
// PENALISES unseen n-grams (set via -logprob, recommended with quintgrams).
bool g_ngram_logprob = false;

// Reversal-invariant n-gram scoring (set via -reversengrams). When true, load_ngrams
// symmetrizes the table so every n-gram and its digit-reversed twin carry the same
// (max) weight -- a plaintext written with any words/segments reversed then scores like
// clean English. Default false => the table and every existing solve are bit-identical.
bool g_ngram_reverse = false;

// Set true once (on the main thread, in solve_cipher, after decode_cipher) iff the
// decoded cipher stream carries NO negative sentinels (spaces/punctuation) -- i.e. the
// scored plaintext is all letters. When true, ngram_score skips its letters-only
// compaction copy (the pass is the identity in that case) and scores the decrypt in
// place -- bit-identical, but one fewer O(n) pass per hill-climb iteration. Default
// false => the safe compaction path, byte-for-byte the historical scorer.
bool g_score_no_sentinel = false;

// Build the index<->char maps and the reindexed monogram table. `excluded` is a
// string of letters to drop from the standard A..Z ordering (NULL/"" => full A..Z).
void init_alphabet(const char *excluded) {
    for (int i = 0; i < 128; i++) g_char_to_idx[i] = -1;
    int pos = 0;
    for (char c = 'A'; c <= 'Z'; c++) {
        bool drop = false;
        if (excluded) {
            for (const char *e = excluded; *e; e++) {
                if (toupper((unsigned char) *e) == c) { drop = true; break; }
            }
        }
        if (drop) continue;
        g_idx_to_char_arr[pos] = c;
        g_char_to_idx[(int) c] = pos;
        g_monograms[pos] = english_monograms[c - 'A'];
        pos++;
    }
    g_idx_to_char_arr[pos] = '\0';
    g_alpha = pos;
}

// Trifid 27-symbol alphabet: the full A..Z (0..25) plus TRIFID_EXTRA_CHAR ('+') at
// index 26, so the 3x3x3 cube has exactly 27 cells. The extra symbol is registered in
// g_char_to_idx (so it decodes from the ciphertext, unlike a sentinel) and given a
// negligible monogram weight -- the cube attack is a bijection and uses no monogram
// penalty, and '+' never occurs in the A..Z plaintext the n-gram score sees.
void init_alphabet_trifid(void) {
    init_alphabet(NULL);                                  // A..Z -> 0..25, g_alpha = 26
    int pos = g_alpha;                                    // == DEFAULT_ALPHABET_SIZE (26)
    g_idx_to_char_arr[pos] = TRIFID_EXTRA_CHAR;
    g_idx_to_char_arr[pos + 1] = '\0';
    g_char_to_idx[(unsigned char) TRIFID_EXTRA_CHAR] = pos;
    g_monograms[pos] = 1e-6;                              // negligible (unused by Trifid)
    g_alpha = pos + 1;                                    // 27
}

// ADFGVX 36-symbol alphabet: the full A..Z (0..25) plus the digits '0'..'9' at indices
// 26..35, so the 6x6 Polybius square has exactly 36 cells. The digits are registered in
// g_char_to_idx (so a digit in the plaintext decodes, unlike a sentinel) and given a
// negligible monogram weight -- the square attack is a bijection and uses no monogram
// penalty, and a digit never occurs in the A..Z n-gram table the score consults.
void init_alphabet_adfgvx(void) {
    init_alphabet(NULL);                                  // A..Z -> 0..25, g_alpha = 26
    int pos = g_alpha;                                    // == DEFAULT_ALPHABET_SIZE (26)
    for (char d = '0'; d <= '9'; d++) {
        g_idx_to_char_arr[pos] = d;
        g_char_to_idx[(unsigned char) d] = pos;
        g_monograms[pos] = 1e-6;                          // negligible (unused by ADFGVX)
        pos++;
    }
    g_idx_to_char_arr[pos] = '\0';
    g_alpha = pos;                                        // 36
}

// Digrafid 27-symbol alphabet: the full A..Z (0..25) plus DIGRAFID_EXTRA_CHAR ('#') at
// index 26, so each 3x9 / 9x3 grid has exactly 27 cells. The '#' is registered in
// g_char_to_idx (so it decodes from the ciphertext, unlike a sentinel) and given a
// negligible monogram weight -- the grid attack is a bijection and uses no monogram
// penalty, and '#' never occurs in the A..Z plaintext the n-gram score sees.
void init_alphabet_digrafid(void) {
    init_alphabet(NULL);                                  // A..Z -> 0..25, g_alpha = 26
    int pos = g_alpha;                                    // == DEFAULT_ALPHABET_SIZE (26)
    g_idx_to_char_arr[pos] = DIGRAFID_EXTRA_CHAR;
    g_idx_to_char_arr[pos + 1] = '\0';
    g_char_to_idx[(unsigned char) DIGRAFID_EXTRA_CHAR] = pos;
    g_monograms[pos] = 1e-6;                              // negligible (unused by Digrafid)
    g_alpha = pos + 1;                                    // 27
}

void init_alphabet_ragbaby(void) {
    init_alphabet("JX");                                  // 24 letters: A..Z minus J and X
    // Re-register the two dropped letters as PAIRS so a plaintext/ciphertext J or X still
    // decodes (ACA Ragbaby pairs I/J and W/X). init_alphabet left both at -1.
    g_char_to_idx[(int) 'J'] = g_char_to_idx[(int) 'I'];
    g_char_to_idx[(int) 'X'] = g_char_to_idx[(int) 'W'];
}

int gcd(int a, int b) {
    while (b) { a %= b; int t = a; a = b; b = t; }
    return a;
}


double entropy(int text[], int len) {
    int frequencies[ALPHABET_SIZE];
    double entropy = 0., freq;

    tally(text, len, frequencies, ALPHABET_SIZE);

    for (int i = 0; i < ALPHABET_SIZE; i++) {
        if (frequencies[i] > 0) {
            freq = ((double) frequencies[i])/len;
            entropy -= freq*log(freq);
        }
    }
    return entropy;
}



double chi_squared(int plaintext[], int len) {
    // counts sized MAX_ALPHABET_SIZE: the loop runs to g_alpha, which can be 27 (Trifid),
    // so a 26-wide array would read one past the end. Bit-identical for alphabets <= 26.
    int i, counts[MAX_ALPHABET_SIZE];
    double frequency, chi2 = 0.;
    tally(plaintext, len, counts, MAX_ALPHABET_SIZE);
    for (i = 0; i < g_alpha; i++) {
        frequency = ((double) counts[i])/len;
        chi2 += pow(frequency - g_monograms[i], 2)/g_monograms[i];
    }
    return chi2;
}



int unique_len(char *str) {
    int len = 0;
    int seen[ALPHABET_SIZE];
    int i, idx;

    // Initialise seen array.
    for (i = 0; i < ALPHABET_SIZE; i++) seen[i] = 0;

    for (i = 0; str[i] != '\0'; i++) {
        unsigned char ch = (unsigned char) toupper((unsigned char) str[i]);
        idx = (ch < 128) ? g_char_to_idx[ch] : -1;
        if (idx >= 0 && idx < ALPHABET_SIZE) {
            if (!seen[idx]) {
                seen[idx] = 1;
                len++;
            }
        }
    }
    return len;
}



double vec_mean(double vec[], int len) {
    int i;
    double total = 0.;
    for (i = 0; i < len; i++) total += vec[i];
    return total/len;
}



double vec_stddev(double vec[], int len) {
    int i;
    double mu, sumdev = 0.;
    mu = vec_mean(vec, len);
    for (i = 0; i < len; i++) sumdev += pow(vec[i] - mu, 2); 
    return sqrt(sumdev/len);
}



void vec_print(int vec[], int len) {
    for (int i = 0; i < len; i++) printf("%d ", vec[i]);
    printf("\n");
}



void print_text(int indices[], int len) {
    for (int i = 0; i < len; i++) printf("%c", index_to_char(indices[i]));
}



void ord(char *text, int indices[]) {
    for (int i = 0; i < strlen(text); i++) {
        unsigned char c = (unsigned char) text[i];
        // Consult g_char_to_idx for any ASCII char (not just isalpha): it is -1 for
        // every unregistered byte, so A..Z inputs are byte-for-byte unchanged, but a
        // registered non-letter (the Trifid '+') maps to its index instead of a
        // sentinel. Letters outside the runtime alphabet and all unregistered
        // non-letters are carried as reversible negative sentinels.
        int v = (c < 128) ? g_char_to_idx[toupper(c)] : -1;
        indices[i] = (v >= 0) ? v : (-(int) c - 1);
    }
}



// Map a single character to a letter index, the way ord() does: A..Z (runtime
// alphabet) -> 0..g_alpha-1, anything else a reversible negative sentinel. Any ASCII
// char registered in g_char_to_idx maps (so the Trifid '+' decodes); unregistered
// bytes stay -1, keeping every A..Z decode byte-for-byte identical.
static int char_to_index(unsigned char c) {
    int v = (c < 128) ? g_char_to_idx[toupper(c)] : -1;
    return (v >= 0) ? v : (-(int) c - 1);
}

// The field separator to use for this decode. An explicit -delimiter wins; otherwise
// a homophonic ciphertext containing commas auto-selects ',' (its conventional form)
// and everything else stays per-character (0 => no delimiter).
static char resolve_delimiter(const ColossusConfig *cfg, const char *text) {
    if (cfg->delimiter_present) return cfg->delimiter;
    if (cfg->cipher_type == HOMOPHONIC && strchr(text, ',')) return ',';
    return 0;
}

// Intern a token into the symbol table, returning its id (new id if unseen).
static int symbol_intern(SymbolTable *tab, const char *tok) {
    for (int i = 0; i < tab->n; i++)
        if (strcmp(tab->tokens[i], tok) == 0) { tab->freq[i]++; return i; }
    if (tab->n >= MAX_SYMBOLS) {
        fprintf(stderr, "\nERROR: ciphertext uses more than %d distinct symbols.\n", MAX_SYMBOLS);
        exit(1);
    }
    int id = tab->n++;
    strncpy(tab->tokens[id], tok, MAX_TOKEN_LEN - 1);
    tab->tokens[id][MAX_TOKEN_LEN - 1] = '\0';
    tab->freq[id] = 1;
    return id;
}

int decode_cipher(const char *text, const ColossusConfig *cfg, int indices[], SymbolTable *tab) {
    char delim = resolve_delimiter(cfg, text);
    int n = 0;

    if (cfg->cipher_type != HOMOPHONIC) {
        // LETTER mapping: each token resolves to a 0..g_alpha-1 letter (or a sentinel).
        if (delim == 0) {
            // Per character -- byte-for-byte ord(); preserves spaces/punctuation as
            // sentinels so the existing letter ciphers stay bit-identical.
            for (int i = 0; text[i] != '\0'; i++)
                indices[n++] = char_to_index((unsigned char) text[i]);
        } else {
            // Delimited letters (e.g. -delimiter , over "H,E,L,L,O"): one field per
            // letter, whitespace trimmed.
            for (int i = 0; ; i++) {
                char ch = text[i];
                if (ch == delim || ch == '\0') {
                    // (a run of separators / a trailing separator yields no field)
                    if (ch == '\0') break;
                } else if (!isspace((unsigned char) ch)) {
                    indices[n++] = char_to_index((unsigned char) ch);
                }
            }
        }
        return n;
    }

    // SYMBOL mapping (homophonic): split into surface tokens and intern them.
    if (tab == NULL) return 0;
    tab->n = 0;
    tab->delimiter = delim;

    if (delim == 0) {
        // One symbol per non-whitespace character.
        for (int i = 0; text[i] != '\0'; i++) {
            unsigned char c = (unsigned char) text[i];
            if (isspace(c)) continue;
            char tok[2] = { (char) c, '\0' };
            indices[n++] = symbol_intern(tab, tok);
        }
    } else {
        char tok[MAX_TOKEN_LEN];
        int t = 0;
        for (int i = 0; ; i++) {
            char ch = text[i];
            if (ch == delim || ch == '\0') {
                if (t > 0) { tok[t] = '\0'; indices[n++] = symbol_intern(tab, tok); t = 0; }
                if (ch == '\0') break;
            } else if (!isspace((unsigned char) ch)) {
                if (t < MAX_TOKEN_LEN - 1) tok[t++] = ch;
            }
        }
    }
    return n;
}

void print_cipher(const int indices[], int len, const SymbolTable *tab) {
    if (tab == NULL) { print_text((int *) indices, len); return; }
    char sep = tab->delimiter ? tab->delimiter : ' ';
    for (int i = 0; i < len; i++) {
        if (i) printf("%c", sep);
        printf("%s", tab->tokens[indices[i]]);
    }
}



void tally(int plaintext[], int len, int frequencies[], int n_frequencies) {
    int i;
    for (i = 0; i < n_frequencies; i++) frequencies[i] = 0;
    // Skip negative sentinels (non-alphabetic characters); they have no histogram
    // bin. For all-letter text this loop is identical to the unguarded version.
    for (i = 0; i < len; i++) {
        if (plaintext[i] >= 0 && plaintext[i] < n_frequencies) frequencies[plaintext[i]]++;
    }
}



float index_of_coincidence(int plaintext[], int len) {
    // Sized MAX_ALPHABET_SIZE so a live symbol id up to g_alpha-1 (the Trifid '+' at 26,
    // present in the ciphertext the period estimator scans via mean_ioc) is counted, not
    // dropped. For any alphabet <= 26 the extra bin stays 0, so the IoC is bit-identical.
    int i, frequencies[MAX_ALPHABET_SIZE];
    double ioc = 0.;
    tally(plaintext, len, frequencies, MAX_ALPHABET_SIZE);
    for (i = 0; i < MAX_ALPHABET_SIZE; i++) ioc += frequencies[i]*(frequencies[i] - 1);
    ioc /= len*(len - 1);
    return ioc;
}



void straight_alphabet(int keyword[], int len) {
    for (int i = 0; i < len; i++) keyword[i] = i;
}

// Build a keyed alphabet (as 0-25 indices) from a keyword string: the keyword's
// distinct letters in order, then the remaining letters in alphabetical order.
void make_keyed_alphabet(char *keyword_str, int *output_indices) {
    int i, char_idx;
    int seen[ALPHABET_SIZE];
    int current_pos = 0;
    int len = strlen(keyword_str);

    for(i = 0; i < ALPHABET_SIZE; i++) seen[i] = 0;

    // Process keyword string.
    for(i = 0; i < len; i++) {
        char_idx = g_char_to_idx[toupper((unsigned char) keyword_str[i]) & 127];
        if(char_idx >= 0 && char_idx < ALPHABET_SIZE) {
            if(!seen[char_idx]) {
                output_indices[current_pos++] = char_idx;
                seen[char_idx] = 1;
            }
        }
    }

    // Fill remaining alphabet.
    for(i = 0; i < g_alpha; i++) {
        if(!seen[i]) {
            output_indices[current_pos++] = i;
        }
    }
}



bool file_exists(const char *filename) {
    FILE *file;
    file = fopen(filename, "r");
    if (file) {
        fclose(file);
        return true;
    }
    return false;
}



void shuffle(int *array, size_t n) {
    if (n > 1) {
        size_t i;
        for (i = 0; i < n - 1; i++) {
          // size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          size_t range = n - i;
          size_t j = i + rand_bounded(range);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}



void vec_copy(int src[], int dest[], int len) {
    for (int i = 0; i < len; i++) dest[i] = src[i]; 
}



int int_pow(int base, int exp) {
    int result = 1;
    while (exp) {
        if (exp % 2) result *= base;
        exp /= 2;
        if (exp) base *= base; // skip the final squaring (would overflow, e.g. 26^4 squared)
    }
    return result;
}
