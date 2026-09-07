// =====================================================================
//  Enigma machine primitive
// =====================================================================
//
// See enigma.h for the model overview. This file owns the reference wiring/notch tables
// and the encipher/stepping math. The three classic correctness traps are all handled and
// tested (tests/test_enigma.c):
//   1. The DOUBLE-STEP anomaly: when the middle rotor sits on its own turnover notch, the
//      next keypress steps BOTH the middle and the left rotor (the middle "double steps"),
//      because the pawl engages the ratchet on both sides.
//   2. MULTI-NOTCH naval rotors VI-VIII turn the next rotor over at TWO window letters
//      (Z and M), not one.
//   3. The ring/notch/position relationship: the turnover fires on the WINDOW letter
//      (`pos`), independent of the ring; the wiring transform uses the offset (pos - ring).
//
// Wiring strings are the A-> mapping read at ETW entry. The reflector is an involution.
//
// PERFORMANCE: the per-char rotor transform folds the whole shift-in / wiring-lookup / shift-out
// chain into precomputed per-offset tables (g_fwd_off/g_inv_off, built in enigma_init), so a rotor
// pass is a single table lookup with NO %26 -- this is the dominant cost of every ciphertext-only
// IoC decrypt. The tables assume the documented contract that letters and pos/ring are in [0,25]
// (the old %26 form silently masked out-of-range values; the table form does not). enigma_encrypt
// is the hot loop; the solver additionally caches whole per-position scramblers (see
// enigma_solver.c / enigma_bombe.c) so a fixed rotor config decodes without re-walking the rotors.

#include <string.h>
#include <strings.h>
#include "enigma.h"

// --- Reference data (Wikipedia "Enigma rotor details" / Crypto Museum) ----------------

static const char *ROTOR_WIRING[ENIGMA_N_ROTORS] = {
    "EKMFLGDQVZNTOWYHXUSPAIBRCJ",  // I
    "AJDKSIRUXBLHWTMCQGZNPYFVOE",  // II
    "BDFHJLCPRTXVZNYEIWGAKMUSQO",  // III
    "ESOVPZJAYQUIRHXLNFTGKDCMWB",  // IV
    "VZBRGITYUPSDNHLXAWMJQOFECK",  // V
    "JPGVOUMFYQBENHZRDKASXLICTW",  // VI
    "NZJHGRCXMYSWBOUFAIVLPEKQDT",  // VII
    "FKQHTLXOCBJSPDZRAMEWNIUYGV",  // VIII
    "LEYJVCNIXWPBQMDRTAKZGFUHOS",  // Beta (Greek)
    "FSOKANUERHMBTIYCWLQPZXVGJD",  // Gamma (Greek)
};

// Turnover WINDOW letters per rotor. I-V: one notch; VI-VIII: two (Z and M); Greek: none.
static const char *ROTOR_NOTCH[ENIGMA_N_ROTORS] = {
    "Q",  "E",  "V",  "J",  "Z",          // I-V
    "ZM", "ZM", "ZM",                     // VI-VIII
    "",   "",                             // Beta, Gamma
};

static const char *ROTOR_NAME[ENIGMA_N_ROTORS] = {
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "Beta", "Gamma",
};

static const char *REFLECTOR_WIRING[ENIGMA_N_REFLECTORS] = {
    "YRUHQSLDPXNGOKMIEBFZCWVJAT",  // UKW-B
    "FVPJIAOYEDRZXWGCTKUQSBNMHL",  // UKW-C
    "ENKQAUYWJICOPBLMDXZVFTHRGS",  // UKW-B thin (M4)
    "RDOBJNTKVEHMLFCWZAXGYIPSUQ",  // UKW-C thin (M4)
};

static const char *REFLECTOR_NAME[ENIGMA_N_REFLECTORS] = {
    "B", "C", "B-thin", "C-thin",
};

// --- Compiled tables (built once by enigma_init) --------------------------------------

static int  g_fwd[ENIGMA_N_ROTORS][26];       // A-> forward permutation (index in)
static int  g_inv[ENIGMA_N_ROTORS][26];       // inverse permutation (return path)
static int  g_notch[ENIGMA_N_ROTORS][26];     // g_notch[r][pos] != 0 => turnover at pos
static int  g_refl[ENIGMA_N_REFLECTORS][26];  // reflector involutions
static bool g_inited = false;

// Offset-folded wiring tables (the per-char hot path). For a rotor at wiring offset
// s = (pos - ring) mod 26 in [0,25], g_fwd_off[r][s][x] == the OLD rotor_forward result
// for every x in [0,25] (and likewise g_inv_off for rotor_backward). Precomputing the whole
// shift-in / lookup / shift-out chain removes all %26 from the encipher loop -- a single
// dependent load per rotor transform. Read-only after enigma_init; ~54 KB total.
static int  g_fwd_off[ENIGMA_N_ROTORS][26][26];
static int  g_inv_off[ENIGMA_N_ROTORS][26][26];

void enigma_init(void) {
    if (g_inited) return;
    for (int r = 0; r < ENIGMA_N_ROTORS; r++) {
        for (int i = 0; i < 26; i++) {
            int m = ROTOR_WIRING[r][i] - 'A';
            g_fwd[r][i] = m;
            g_inv[r][m] = i;
        }
        for (int i = 0; i < 26; i++) g_notch[r][i] = 0;
        for (const char *n = ROTOR_NOTCH[r]; *n; n++) g_notch[r][*n - 'A'] = 1;
    }
    for (int f = 0; f < ENIGMA_N_REFLECTORS; f++)
        for (int i = 0; i < 26; i++) g_refl[f][i] = REFLECTOR_WIRING[f][i] - 'A';
    // Fold the shift math into per-offset lookup tables (see the g_fwd_off comment).
    for (int r = 0; r < ENIGMA_N_ROTORS; r++)
        for (int s = 0; s < 26; s++)
            for (int x = 0; x < 26; x++) {
                g_fwd_off[r][s][x] = ((g_fwd[r][(x + s) % 26] - s) % 26 + 26) % 26;
                g_inv_off[r][s][x] = ((g_inv[r][(x + s) % 26] - s) % 26 + 26) % 26;
            }
    g_inited = true;
}

// --- Wiring transforms -----------------------------------------------------------------
//
// pos and ring are invariants in [0,25] (maintained by enigma_step's %26 and by init), so
// s = pos - ring is in [-25,25]; the single conditional lands it in [0,25] = the table index.
// Bit-identical to the historical ((pos-ring)%26 ...) form for every reachable argument.

static inline int rotor_forward(int rotor_id, int pos, int ring, int x) {
    int s = pos - ring; if (s < 0) s += 26;
    return g_fwd_off[rotor_id][s][x];
}

static inline int rotor_backward(int rotor_id, int pos, int ring, int x) {
    int s = pos - ring; if (s < 0) s += 26;
    return g_inv_off[rotor_id][s][x];
}

// --- Stepping --------------------------------------------------------------------------

static inline int step26(int x) { return x == 25 ? 0 : x + 1; }   // (x+1)%26, x in [0,25]

void enigma_step(EnigmaKey *key) {
    int n = key->n_wheels;
    int r = n - 1, m = n - 2, l = n - 3;   // right(fast), middle, left(slow) STEPPING wheels
    int rr = key->rotor[r], rm = key->rotor[m];
    int mid_at_notch   = g_notch[rm][key->pos[m]];
    int right_at_notch = g_notch[rr][key->pos[r]];
    if (mid_at_notch) {                    // double step: middle carries the left, and itself
        key->pos[m] = step26(key->pos[m]);
        key->pos[l] = step26(key->pos[l]);
    } else if (right_at_notch) {
        key->pos[m] = step26(key->pos[m]);
    }
    key->pos[r] = step26(key->pos[r]);     // the fast rotor always steps
}

// --- Encipher --------------------------------------------------------------------------

int enigma_encipher_letter(const EnigmaKey *key, int c) {
    int n = key->n_wheels;
    c = key->plug[c];
    for (int i = n - 1; i >= 0; i--)                        // right -> left, into reflector
        c = rotor_forward(key->rotor[i], key->pos[i], key->ring[i], c);
    c = g_refl[key->reflector][c];
    for (int i = 0; i < n; i++)                             // left -> right, back out
        c = rotor_backward(key->rotor[i], key->pos[i], key->ring[i], c);
    c = key->plug[c];
    return c;
}

void enigma_encrypt(const int *in, int n, const EnigmaKey *key, int *out) {
    EnigmaKey k = *key;                    // local copy: advance pos without touching caller
    for (int i = 0; i < n; i++) {
        enigma_step(&k);
        out[i] = enigma_encipher_letter(&k, in[i]);
    }
}

void enigma_scrambler_at(const EnigmaKey *base, int offset, int perm[26]) {
    EnigmaKey k = *base;
    for (int s = 0; s <= offset; s++) enigma_step(&k);     // step (offset+1) times
    int save[26];
    memcpy(save, k.plug, sizeof(save));
    enigma_plug_identity(k.plug);                          // scrambler excludes the plugboard
    for (int c = 0; c < 26; c++) perm[c] = enigma_encipher_letter(&k, c);
    memcpy(k.plug, save, sizeof(save));
}

// --- Plugboard -------------------------------------------------------------------------

void enigma_plug_identity(int plug[26]) {
    for (int i = 0; i < 26; i++) plug[i] = i;
}

int enigma_plug_set_pair(int plug[26], int a, int b) {
    if (a == b || a < 0 || b < 0 || a >= 26 || b >= 26) return -1;
    if (plug[a] != a || plug[b] != b) return -1;           // one endpoint already steckered
    plug[a] = b;
    plug[b] = a;
    return 0;
}

// --- Names -----------------------------------------------------------------------------

int enigma_rotor_from_name(const char *s) {
    for (int r = 0; r < ENIGMA_N_ROTORS; r++)
        if (strcasecmp(s, ROTOR_NAME[r]) == 0) return r;
    return -1;
}

int enigma_reflector_from_name(const char *s) {
    if (strcasecmp(s, "B") == 0) return ENIGMA_UKW_B;
    if (strcasecmp(s, "C") == 0) return ENIGMA_UKW_C;
    if (strcasecmp(s, "BTHIN") == 0 || strcasecmp(s, "B-THIN") == 0 ||
        strcasecmp(s, "BDUNN") == 0) return ENIGMA_UKW_B_THIN;
    if (strcasecmp(s, "CTHIN") == 0 || strcasecmp(s, "C-THIN") == 0 ||
        strcasecmp(s, "CDUNN") == 0) return ENIGMA_UKW_C_THIN;
    return -1;
}

const char *enigma_rotor_name(int rotor_id) {
    if (rotor_id < 0 || rotor_id >= ENIGMA_N_ROTORS) return "?";
    return ROTOR_NAME[rotor_id];
}

const char *enigma_reflector_name(int reflector_id) {
    if (reflector_id < 0 || reflector_id >= ENIGMA_N_REFLECTORS) return "?";
    return REFLECTOR_NAME[reflector_id];
}

int enigma_num_notches(int rotor_id) {
    if (rotor_id < 0 || rotor_id >= ENIGMA_N_ROTORS) return 0;
    return (int) strlen(ROTOR_NOTCH[rotor_id]);
}

const char *enigma_rotor_wiring(int rotor_id) {
    if (rotor_id < 0 || rotor_id >= ENIGMA_N_ROTORS) return "";
    return ROTOR_WIRING[rotor_id];
}

const char *enigma_reflector_wiring(int reflector_id) {
    if (reflector_id < 0 || reflector_id >= ENIGMA_N_REFLECTORS) return "";
    return REFLECTOR_WIRING[reflector_id];
}

// --- CLI parse helpers -----------------------------------------------------------------

#include <ctype.h>
#include <stdlib.h>

int enigma_parse_setting(const char *tok) {
    if (!tok || !*tok) return -1;
    if (isalpha((unsigned char) tok[0]) && tok[1] == '\0')
        return toupper((unsigned char) tok[0]) - 'A';
    // 1-based number 1..26
    char *end;
    long v = strtol(tok, &end, 10);
    if (*end != '\0' || v < 1 || v > 26) return -1;
    return (int) (v - 1);
}

int enigma_parse_settings(const char *s, int out[], int max) {
    char buf[256];
    strncpy(buf, s, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    int n = 0;
    for (char *tok = strtok(buf, " ,\t"); tok && n < max; tok = strtok(NULL, " ,\t")) {
        int v = enigma_parse_setting(tok);
        if (v < 0) return n;    // stop at the first bad token
        out[n++] = v;
    }
    return n;
}

int enigma_parse_rotors(const char *s, int rotor_ids[], int max) {
    char buf[256];
    strncpy(buf, s, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    int n = 0;
    for (char *tok = strtok(buf, " ,\t"); tok && n < max; tok = strtok(NULL, " ,\t")) {
        int id = enigma_rotor_from_name(tok);
        if (id < 0) return -1;
        rotor_ids[n++] = id;
    }
    return n;
}

int enigma_parse_plugs(const char *s, int plug[26]) {
    enigma_plug_identity(plug);
    char buf[256];
    strncpy(buf, s, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    int npairs = 0;
    for (char *tok = strtok(buf, " ,\t"); tok; tok = strtok(NULL, " ,\t")) {
        if (!isalpha((unsigned char) tok[0]) || !isalpha((unsigned char) tok[1]) || tok[2] != '\0')
            return -1;
        int a = toupper((unsigned char) tok[0]) - 'A';
        int b = toupper((unsigned char) tok[1]) - 'A';
        if (enigma_plug_set_pair(plug, a, b) != 0) return -1;   // conflict / self-pair
        npairs++;
    }
    return npairs;
}
