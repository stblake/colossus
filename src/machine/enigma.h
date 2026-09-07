// =====================================================================
//  Enigma machine primitive (TYPE enigma)
// =====================================================================
//
// The German Enigma: a rotor stream cipher. Signal path per keypress is
//   plugboard -> (right..left rotors) -> reflector -> (left..right rotors) -> plugboard.
// It is SELF-RECIPROCAL (encrypt == decrypt from the same start key) and has NO fixed
// point (the reflector guarantees no letter enciphers to itself).
//
// This module models the Services Enigma I (rotors I-V, reflectors UKW-B/C) and the naval
// M3/M4 (rotors VI-VIII with two notches, the non-stepping Greek 4th rotor Beta/Gamma, and
// the thin reflectors B/C). ETW (entry wheel) is the identity. All wiring/notch strings are
// the standard reference values (Wikipedia "Enigma rotor details" / Crypto Museum) and are
// pinned as KAT vectors in tests/test_enigma.c.
//
// Convention: letters are 0..25 (A=0). Ring (Ringstellung) and window position
// (Grundstellung / message key) are stored 0-based here; the CLI/report convert to the
// 1-based (A=1) or letter forms the historical literature uses.

#ifndef ENIGMA_H
#define ENIGMA_H

#include <stdbool.h>

// Rotor ids. 0..7 = I..VIII; 8/9 = Greek Beta/Gamma (M4 4th wheel, non-stepping).
enum {
    ENIGMA_I = 0, ENIGMA_II, ENIGMA_III, ENIGMA_IV, ENIGMA_V,
    ENIGMA_VI, ENIGMA_VII, ENIGMA_VIII,
    ENIGMA_BETA, ENIGMA_GAMMA,
    ENIGMA_N_ROTORS
};

// Reflector ids.
enum {
    ENIGMA_UKW_B = 0, ENIGMA_UKW_C,
    ENIGMA_UKW_B_THIN, ENIGMA_UKW_C_THIN,
    ENIGMA_N_REFLECTORS
};

#define ENIGMA_MAX_WHEELS 4   // M4: Greek + 3 stepping rotors

// A full Enigma key: machine geometry + settings. `n_wheels` is 3 (M3/Enigma I) or 4 (M4).
// rotor[0] is the LEFTMOST wheel; for M4 that is the non-stepping Greek rotor and the three
// stepping rotors are rotor[1..3]. ring/pos are 0..25 per wheel; plug is a 26-entry
// involution (plug[i]==i means letter i is unsteckered).
typedef struct {
    int reflector;                 // ENIGMA_UKW_*
    int n_wheels;                  // 3 or 4
    int rotor[ENIGMA_MAX_WHEELS];  // rotor ids, [0]=leftmost
    int ring[ENIGMA_MAX_WHEELS];   // 0..25
    int pos[ENIGMA_MAX_WHEELS];    // 0..25 (window position; advances during encrypt)
    int plug[26];                  // involution over 0..25
} EnigmaKey;

// Build the internal wiring/notch tables from the reference strings. Idempotent; call once
// on the main thread before any encrypt (solve_enigma, the generator, and tests all do).
void enigma_init(void);

// Encipher `in[0..n-1]` (0..25) into `out[0..n-1]`. Steps the rotors BEFORE each letter.
// Does not mutate `*key`. Decryption is the same call with the same start key.
void enigma_encrypt(const int *in, int n, const EnigmaKey *key, int *out);

// Advance the rotor odometer one step (double-step anomaly + multi-notch aware). Exposed
// for tests; enigma_encrypt calls it internally.
void enigma_step(EnigmaKey *key);

// Encipher a single letter through the machine in its CURRENT position (no stepping).
int  enigma_encipher_letter(const EnigmaKey *key, int c);

// Fill perm[0..25] with the full scrambler permutation (plugboard EXCLUDED) at menu
// position `offset`: the machine `*base` stepped (offset+1) times, then each of the 26
// letters enciphered through the rotors+reflector only. Used by the Bombe. perm is an
// involution with no fixed point.
void enigma_scrambler_at(const EnigmaKey *base, int offset, int perm[26]);

// Plugboard helpers.
void enigma_plug_identity(int plug[26]);
int  enigma_plug_set_pair(int plug[26], int a, int b);  // 0 ok; -1 if a==b or a conflict

// Name <-> id helpers for the CLI/generator. Rotor names: "I".."VIII","BETA","GAMMA"
// (case-insensitive; also accepts "6".."8" style is NOT supported — use roman). Reflector
// names: "B","C","BTHIN"/"B-THIN","CTHIN"/"C-THIN".
int  enigma_rotor_from_name(const char *s);       // -1 on failure
int  enigma_reflector_from_name(const char *s);   // -1 on failure
const char *enigma_rotor_name(int rotor_id);
const char *enigma_reflector_name(int reflector_id);

// Number of turnover notches on a rotor (0 for Greek, 1 for I-V, 2 for VI-VIII).
int  enigma_num_notches(int rotor_id);

// CLI parse helpers (used by main()). A ring/position setting token is a single letter
// A..Z (A=0) or a 1-based number 1..26. Lists are space/comma-separated.
int  enigma_parse_setting(const char *tok);                    // 0..25, or -1
int  enigma_parse_settings(const char *s, int out[], int max); // count parsed (each 0..25)
int  enigma_parse_rotors(const char *s, int rotor_ids[], int max); // count, or -1 on bad name
int  enigma_parse_plugs(const char *s, int plug[26]);          // pair count, or -1

// Raw reference strings, exposed so tests can KAT-pin them. Rotor wiring is the A-> mapping
// (26 chars); reflector wiring is the involution (26 chars).
const char *enigma_rotor_wiring(int rotor_id);
const char *enigma_reflector_wiring(int reflector_id);

#endif // ENIGMA_H
