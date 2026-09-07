// Standalone Enigma test-data generator. Links the real enigma.c + utils.c so the generator
// and the solver can never drift in convention.
//
//   make enigma_gen
//   ./tools/enigma_gen plaintext.txt B II,I,III A,A,A A,B,C "EZ RW MV" >cipher.txt 2>solution.txt
//
// argv: <plaintext|-> <reflector> <rotors> <rings> <startpos> [plugs]
//   reflector : B | C | Bthin | Cthin
//   rotors    : comma/space list of 3 (M3) or 4 (M4, Greek first) names, e.g. II,I,III
//   rings     : comma/space list of settings (letters A..Z or 1-based numbers), one per wheel
//   startpos  : comma/space list of window positions, one per wheel
//   plugs     : optional pairs, e.g. "EZ RW MV"  (quote if space-separated)
// stdout: the ciphertext (letters only). stderr: the folded plaintext (the solution).

#include <stdio.h>
#include <string.h>
#include "colossus.h"
#include "enigma.h"

#define MAXLEN 20000

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s <plaintext|-> <reflector> <rotors> <rings> <startpos> [plugs]\n",
                argv[0]);
        return 1;
    }
    init_alphabet(NULL);
    enigma_init();

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s\n", argv[1]); return 1; }
    static int raw[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n' && n < MAXLEN) {
        if (ch >= 'A' && ch <= 'Z') raw[n++] = ch - 'A';
        else if (ch >= 'a' && ch <= 'z') raw[n++] = ch - 'a';
    }
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "ERROR: no letters in plaintext\n"); return 1; }

    EnigmaKey k;
    memset(&k, 0, sizeof(k));
    k.reflector = enigma_reflector_from_name(argv[2]);
    if (k.reflector < 0) { fprintf(stderr, "ERROR: bad reflector \"%s\"\n", argv[2]); return 1; }

    int ids[4];
    int nr = enigma_parse_rotors(argv[3], ids, 4);
    if (nr != 3 && nr != 4) { fprintf(stderr, "ERROR: -rotors needs 3 or 4 names\n"); return 1; }
    k.n_wheels = nr;
    for (int i = 0; i < nr; i++) k.rotor[i] = ids[i];

    int rings[4], pos[4];
    if (enigma_parse_settings(argv[4], rings, 4) != nr ||
        enigma_parse_settings(argv[5], pos, 4) != nr) {
        fprintf(stderr, "ERROR: rings/startpos must each list %d settings\n", nr);
        return 1;
    }
    for (int i = 0; i < nr; i++) { k.ring[i] = rings[i]; k.pos[i] = pos[i]; }

    enigma_plug_identity(k.plug);
    if (argc >= 7 && enigma_parse_plugs(argv[6], k.plug) < 0) {
        fprintf(stderr, "ERROR: bad plugboard \"%s\"\n", argv[6]);
        return 1;
    }

    static int ct[MAXLEN];
    enigma_encrypt(raw, n, &k, ct);
    for (int i = 0; i < n; i++) putchar(index_to_char(ct[i]));
    putchar('\n');
    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}
