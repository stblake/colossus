//
//  Primitive known-answer tests for the Key Phrase cipher (type 86). Framework-free;
//  built by `make test` (links utils.c + keyphrase.c). A 26-letter phrase IS the cipher
//  alphabet: encryption maps plaintext letter i -> key[i]. With a phrase that is a
//  permutation of A..Z the map is a bijection (invertible); a phrase with repeats is
//  many-to-one (the ambiguity the solver must resolve by context).
//

#include "colossus.h"
#include "keyphrase.h"

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

int main(void) {
    init_alphabet(NULL);                            // 26-letter default; populates g_char_to_idx

    // A permutation phrase (the QWERTY keyboard = 26 distinct letters).
    int key[KEYPHRASE_LEN];
    int n = keyphrase_build_key("QWERTYUIOPASDFGHJKLZXCVBNM", key);
    CHECK(n == KEYPHRASE_LEN, "keyphrase key length %d != %d", n, KEYPHRASE_LEN);

    // Every key letter is distinct (a permutation), so the map is a bijection.
    int seen[26]; for (int i = 0; i < 26; i++) seen[i] = 0;
    for (int i = 0; i < KEYPHRASE_LEN; i++) seen[key[i]]++;
    int perm = 1; for (int i = 0; i < 26; i++) if (seen[i] != 1) perm = 0;
    CHECK(perm, "keyphrase QWERTY phrase should be a permutation of A..Z");

    // Encryption maps plaintext i -> key[i]; check the whole alphabet + a KAT (A -> Q).
    int pt[26], ct[26]; for (int i = 0; i < 26; i++) pt[i] = i;
    keyphrase_encrypt(pt, 26, key, ct);
    int ok = 1; for (int i = 0; i < 26; i++) if (ct[i] != key[i]) ok = 0;
    CHECK(ok, "keyphrase_encrypt must map i -> key[i]");
    CHECK(key[0] == g_char_to_idx[(unsigned char) 'Q'], "keyphrase KAT: A must map to Q");

    // A phrase with repeated letters yields a many-to-one (non-injective) key.
    int key2[KEYPHRASE_LEN];
    int n2 = keyphrase_build_key("HELLOWORLDXYZABCFGIJKMNPQS", key2);
    CHECK(n2 == KEYPHRASE_LEN, "keyphrase key2 length %d != %d", n2, KEYPHRASE_LEN);
    int distinct = 0, mark[26]; for (int i = 0; i < 26; i++) mark[i] = 0;
    for (int i = 0; i < KEYPHRASE_LEN; i++) if (!mark[key2[i]]++) distinct++;
    CHECK(distinct < KEYPHRASE_LEN, "keyphrase: a repeating phrase must be many-to-one");

    printf("test_keyphrase: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
