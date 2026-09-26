/* loopgen/lk_table.h -- two-phase L(k) table (POD machinery).
 *
 * Port of RubyRef/emol_expression.rb `module LkTable` (:244-292) plus its
 * two-phase driver (:831-835) and warm_lk_table (:853-860). L(k) =
 * hmlkva(gf,ff,hf,psi,phi,theta,theta_prev) is a PURE function of its 7 integer
 * arguments, so we compute each distinct value once (generate pass) and read it
 * on the hot line (search pass).
 *
 *   TWO EXPLICIT PASSES the reader can see (never compute-on-miss as the primary
 *   search path):
 *     pass 1 (generate): warm_lk_table walks the same DRT ops `all` walks; the
 *                        hot line calls lk_value in GENERATE phase, which computes
 *                        L(k) via hmlkva->matrix_el+wigner_9j and fills the table.
 *                        This is the ONLY phase hmlkva runs on the whole path.
 *     pass 2 (search):   `all` walks again; the hot line is a pure table lookup.
 *                        A miss is a counted safety net (assert 0), not the path.
 *
 * This header holds the pure-POD machinery only (lk_table_t, lk_pack, lk_get,
 * lk_put, lk_clear, the phase enum). The parts that call hmlkva / iterate
 * BROOKS_CASES (lk_value, phase control, warm_lk_table, the search-phase
 * counters and their getters) live in expr_tree.c / expr.c and are wired there.
 */
#ifndef HEXPR_LK_TABLE_H
#define HEXPR_LK_TABLE_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* ----------------------------------------------------------------------
 * 1. KEY PACKING  (port of LkTable.pack, emol_expression.rb:252-268)
 *
 * Mixed radix, Ruby radices verbatim:
 *   R_GF=8, R_FF=8, R_HF=16, R_PSI=16, R_PHI=16, R_TH=8, R_THP=8
 * Total key space = 8*8*16*16*16*8*8 = 2^24 = 16,777,216.
 * Max packed key   = 2^24 - 1        = 16,777,215.
 * Key type uint32_t: fits with ~256x margin (UINT32_MAX = 4.29e9); even
 * int32_t would fit (~128x).  Real component ranges are far inside the radices
 *   gf,ff in 1..4 (<8), hf in 1..13 (<16), psi,phi = 2S+1 (<=7 here, <16),
 *   theta,theta_prev in 1..3 (<8)  -- so every component clears its radix.
 * Out-of-range component => HARD ERROR (abort), mirroring Ruby's raise; never
 * a silent alias.  With the margins above it cannot fire in practice; it is a
 * bug trap for a future caller that violates the range invariant.
 * -------------------------------------------------------------------- */
enum { LK_R_GF = 8, LK_R_FF = 8, LK_R_HF = 16, LK_R_PSI = 16,
       LK_R_PHI = 16, LK_R_TH = 8, LK_R_THP = 8 };
#define LK_KEY_SPACE (LK_R_GF*LK_R_FF*LK_R_HF*LK_R_PSI*LK_R_PHI*LK_R_TH*LK_R_THP)  /* 2^24 */
#define LK_KEY_EMPTY 0xFFFFFFFFu   /* sentinel: > any valid key (2^24-1) */

static inline uint32_t lk_pack(int gf, int ff, int hf, int psi, int phi,
                               int theta, int theta_prev)
{
    if (gf  < 0 || gf  >= LK_R_GF  || ff    < 0 || ff    >= LK_R_FF  ||
        hf  < 0 || hf  >= LK_R_HF  || psi   < 0 || psi   >= LK_R_PSI ||
        phi < 0 || phi >= LK_R_PHI || theta < 0 || theta >= LK_R_TH  ||
        theta_prev < 0 || theta_prev >= LK_R_THP) {
        fprintf(stderr, "lk_pack: component out of radix range: "
                "gf=%d ff=%d hf=%d psi=%d phi=%d theta=%d theta_prev=%d\n",
                gf, ff, hf, psi, phi, theta, theta_prev);
        abort();                                  /* hard error, not silent alias */
    }
    return (uint32_t)(((((((unsigned)gf * LK_R_FF + ff) * LK_R_HF + hf)
                        * LK_R_PSI + psi) * LK_R_PHI + phi) * LK_R_TH + theta)
                        * LK_R_THP + theta_prev);
}

/* ----------------------------------------------------------------------
 * 2. TABLE STORAGE  (open-addressing hash: packed uint32 key -> double)
 *
 * WHY open addressing (not direct-index, not sorted+bsearch):
 *   - direct-index over the 2^24 Ruby key space would be 2^24 doubles = 128 MB
 *     for ~500 live entries -- absurd density.
 *   - sorted array + bsearch needs a post-generate sort and cannot dedup /
 *     look up cheaply *during* the fill pass (generate hits its own keys).
 *   - open addressing gives O(1) explicit fill AND O(1) explicit lookup on the
 *     packed int, plain POD, no function pointers, tiny fixed footprint.
 * Measured live sizes (Ruby): whole-gen key set 458, pair 423, ~62 distinct
 * values; the two paths SHARE one table (keys are DRT-independent -- psi/phi are
 * spin multiplicities, not node ids), so the union is well under 1000.
 * LK_CAP = 4096 (power of two) => load factor < 0.25 even at 1000 keys.
 * 4096 * (4B key + 8B val) = 48 KB.  Overfill beyond 3/4 => hard error (asserts
 * the size model, same spirit as the range raise).
 * -------------------------------------------------------------------- */
#define LK_CAP       4096u          /* power of two */
#define LK_CAP_MASK  (LK_CAP - 1u)
#define LK_MAX_LOAD  (LK_CAP * 3u / 4u)

typedef struct {
    uint32_t keys[LK_CAP];   /* LK_KEY_EMPTY = free */
    double   vals[LK_CAP];
    int      size;
    /* phase + instrumentation (used by the wiring in expr_tree.c) */
    int      phase;          /* LK_OFF / LK_GENERATE / LK_SEARCH */
    long     miss;           /* search-phase table misses (assert 0) */
    long     search_wig_calls; /* wigner_9j calls during SEARCH (assert 0) */
    long     search_mel_calls; /* matrix_el  calls during SEARCH (assert 0) */
    long     off_fills;        /* OFF-phase (pair) fill-on-miss count (informational) */
} lk_table_t;

enum { LK_OFF = 0, LK_GENERATE = 1, LK_SEARCH = 2 };

/* Knuth multiplicative hash, take the high bits. */
static inline uint32_t lk_slot0(uint32_t key)
{
    return (uint32_t)((key * 2654435761u) >> (32 - 12));   /* 12 = log2(LK_CAP) */
}

/* Pure lookup. Returns 1 and sets *out on hit; returns 0 on miss. */
static inline int lk_get(const lk_table_t *t, uint32_t key, double *out)
{
    uint32_t i = lk_slot0(key);
    for (;;) {
        uint32_t k = t->keys[i];
        if (k == key)         { *out = t->vals[i]; return 1; }
        if (k == LK_KEY_EMPTY) return 0;             /* empty slot => not present */
        i = (i + 1) & LK_CAP_MASK;                   /* linear probe */
    }
}

/* Insert (key,val) if absent. Overfill beyond 3/4 is a hard error (mirrors the
 * range raise). No update path: L(k) is a pure fn, so a present key already
 * holds the identical value -- and we VERIFY that: if a present key carries a
 * different value, abort. That can only happen on a pack collision or a radix
 * error, so it is a cheap correctness sentinel, never expected to fire. */
static inline void lk_put(lk_table_t *t, uint32_t key, double val)
{
    if ((uint32_t)t->size >= LK_MAX_LOAD) {
        fprintf(stderr, "lk_put: table over capacity (%d/%u); size model wrong\n",
                t->size, LK_CAP);
        abort();
    }
    uint32_t i = lk_slot0(key);
    for (;;) {
        uint32_t k = t->keys[i];
        if (k == key) {                              /* already present */
            if (t->vals[i] != val) {                 /* pure fn => must be identical */
                fprintf(stderr, "lk_put: key %u present with %.17g but incoming %.17g "
                        "(pack collision or radix error)\n", key, t->vals[i], val);
                abort();
            }
            return;
        }
        if (k == LK_KEY_EMPTY) { t->keys[i] = key; t->vals[i] = val; t->size++; return; }
        i = (i + 1) & LK_CAP_MASK;
    }
}

/* Reset all slots to empty (called once before first generate, if desired). The
 * table is otherwise a persistent process-lifetime cache: a warmed key stays
 * warm across builds, which is safe because L(k) is a pure function. */
static inline void lk_clear(lk_table_t *t)
{
    for (uint32_t i = 0; i < LK_CAP; i++) t->keys[i] = LK_KEY_EMPTY;
    t->size = 0; t->phase = LK_OFF; t->miss = 0;
    t->search_wig_calls = 0; t->search_mel_calls = 0; t->off_fills = 0;
}

#endif  /* HEXPR_LK_TABLE_H */
