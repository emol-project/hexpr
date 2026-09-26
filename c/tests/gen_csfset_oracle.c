/* gen_csfset_oracle.c -- generate the frozen CSFSET pair-eval oracle.
 *
 * Per docs/SPEC_CSFSET_ORACLE.md. ADD-ONLY tooling: builds each CSFSET
 * sub-DRT via the public hexpr_drt_from_csfset and writes the current
 * (cb8a3f3) sub-DRT path's subspace_eval output to reference_csfset/<name>.txt
 * for which = full and (where applicable) one. The committed files ARE the
 * oracle; this generator is for transparency / regeneration only.
 *
 * Usage:  ./gen_csfset_oracle [outdir]
 *   outdir defaults to ../../reference_csfset (the committed location).
 *
 * Builds the source CSFs from the committed reference cases by row index, so
 * the oracle stays anchored to s_full_4o4e / s_foci_mid / s_soci_mid.
 */
#include "hexpr.h"
#include "csfset_oracle_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ------------------------------------------------------------------ */
/* Oracle case table                                                   */
/* ------------------------------------------------------------------ */

#define MAX_ROWS 32

typedef struct {
    const char *name;
    /* Source DRT the rows are drawn from (NULL => literal rows below). */
    hexpr_drt_kind_t kind;
    int nve, spin, ncore, nval, next;
    /* Row selection. If src_norb==0 the rows are taken from the source DRT's
     * canonical CSF list by index (rows[]); otherwise rows[] is interpreted
     * as packed step bytes (src_norb each) given literally. */
    int src_norb;             /* 0 => use source DRT rows by index */
    int rows[MAX_ROWS];       /* source-row indices, or literal step bytes */
    int nrows;                /* number of source rows (or literal CSFs) */
    int do_one;               /* also emit which=one */
} ocase_t;

/* Literal-row helper cases (norb 4). */
static const ocase_t CASES[] = {
    /* O1: all 20 CSFs of s_full_4o4e in canonical order; full + one. */
    { "O1", HEXPR_DRT_SOCI, 4, 0, 0, 4, 0, 0,
      {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}, 20, 1 },

    /* O2: [uefd, eudf] (norb 4) -- non-canonical order + superset M>N. */
    { "O2", HEXPR_DRT_CSFSET, 0,0,0,0,0, 4,
      { 1,0,3,2,  0,1,2,3 }, 2, 1 },

    /* O5: [eeff, eeff] (norb 4) -- duplicate input (N>M dedup). */
    { "O5", HEXPR_DRT_CSFSET, 0,0,0,0,0, 4,
      { 0,0,3,3,  0,0,3,3 }, 2, 0 },

    /* OF: representative subset of s_foci_mid (core2/act2/ext2, nve6, FOCI),
     * input order. Rows span core-hole-at-orb0 (uf..), core-hole-at-orb1
     * (fu..) and cores-full (ff..) sectors, with external particles. */
    { "OF", HEXPR_DRT_FOCI, 6, 0, 2, 2, 2, 0,
      { 0, 5, 13, 16, 21, 29, 32, 37, 40 }, 9, 1 },

    /* OFs: the SAME OF rows in a deliberately non-canonical (shuffled) order;
     * input-order contract check (OFs == OF reordered). */
    { "OFs", HEXPR_DRT_FOCI, 6, 0, 2, 2, 2, 0,
      { 40, 16, 5, 32, 0, 29, 13, 37, 21 }, 9, 0 },

    /* OS: representative subset of s_soci_mid (core2/act2/ext2, nve6, SOCI),
     * input order. Includes the double-core-hole sector (uu.., ud..). */
    { "OS", HEXPR_DRT_SOCI, 6, 0, 2, 2, 2, 0,
      { 0, 3, 15, 16, 26, 29, 59, 92, 101 }, 9, 1 },
};
static const int NCASES = (int)(sizeof(CASES) / sizeof(CASES[0]));

/* ------------------------------------------------------------------ */
/* Sortable term rows                                                  */
/* ------------------------------------------------------------------ */

typedef struct { int bp, kp, ij, pqrs, seq; double coef; } term_t;

static int term_cmp(const void *a, const void *b)
{
    const term_t *x = a, *y = b;
    if (x->bp   != y->bp)   return x->bp   - y->bp;
    if (x->kp   != y->kp)   return x->kp   - y->kp;
    if (x->ij   != y->ij)   return x->ij   - y->ij;
    if (x->pqrs != y->pqrs) return x->pqrs - y->pqrs;
    return x->seq - y->seq;   /* stable tiebreak preserves emit order */
}

/* Build the packed input rows for a case; returns malloc'd buffer or NULL.
 * Sets *out_norb and copies the n_csfs * norb step bytes. */
static unsigned char *build_input_rows(const ocase_t *c, int *out_norb,
                                       int *out_ncsfs)
{
    if (c->src_norb > 0) {
        /* literal step bytes */
        int norb = c->src_norb;
        int n = c->nrows;
        unsigned char *buf = malloc((size_t)n * (size_t)norb);
        for (int i = 0; i < n * norb; i++)
            buf[i] = (unsigned char)c->rows[i];
        *out_norb = norb; *out_ncsfs = n;
        return buf;
    }
    /* draw rows from a source DRT by canonical index */
    hexpr_drt_t *src = hexpr_drt_create(c->kind, c->nve, c->spin,
                                        c->ncore, c->nval, c->next);
    if (!src) { fprintf(stderr, "  %s: source DRT create failed\n", c->name); return NULL; }
    int norb  = hexpr_drt_norb(src);
    int nsrc  = hexpr_drt_ncsfs(src);
    const unsigned char *sc = hexpr_drt_csfs(src);
    int n = c->nrows;
    unsigned char *buf = malloc((size_t)n * (size_t)norb);
    for (int i = 0; i < n; i++) {
        int r = c->rows[i];
        if (r < 0 || r >= nsrc) {
            fprintf(stderr, "  %s: row %d out of range (nsrc=%d)\n", c->name, r, nsrc);
            free(buf); hexpr_drt_destroy(src); return NULL;
        }
        memcpy(buf + (size_t)i * (size_t)norb,
               sc  + (size_t)r * (size_t)norb, (size_t)norb);
    }
    *out_norb = norb; *out_ncsfs = n;
    hexpr_drt_destroy(src);
    return buf;
}

static int emit_which(FILE *f, hexpr_drt_t *sub, const int *idx, int n,
                      hexpr_expr_which_t which, const char *wname)
{
    int cnt = 0;
    if (hexpr_subspace_count(sub, idx, n, which, &cnt) != HEXPR_OK) return 1;
    int32_t *bp = malloc((size_t)(cnt+1) * sizeof(int32_t));
    int32_t *kp = malloc((size_t)(cnt+1) * sizeof(int32_t));
    int32_t *ij = malloc((size_t)(cnt+1) * sizeof(int32_t));
    int32_t *pq = malloc((size_t)(cnt+1) * sizeof(int32_t));
    double  *co = malloc((size_t)(cnt+1) * sizeof(double));
    if (hexpr_subspace_eval(sub, idx, n, which, &cnt, bp, kp, ij, pq, co)
        != HEXPR_OK) {
        free(bp); free(kp); free(ij); free(pq); free(co); return 1;
    }
    term_t *T = malloc((size_t)(cnt+1) * sizeof(term_t));
    for (int k = 0; k < cnt; k++)
        T[k] = (term_t){ bp[k], kp[k], ij[k], pq[k], k, co[k] };
    qsort(T, (size_t)cnt, sizeof(term_t), term_cmp);

    fprintf(f, "which %s\n", wname);
    fprintf(f, "nterms %d\n", cnt);
    fprintf(f, "# bra_pos ket_pos ij pqrs coef\n");
    for (int k = 0; k < cnt; k++)
        fprintf(f, "%d %d %d %d %+.16e\n",
                T[k].bp, T[k].kp, T[k].ij, T[k].pqrs, T[k].coef);

    free(T); free(bp); free(kp); free(ij); free(pq); free(co);
    return 0;
}

static int gen_case(const ocase_t *c, const char *outdir)
{
    int norb = 0, n = 0;
    unsigned char *buf = build_input_rows(c, &norb, &n);
    if (!buf) return 1;

    hexpr_drt_t *sub = hexpr_drt_from_csfset(norb, buf, n);
    if (!sub) {
        fprintf(stderr, "  %s: from_csfset failed: %s\n", c->name, hexpr_last_error());
        free(buf); return 1;
    }

    /* The eval API takes INPUT-ORDER indices, so this is 0..n-1 and
     * bra_pos/ket_pos come back as input positions. Earlier versions mapped
     * each row to its canonical index here; the library does that now, and
     * the emitted files are unchanged by the move. */
    int *idx = malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) idx[i] = i;

    char path[1024];
    snprintf(path, sizeof(path), "%s/%s.txt", outdir, c->name);
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); free(idx); free(buf); hexpr_drt_destroy(sub); return 1; }

    fprintf(f, "# CSFSET pair-eval oracle -- frozen from the sub-DRT path (cb8a3f3).\n");
    fprintf(f, "# Generated by c/tests/gen_csfset_oracle.c. Committed files are the\n");
    fprintf(f, "# oracle; see reference_csfset/README. case %s\n", c->name);
    fprintf(f, "norb %d\n", norb);
    fprintf(f, "ncsfs %d\n", n);
    for (int i = 0; i < n; i++) {
        fprintf(f, "csf %d ", i);
        for (int o = 0; o < norb; o++)
            fprintf(f, " %d", buf[(size_t)i * (size_t)norb + o]);
        fprintf(f, "\n");
    }

    int rc = emit_which(f, sub, idx, n, HEXPR_EXPR_FULL, "full");
    if (!rc && c->do_one)
        rc = emit_which(f, sub, idx, n, HEXPR_EXPR_ONE, "one");

    fclose(f);
    free(idx); free(buf); hexpr_drt_destroy(sub);
    if (!rc) fprintf(stderr, "  wrote %s\n", path);
    return rc;
}

int main(int argc, char **argv)
{
    const char *outdir = (argc > 1) ? argv[1] : CSFSET_ORACLE_DIR;
    if (hexpr_init() != HEXPR_OK) {
        fprintf(stderr, "hexpr_init failed\n"); return 1;
    }
    int rc = 0;
    for (int i = 0; i < NCASES; i++)
        rc |= gen_case(&CASES[i], outdir);
    hexpr_shutdown();
    if (rc) fprintf(stderr, "gen_csfset_oracle: FAILED\n");
    else    fprintf(stderr, "gen_csfset_oracle: OK (%d cases)\n", NCASES);
    return rc;
}
