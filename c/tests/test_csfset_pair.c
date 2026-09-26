/* test_csfset_pair.c -- consumer for the frozen CSFSET pair-eval oracle.
 *
 * Per docs/SPEC_CSFSET_ORACLE.md. For each committed reference_csfset/<name>.txt
 * it rebuilds the CSFSET sub-DRT via hexpr_drt_from_csfset, re-runs the current
 * library's subspace_eval (and pair_eval), and asserts the output matches the
 * committed file exactly: bra_pos/ket_pos and ij/pqrs integer-exact, coef bit-
 * exact (same code path now; %.16e round-trips a double losslessly). It also:
 *   - cross-checks O1 against reference_pair/s_full_4o4e (Hermitian first-wins);
 *   - confirms OFs is OF reordered (input-order contract: identical canonical
 *     (ij,pqrs,coef) multiset under the input-order remap).
 *
 * Now: confirms the oracle round-trips against current code. Later: this SAME
 * unchanged test becomes the guard that the deferred direct-per-pair rewrite
 * reproduces the frozen baseline.
 */
#include "hexpr.h"
#include "csfset_oracle_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define MAX_NORB 8
#define MAX_ROWS 32

static int failures = 0;
static int checks   = 0;
static void check(const char *label, int ok)
{
    checks++;
    printf("  %-58s %s\n", label, ok ? "PASS" : "FAIL");
    if (!ok) failures++;
}

/* ------------------------------------------------------------------ */
/* Parsed oracle file                                                  */
/* ------------------------------------------------------------------ */
typedef struct { int bp, kp, ij, pqrs, seq; double coef; } term_t;

typedef struct {
    char name[16];
    int  norb, ncsfs;
    unsigned char rows[MAX_ROWS * MAX_NORB];
    term_t *full; int nfull;
    term_t *one;  int none_; int has_one;
} ofile_t;

static int term_cmp(const void *a, const void *b)   /* (bp,kp,ij,pqrs,seq) */
{
    const term_t *x = a, *y = b;
    if (x->bp   != y->bp)   return x->bp   - y->bp;
    if (x->kp   != y->kp)   return x->kp   - y->kp;
    if (x->ij   != y->ij)   return x->ij   - y->ij;
    if (x->pqrs != y->pqrs) return x->pqrs - y->pqrs;
    return x->seq - y->seq;
}
static int term_cmp_canon(const void *a, const void *b) /* (ij,pqrs,coef) */
{
    const term_t *x = a, *y = b;
    if (x->ij   != y->ij)   return x->ij   - y->ij;
    if (x->pqrs != y->pqrs) return x->pqrs - y->pqrs;
    if (x->coef < y->coef)  return -1;
    if (x->coef > y->coef)  return  1;
    return 0;
}

/* parse one oracle file; returns 0 on success */
static int parse_ofile(const char *path, ofile_t *o)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return 1; }
    memset(o, 0, sizeof(*o));
    char line[512];
    enum { M_NONE, M_FULL, M_ONE } mode = M_NONE;
    int curn = 0, curseq = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        if (strncmp(line, "norb ", 5) == 0) { sscanf(line, "norb %d", &o->norb); continue; }
        if (strncmp(line, "ncsfs ", 6) == 0) { sscanf(line, "ncsfs %d", &o->ncsfs); continue; }
        if (strncmp(line, "csf ", 4) == 0) {
            /* csf <row_index> b0 b1 ... b{norb-1} */
            char *q = line + 4, *end;
            int ri = (int)strtol(q, &end, 10); q = end;
            for (int b = 0; b < o->norb; b++) {
                long v = strtol(q, &end, 10);
                if (end == q) break;
                o->rows[(size_t)ri * (size_t)o->norb + b] = (unsigned char)v;
                q = end;
            }
            continue;
        }
        if (strncmp(line, "which ", 6) == 0) {
            char w[16] = {0}; sscanf(line, "which %15s", w);
            if (strcmp(w, "full") == 0)      mode = M_FULL;
            else if (strcmp(w, "one") == 0) { mode = M_ONE; o->has_one = 1; }
            else                             mode = M_NONE;
            curn = 0; curseq = 0;
            continue;
        }
        if (strncmp(line, "nterms ", 7) == 0) {
            int nt = 0; sscanf(line, "nterms %d", &nt);
            term_t *buf = malloc((size_t)(nt > 0 ? nt : 1) * sizeof(term_t));
            if (mode == M_FULL)      { o->full = buf; o->nfull = 0; }
            else if (mode == M_ONE) { o->one  = buf; o->none_ = 0; }
            curn = 0;
            continue;
        }
        /* data row: bp kp ij pqrs coef */
        {
            term_t t; t.seq = curseq++;
            if (sscanf(line, "%d %d %d %d %lf",
                       &t.bp, &t.kp, &t.ij, &t.pqrs, &t.coef) == 5) {
                if (mode == M_FULL && o->full) { o->full[curn++] = t; o->nfull = curn; }
                else if (mode == M_ONE && o->one) { o->one[curn++] = t; o->none_ = curn; }
            }
        }
    }
    fclose(f);
    return 0;
}

/* run subspace_eval for a which, return sorted term array (caller frees) */
static term_t *run_subspace(hexpr_drt_t *sub, const int *idx, int n,
                            hexpr_expr_which_t which, int *out_n)
{
    int cnt = 0;
    if (hexpr_subspace_count(sub, idx, n, which, &cnt) != HEXPR_OK) { *out_n = -1; return NULL; }
    int32_t *bp = malloc((size_t)(cnt+1)*sizeof(int32_t));
    int32_t *kp = malloc((size_t)(cnt+1)*sizeof(int32_t));
    int32_t *ij = malloc((size_t)(cnt+1)*sizeof(int32_t));
    int32_t *pq = malloc((size_t)(cnt+1)*sizeof(int32_t));
    double  *co = malloc((size_t)(cnt+1)*sizeof(double));
    hexpr_subspace_eval(sub, idx, n, which, &cnt, bp, kp, ij, pq, co);
    term_t *T = malloc((size_t)(cnt+1)*sizeof(term_t));
    for (int k = 0; k < cnt; k++) T[k] = (term_t){bp[k],kp[k],ij[k],pq[k],k,co[k]};
    qsort(T, (size_t)cnt, sizeof(term_t), term_cmp);
    free(bp); free(kp); free(ij); free(pq); free(co);
    *out_n = cnt;
    return T;
}

/* compare two sorted term arrays exactly (ints exact, coef bit-exact) */
static int terms_equal(const term_t *a, int na, const term_t *b, int nb,
                       const char *label)
{
    if (na != nb) {
        printf("    %s: count %d != %d\n", label, na, nb);
        return 0;
    }
    for (int k = 0; k < na; k++) {
        if (a[k].bp != b[k].bp || a[k].kp != b[k].kp ||
            a[k].ij != b[k].ij || a[k].pqrs != b[k].pqrs ||
            a[k].coef != b[k].coef) {
            printf("    %s: row %d differs: "
                   "(%d %d %d %d %.17e) vs (%d %d %d %d %.17e)\n", label, k,
                   a[k].bp,a[k].kp,a[k].ij,a[k].pqrs,a[k].coef,
                   b[k].bp,b[k].kp,b[k].ij,b[k].pqrs,b[k].coef);
            return 0;
        }
    }
    return 1;
}

/* per-pair cross-check via hexpr_pair_eval ("pair" API coverage) */
static int check_pair_api(hexpr_drt_t *sub, const int *idx, int n,
                          const term_t *file, int nfile, hexpr_expr_which_t which)
{
    int ok = 1;
    for (int i = 0; i < n && ok; i++) {
        for (int j = 0; j < n && ok; j++) {
            int pc = 0;
            hexpr_pair_eval(sub, idx[i], idx[j], which, &pc, NULL, NULL, NULL);
            int32_t *ij = malloc((size_t)(pc+1)*sizeof(int32_t));
            int32_t *pq = malloc((size_t)(pc+1)*sizeof(int32_t));
            double  *co = malloc((size_t)(pc+1)*sizeof(double));
            hexpr_pair_eval(sub, idx[i], idx[j], which, &pc, ij, pq, co);
            term_t *P = malloc((size_t)(pc+1)*sizeof(term_t));
            for (int k = 0; k < pc; k++) P[k] = (term_t){i,j,ij[k],pq[k],k,co[k]};
            qsort(P, (size_t)pc, sizeof(term_t), term_cmp);
            /* gather file rows with bp==i,kp==j (file is globally sorted, so
             * these are contiguous and already in (ij,pqrs,seq) order) */
            int start = -1, cnt = 0;
            for (int k = 0; k < nfile; k++)
                if (file[k].bp == i && file[k].kp == j) {
                    if (start < 0) start = k;
                    cnt++;
                }
            if (cnt != pc) { ok = 0; printf("    pair(%d,%d): count %d != file %d\n", i, j, pc, cnt); }
            for (int k = 0; k < pc && ok; k++) {
                const term_t *fr = &file[start + k];
                if (fr->ij != P[k].ij || fr->pqrs != P[k].pqrs || fr->coef != P[k].coef)
                    { ok = 0; printf("    pair(%d,%d) row %d differs\n", i, j, k); }
            }
            free(ij); free(pq); free(co); free(P);
        }
    }
    return ok;
}

/* rebuild sub-DRT + index array from a parsed file.
 *
 * The eval API takes INPUT-ORDER indices, so the array is 0..ncsfs-1 and
 * bra_pos/ket_pos come back as those same input positions. This used to run
 * each row through csfset_canon_index() to reach the canonical superset; the
 * library does that translation now. */
static hexpr_drt_t *rebuild(const ofile_t *o, int *idx)
{
    hexpr_drt_t *sub = hexpr_drt_from_csfset(o->norb, o->rows, o->ncsfs);
    if (!sub) return NULL;
    for (int i = 0; i < o->ncsfs; i++)
        idx[i] = i;
    return sub;
}

/* ------------------------------------------------------------------ */
/* O1 cross-check vs reference_pair/s_full_4o4e (Hermitian first-wins) */
/* ------------------------------------------------------------------ */
typedef struct { int ij_key, pqrs, seq; double coef; } row_t;
typedef struct { row_t *d; int n, cap; } rvec_t;
static void rv_push(rvec_t *v, int ij, int pqrs, int seq, double c) {
    if (v->n == v->cap) { v->cap = v->cap ? v->cap*2 : 4096; v->d = realloc(v->d, v->cap*sizeof(row_t)); }
    v->d[v->n++] = (row_t){ij,pqrs,seq,c};
}
static int row_cmp(const void *a, const void *b) {
    const row_t *x=a,*y=b;
    if (x->ij_key!=y->ij_key) return x->ij_key-y->ij_key;
    if (x->pqrs!=y->pqrs) return x->pqrs-y->pqrs;
    return x->seq-y->seq;
}
static int sort_dedup(rvec_t *v) {
    if (!v->n) return 0;
    qsort(v->d, v->n, sizeof(row_t), row_cmp);
    int out=0;
    for (int i=0;i<v->n;) { v->d[out++]=v->d[i]; int j=i+1;
        while (j<v->n && v->d[j].ij_key==v->d[i].ij_key && v->d[j].pqrs==v->d[i].pqrs) j++;
        i=j; }
    v->n=out; return out;
}
static int tri(int a,int b){int mx=a>b?a:b,mn=a<b?a:b;return mx*(mx+1)/2+mn;}

static void o1_cross_check(hexpr_drt_t *sub, const int *idx, int n)
{
    /* parse reference_pair pair.txt */
    const char *rp = "../../reference_pair/s_full_4o4e/pair.txt";
    FILE *f = fopen(rp, "r");
    if (!f) { check("O1 vs reference_pair/s_full_4o4e (open)", 0); return; }
    char line[512];
    while (fgets(line,sizeof(line),f)) if (strstr(line,"Energy Expression")) break;
    fgets(line,sizeof(line),f);  /* header */
    rvec_t ref={0}; int I,J,pqrs; double c; int seq=0;
    while (fscanf(f,"%d %d %d %lf",&I,&J,&pqrs,&c)==4)
        rv_push(&ref, tri(I-1,J-1), pqrs, seq++, c);
    fclose(f);
    int refn = sort_dedup(&ref);

    /* aggregate O1 full subspace (emit order) with first-wins */
    int cnt=0; hexpr_subspace_count(sub, idx, n, HEXPR_EXPR_FULL, &cnt);
    int32_t *bp=malloc((size_t)cnt*4),*kp=malloc((size_t)cnt*4),
            *ij=malloc((size_t)cnt*4),*pq=malloc((size_t)cnt*4);
    double *co=malloc((size_t)cnt*sizeof(double));
    hexpr_subspace_eval(sub, idx, n, HEXPR_EXPR_FULL, &cnt, bp,kp,ij,pq,co);
    rvec_t got={0};
    for (int t=0;t<cnt;t++) rv_push(&got, ij[t], pq[t], t, co[t]);
    int gotn = sort_dedup(&got);
    free(bp);free(kp);free(ij);free(pq);free(co);

    int mism=0; double maxd=0; int ri=0,gi=0;
    while (ri<refn && gi<gotn) {
        int cmp = (ref.d[ri].ij_key!=got.d[gi].ij_key)
                  ? ref.d[ri].ij_key-got.d[gi].ij_key
                  : ref.d[ri].pqrs-got.d[gi].pqrs;
        if (cmp<0){mism++;ri++;} else if (cmp>0){mism++;gi++;}
        else { double d=fabs(ref.d[ri].coef-got.d[gi].coef); if(d>maxd)maxd=d; if(d>1e-10)mism++; ri++;gi++; }
    }
    mism += (refn-ri)+(gotn-gi);
    printf("    O1 unique keys: ref=%d got=%d  max|dcoef|=%.3e\n", refn, gotn, maxd);
    check("O1 vs reference_pair/s_full_4o4e (Hermitian first-wins)", mism==0 && refn==gotn);
    free(ref.d); free(got.d);
}

/* ------------------------------------------------------------------ */
int main(void)
{
    printf("=== test_csfset_pair (frozen CSFSET oracle) ===\n");
    if (hexpr_init() != HEXPR_OK) { printf("FAIL: hexpr_init\n"); return 1; }

    ofile_t parsed[CSFSET_ORACLE_NCASES];
    int have[CSFSET_ORACLE_NCASES]; memset(have,0,sizeof(have));

    for (int ci = 0; ci < CSFSET_ORACLE_NCASES; ci++) {
        const char *nm = CSFSET_ORACLE_CASES[ci];
        char path[512];
        snprintf(path,sizeof(path), "%s/%s.txt", CSFSET_ORACLE_DIR, nm);
        printf("\n--- %s ---\n", nm);
        ofile_t *o = &parsed[ci];
        if (parse_ofile(path, o) != 0) { check("parse file", 0); continue; }
        snprintf(o->name,sizeof(o->name),"%s",nm);
        have[ci] = 1;

        int idx[MAX_ROWS];
        hexpr_drt_t *sub = rebuild(o, idx);
        char lab[128];
        snprintf(lab,sizeof(lab),"%s: from_csfset rebuild", nm);
        if (!sub) { check(lab, 0); continue; }
        check(lab, 1);
        /* Every input row must be reachable through the eval API by its own
         * input position. A row that the canonical walk does not contain
         * would come back as HEXPR_ERR_RANGE here. */
        int bad_idx = 0;
        for (int i=0;i<o->ncsfs;i++) {
            int c0 = 0;
            if (hexpr_pair_count(sub, idx[i], idx[i], HEXPR_EXPR_FULL, &c0)
                != HEXPR_OK) bad_idx = 1;
        }
        snprintf(lab,sizeof(lab),"%s: every input row is addressable", nm);
        check(lab, !bad_idx);

        /* FULL */
        int gn=0; term_t *g = run_subspace(sub, idx, o->ncsfs, HEXPR_EXPR_FULL, &gn);
        snprintf(lab,sizeof(lab),"%s: subspace full matches committed file", nm);
        check(lab, terms_equal(g, gn, o->full, o->nfull, "full"));
        snprintf(lab,sizeof(lab),"%s: pair_eval full == file blocks", nm);
        check(lab, check_pair_api(sub, idx, o->ncsfs, o->full, o->nfull, HEXPR_EXPR_FULL));
        free(g);

        /* ONE (if present) */
        if (o->has_one) {
            int on=0; term_t *g1 = run_subspace(sub, idx, o->ncsfs, HEXPR_EXPR_ONE, &on);
            snprintf(lab,sizeof(lab),"%s: subspace one matches committed file", nm);
            check(lab, terms_equal(g1, on, o->one, o->none_, "one"));
            snprintf(lab,sizeof(lab),"%s: pair_eval one == file blocks", nm);
            check(lab, check_pair_api(sub, idx, o->ncsfs, o->one, o->none_, HEXPR_EXPR_ONE));
            free(g1);
        }

        /* O1 cross-check vs reference_pair */
        if (strcmp(nm,"O1")==0) o1_cross_check(sub, idx, o->ncsfs);

        hexpr_drt_destroy(sub);
    }

    /* OF vs OFs input-order contract: identical canonical (ij,pqrs,coef)
     * multiset (only bra_pos/ket_pos labels differ between the two orders). */
    printf("\n--- OF vs OFs input-order contract ---\n");
    int iOF=-1,iOFs=-1;
    for (int ci=0;ci<CSFSET_ORACLE_NCASES;ci++){
        if (strcmp(CSFSET_ORACLE_CASES[ci],"OF")==0) iOF=ci;
        if (strcmp(CSFSET_ORACLE_CASES[ci],"OFs")==0) iOFs=ci;
    }
    if (iOF>=0 && iOFs>=0 && have[iOF] && have[iOFs]) {
        ofile_t *a=&parsed[iOF], *b=&parsed[iOFs];
        check("OF/OFs same CSF count", a->ncsfs==b->ncsfs);
        /* same CSF SET (as multisets of rows) */
        int same_set = (a->ncsfs==b->ncsfs);
        for (int i=0;i<a->ncsfs && same_set;i++){
            int found=0;
            for (int j=0;j<b->ncsfs;j++)
                if (memcmp(a->rows+(size_t)i*a->norb, b->rows+(size_t)j*b->norb, a->norb)==0){found=1;break;}
            if(!found) same_set=0;
        }
        check("OFs is a reordering of OF (same CSF set)", same_set);
        /* OFs must NOT be byte-identical order to OF (genuinely shuffled) */
        int reordered = (memcmp(a->rows, b->rows, (size_t)a->ncsfs*a->norb)!=0);
        check("OFs order differs from OF (non-canonical)", reordered);
        /* identical canonical (ij,pqrs,coef) multiset */
        term_t *ca=malloc((size_t)a->nfull*sizeof(term_t));
        term_t *cb=malloc((size_t)b->nfull*sizeof(term_t));
        memcpy(ca,a->full,(size_t)a->nfull*sizeof(term_t));
        memcpy(cb,b->full,(size_t)b->nfull*sizeof(term_t));
        qsort(ca,(size_t)a->nfull,sizeof(term_t),term_cmp_canon);
        qsort(cb,(size_t)b->nfull,sizeof(term_t),term_cmp_canon);
        int eq = (a->nfull==b->nfull);
        for (int k=0;k<a->nfull && eq;k++)
            if (ca[k].ij!=cb[k].ij||ca[k].pqrs!=cb[k].pqrs||ca[k].coef!=cb[k].coef) eq=0;
        check("OF/OFs identical canonical (ij,pqrs,coef) multiset", eq);
        free(ca); free(cb);
    } else {
        check("OF/OFs both present", 0);
    }

    /* free parsed */
    for (int ci=0;ci<CSFSET_ORACLE_NCASES;ci++) if (have[ci]) {
        free(parsed[ci].full); if (parsed[ci].has_one) free(parsed[ci].one);
    }

    hexpr_shutdown();
    printf("\n=== test_csfset_pair: %d/%d checks passed ===\n", checks-failures, checks);
    if (failures==0) printf("ALL TESTS PASSED\n");
    else             printf("%d TEST(S) FAILED\n", failures);
    return failures ? 1 : 0;
}
