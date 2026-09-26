/* run_reference.c -- CLI driver mirroring RubyRef/run_reference.rb
 *
 * Usage:
 *   run_reference [-o <outdir>] [--pair] [--pair-outdir <pair_outdir>]
 *                 [case_name ...]
 *
 * With no case_name args all 12 cases are processed.
 * --pair enables per-pair output (pair.txt, pair_one.txt, pair_summary.txt).
 * Progress lines go to stderr; stdout stays clean.
 *
 * _POSIX_C_SOURCE must precede all headers for clock_gettime. */
/**
 * @file run_reference.c
 * @brief Internal CLI reference driver (mirrors @c RubyRef/run_reference.rb):
 *        a 14-entry case table, per-case output generators, and an arg-parsing
 *        @c main. It computes no physics -- it only drives the API and writes
 *        text files.
 *
 * @note INTERNAL driver, NOT a copy-paste template for public-API users. From
 *       this file's own includes/calls: it pulls in the INTERNAL headers
 *       @c loopgen/drt.h, @c loopgen/expr.h and @c fmt.h (not only the public
 *       @c include/hexpr.h), and it calls entry points ABSENT from the public
 *       @c hexpr.h -- @c hexpr_drt_free and @c hexpr_expr_free (the public API
 *       exposes @c hexpr_drt_destroy / @c hexpr_expr_destroy instead; the
 *       free/destroy two-layer split and its ownership contract are documented
 *       upstream in @c loopgen/drt.c and @c loopgen/expr.c -- read those, not this
 *       driver), plus @c hexpr_drt_mk_csf_list and @c hexpr_print_csf. A public
 *       consumer should use the public destroy entry points and headers.
 */
#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../include/hexpr.h"
#include "loopgen/drt.h"
#include "loopgen/expr.h"
#include "fmt.h"

/* ------------------------------------------------------------------ */
/* Case table (mirrors run_reference.rb CASES)                        */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *name;
    int nve, spin, ncore, nval, next;
    hexpr_drt_kind_t kind;
} case_def_t;

static const case_def_t ALL_CASES[] = {
    {"s_full_4o4e",   4, 0, 0, 4, 0, HEXPR_DRT_SOCI},
    {"t_full_4o4e",   4, 2, 0, 4, 0, HEXPR_DRT_SOCI},
    {"s_full_6o6e",   6, 0, 0, 6, 0, HEXPR_DRT_SOCI},
    {"t_full_6o6e",   6, 2, 0, 6, 0, HEXPR_DRT_SOCI},
    {"s_soci_small",  6, 0, 1, 4, 1, HEXPR_DRT_SOCI},
    {"s_foci_small",  6, 0, 1, 4, 1, HEXPR_DRT_FOCI},
    {"t_soci_small",  6, 2, 1, 4, 1, HEXPR_DRT_SOCI},
    {"t_foci_small",  6, 2, 1, 4, 1, HEXPR_DRT_FOCI},
    {"s_foci_mid",    6, 0, 2, 2, 2, HEXPR_DRT_FOCI},
    {"s_soci_mid",    6, 0, 2, 2, 2, HEXPR_DRT_SOCI},
    {"s_soci_med",   10, 0, 2, 4, 4, HEXPR_DRT_SOCI},
    {"s_foci_med",   10, 0, 2, 4, 4, HEXPR_DRT_FOCI},
    {"t_soci_med",   10, 2, 2, 4, 4, HEXPR_DRT_SOCI},
    {"t_foci_med",   10, 2, 2, 4, 4, HEXPR_DRT_FOCI},
};
static const int NCASES = (int)(sizeof(ALL_CASES) / sizeof(ALL_CASES[0]));

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
/* ------------------------------------------------------------------ */

/** @brief Return the @c st_size of @p path in bytes, or @c -1 if @c stat() fails.
 *  @param path filesystem path to stat. @return size in bytes, or -1 on error. */
static long path_size(const char *path)
{
    struct stat st;
    return (stat(path, &st) == 0) ? (long)st.st_size : -1L;
}

/**
 * @brief Monotonic wall-clock seconds elapsed since @p t0.
 * @param t0 a start timestamp the caller filled via @c clock_gettime(CLOCK_MONOTONIC).
 * @return elapsed seconds as a @c double.
 * @note Reads @c CLOCK_MONOTONIC, so it is unaffected by wall-clock (time-of-day)
 *       adjustments.
 */
static double wall_elapsed(const struct timespec *t0)
{
    struct timespec t1;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    return (double)(t1.tv_sec  - t0->tv_sec)
         + (double)(t1.tv_nsec - t0->tv_nsec) * 1e-9;
}

/* ------------------------------------------------------------------ */
/* Per-case generator                                                  */
/* ------------------------------------------------------------------ */

/**
 * @brief Generate the per-case @c .drt.txt / @c .expr.txt / @c .expr_one.txt
 *        outputs for one case into @p outdir.
 * @param c      the case definition (parameters + name) to run.
 * @param outdir destination directory for the three text files.
 * @return 0 on success, 1 on any allocation / I/O / build error.
 *
 * @note Lifetime order, as actually written in the code below: the DRT is
 *       created, the @c .drt.txt file (DRT table + CSF list) is written while the
 *       DRT is live, then @c hexpr_expr_build(drt) builds the expression FROM the
 *       still-live DRT. Only AFTER that build returns is @c hexpr_drt_free(drt)
 *       called -- the DRT is NOT freed before the expression is built. The
 *       expression then outlives the DRT and is released last with
 *       @c hexpr_expr_free, after the full and one-electron files are written.
 *       (For the destroy/free ownership contract, see @c loopgen/drt.c /
 *       @c loopgen/expr.c per the @ref run_reference.c "file note".)
 * @note CSF-list ownership: @c hexpr_drt_mk_csf_list returns a caller-owned
 *       @c char** -- each @c csfs[i] is freed in the print loop and the @c csfs
 *       array itself is freed afterwards.
 */
static int generate_case(const case_def_t *c, const char *outdir)
{
    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    fprintf(stderr, "[%s] nve=%d spin=%d ncore=%d nval=%d next=%d type=%s\n",
            c->name, c->nve, c->spin, c->ncore, c->nval, c->next,
            c->kind == HEXPR_DRT_SOCI ? "soci" : "foci");

    /* --- (1) DRT + CSF list --- */
    hexpr_drt_t *drt = hexpr_drt_create(c->kind, c->nve, c->spin,
                                        c->ncore, c->nval, c->next);
    if (!drt) {
        fprintf(stderr, "  ERROR: hexpr_drt_create failed\n");
        return 1;
    }

    char path[1024];
    snprintf(path, sizeof(path), "%s/%s.drt.txt", outdir, c->name);
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); hexpr_drt_free(drt); return 1; }

    hexpr_drt_show(drt, f);
    fprintf(f, "\n CSF List \n");

    char **csfs = NULL;
    int ncsf = hexpr_drt_mk_csf_list(drt, &csfs);
    int norb  = hexpr_drt_norb(drt);
    for (int i = 0; i < ncsf; i++) {
        hexpr_print_csf(f, csfs[i], norb);
        free(csfs[i]);
    }
    free(csfs);
    fclose(f);
    fprintf(stderr, "  wrote %s (%ld bytes)\n", path, path_size(path));

    /* --- (2) Full expression --- */
    /* Build the expression from the still-live DRT, THEN free the DRT (the build
     * has already consumed everything it needs from drt). */
    hexpr_expr_t *expr = hexpr_expr_build(drt);
    hexpr_drt_free(drt);
    if (!expr) {
        fprintf(stderr, "  ERROR: hexpr_expr_build failed\n");
        return 1;
    }

    snprintf(path, sizeof(path), "%s/%s.expr.txt", outdir, c->name);
    f = fopen(path, "w");
    if (!f) { perror(path); hexpr_expr_free(expr); return 1; }
    hexpr_expr_show_full(expr, f);
    fclose(f);
    fprintf(stderr, "  wrote %s (%ld bytes)\n", path, path_size(path));

    /* --- (3) 1-electron expression --- */
    snprintf(path, sizeof(path), "%s/%s.expr_one.txt", outdir, c->name);
    f = fopen(path, "w");
    if (!f) { perror(path); hexpr_expr_free(expr); return 1; }
    hexpr_expr_show_one(expr, f);
    fclose(f);
    fprintf(stderr, "  wrote %s (%ld bytes)\n", path, path_size(path));

    hexpr_expr_free(expr);

    fprintf(stderr, "  elapsed: %.3f sec\n", wall_elapsed(&t0));

    return 0;
}

/* ------------------------------------------------------------------ */
/* Per-case pair output generator                                      */
/* ------------------------------------------------------------------ */

/**
 * @brief Generate the per-pair @c pair.txt / @c pair_one.txt / @c pair_summary.txt
 *        outputs for one case under @c pair_outdir/<case>/.
 * @param c           the case definition to run.
 * @param pair_outdir parent directory; a per-case subdirectory is created in it.
 * @return 0 on success, 1 on any allocation / I/O / build error (like
 *         @ref generate_case).
 *
 * @note Grow-on-demand buffers: @c pqrs_buf / @c coef_buf start at capacity 65;
 *       whenever @c hexpr_pair_count reports more terms than the current cap they
 *       are @c realloc'd to @c needed*2 BEFORE the matching @c hexpr_pair_eval
 *       fills them.
 * @note On-disk 1-based indexing: @c bra and @c ket drive the loop and the API
 *       calls as 0-based indices, but are written to the files as @c bra+1 /
 *       @c ket+1 because the reference file format is 1-indexed. (Output
 *       bookkeeping only.)
 * @note Directory creation uses @c system("mkdir -p ...") -- a deliberate
 *       shell-out and the only process spawn in this file (a driver convenience).
 */
static int generate_pair_outputs(const case_def_t *c, const char *pair_outdir)
{
    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    fprintf(stderr, "[%s] pair output\n", c->name);

    hexpr_drt_t *drt = hexpr_drt_create(c->kind, c->nve, c->spin,
                                        c->ncore, c->nval, c->next);
    if (!drt) {
        fprintf(stderr, "  ERROR: hexpr_drt_create failed\n");
        return 1;
    }

    int ncsfs = hexpr_drt_ncsfs(drt);
    int norb  = hexpr_drt_norb(drt);

    /* mkdir -p pair_outdir/case */
    char dir[1024], path[2048];
    snprintf(dir, sizeof(dir), "%s/%s", pair_outdir, c->name);
    {
        char cmd[1280];
        snprintf(cmd, sizeof(cmd), "mkdir -p \"%s\"", dir);
        if (system(cmd) != 0) {
            fprintf(stderr, "ERROR: could not create directory %s\n", dir);
            hexpr_drt_free(drt);
            return 1;
        }
    }

    /* Open output files */
    snprintf(path, sizeof(path), "%s/pair.txt", dir);
    FILE *f_full = fopen(path, "w");
    if (!f_full) { perror(path); hexpr_drt_free(drt); return 1; }

    snprintf(path, sizeof(path), "%s/pair_one.txt", dir);
    FILE *f_one = fopen(path, "w");
    if (!f_one) { perror(path); fclose(f_full); hexpr_drt_free(drt); return 1; }

    snprintf(path, sizeof(path), "%s/pair_summary.txt", dir);
    FILE *f_sum = fopen(path, "w");
    if (!f_sum) { perror(path); fclose(f_full); fclose(f_one); hexpr_drt_free(drt); return 1; }

    /* CSF list header for pair.txt and pair_one.txt */
    char **csfs = NULL;
    int ncsf_chk = hexpr_drt_mk_csf_list(drt, &csfs);
    (void)ncsf_chk;

    fprintf(f_full, " CSF List \n");
    fprintf(f_one,  " CSF List \n");
    for (int i = 0; i < ncsfs; i++) {
        hexpr_print_csf(f_full, csfs[i], norb);
        hexpr_print_csf(f_one,  csfs[i], norb);
        free(csfs[i]);
    }
    free(csfs);

    fprintf(f_full, " Energy Expression \n");
    fprintf(f_full, "       I       J       pqrs     Coef \n");
    fprintf(f_one,  " Energy Expression \n");
    fprintf(f_one,  "       I       J       pqrs     Coef \n");

    /* Summary header */
    fprintf(f_sum, "       I       J     n_full      n_one \n");

    /* Pair evaluation buffers -- grow on demand */
    int buf_cap = 65;
    int32_t *pqrs_buf = malloc(buf_cap * sizeof(int32_t));
    double  *coef_buf = malloc(buf_cap * sizeof(double));
    if (!pqrs_buf || !coef_buf) {
        fprintf(stderr, "ERROR: OOM allocating pair buffers\n");
        fclose(f_full); fclose(f_one); fclose(f_sum);
        free(pqrs_buf); free(coef_buf);
        hexpr_drt_free(drt);
        return 1;
    }

    for (int bra = 0; bra < ncsfs; bra++) {
        for (int ket = 0; ket < ncsfs; ket++) {
            int needed, n_full, n_one;

            /* Full expression */
            hexpr_pair_count(drt, bra, ket, HEXPR_EXPR_FULL, &needed);
            if (needed > buf_cap) {
                buf_cap = needed * 2;
                pqrs_buf = realloc(pqrs_buf, buf_cap * sizeof(int32_t));
                coef_buf = realloc(coef_buf, buf_cap * sizeof(double));
            }
            hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_FULL,
                            &n_full, NULL, pqrs_buf, coef_buf);
            for (int k = 0; k < n_full; k++) {
                fprintf(f_full, " %7d %7d    %10d    %24.16e \n",
                        bra + 1, ket + 1, (int)pqrs_buf[k], coef_buf[k]);
            }

            /* One-electron expression */
            hexpr_pair_count(drt, bra, ket, HEXPR_EXPR_ONE, &needed);
            if (needed > buf_cap) {
                buf_cap = needed * 2;
                pqrs_buf = realloc(pqrs_buf, buf_cap * sizeof(int32_t));
                coef_buf = realloc(coef_buf, buf_cap * sizeof(double));
            }
            hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_ONE,
                            &n_one, NULL, pqrs_buf, coef_buf);
            for (int k = 0; k < n_one; k++) {
                fprintf(f_one, " %7d %7d    %10d    %24.16e \n",
                        bra + 1, ket + 1, (int)pqrs_buf[k], coef_buf[k]);
            }

            /* Summary row */
            fprintf(f_sum, " %7d %7d    %7d    %7d \n",
                    bra + 1, ket + 1, n_full, n_one);
        }
    }

    fclose(f_full);
    fclose(f_one);
    fclose(f_sum);
    free(pqrs_buf);
    free(coef_buf);
    hexpr_drt_free(drt);

    fprintf(stderr, "  pair output: %s (elapsed: %.3f sec)\n",
            dir, wall_elapsed(&t0));
    return 0;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

/**
 * @brief Entry point: parse args, optionally filter the case table, run the
 *        generators, and report the processed count.
 *
 * Usage: @c [-o outdir] @c [--pair] @c [--pair-outdir d] @c [case ...]. With no
 * case names, all table entries run; @c --pair additionally runs
 * @ref generate_pair_outputs for each selected case.
 * @param argc,argv standard CLI arguments.
 * @return 0 on success, 1 on error (directory creation, a generator failure, or
 *         a case filter that matched nothing).
 *
 * @note Progress goes to @c stderr so @c stdout stays clean (mirrors the @c .rb
 *       driver).
 * @note The case filter matches @c case_def_t.name exactly; if filter names were
 *       given but none matched, it prints the available-case list to @c stderr
 *       and returns 1.
 */
int main(int argc, char *argv[])
{
    hexpr_init();

    const char *outdir      = "./reference_c";
    const char *pair_outdir = "./reference_pair_c";
    int pair_enabled = 0;
    const char *filter[64];
    int nfilter = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            outdir = argv[++i];
        } else if (strcmp(argv[i], "--pair") == 0) {
            pair_enabled = 1;
        } else if (strcmp(argv[i], "--pair-outdir") == 0 && i + 1 < argc) {
            pair_outdir = argv[++i];
        } else if (argv[i][0] != '-') {
            if (nfilter < 64)
                filter[nfilter++] = argv[i];
        }
    }

    /* mkdir -p outdir */
    {
        char cmd[1280];
        snprintf(cmd, sizeof(cmd), "mkdir -p \"%s\"", outdir);
        if (system(cmd) != 0) {
            fprintf(stderr, "ERROR: could not create output directory %s\n", outdir);
            return 1;
        }
    }

    int processed = 0;
    for (int i = 0; i < NCASES; i++) {
        const case_def_t *c = &ALL_CASES[i];

        /* apply filter */
        if (nfilter > 0) {
            int match = 0;
            for (int j = 0; j < nfilter; j++)
                if (strcmp(c->name, filter[j]) == 0) { match = 1; break; }
            if (!match) continue;
        }

        if (generate_case(c, outdir) != 0)
            return 1;
        processed++;
    }

    if (processed == 0 && nfilter > 0) {
        fprintf(stderr, "No case matched. Available cases:\n");
        for (int i = 0; i < NCASES; i++)
            fprintf(stderr, "  %s\n", ALL_CASES[i].name);
        return 1;
    }

    if (pair_enabled) {
        /* mkdir -p pair_outdir */
        {
            char cmd[1280];
            snprintf(cmd, sizeof(cmd), "mkdir -p \"%s\"", pair_outdir);
            if (system(cmd) != 0) {
                fprintf(stderr, "ERROR: could not create pair output directory %s\n",
                        pair_outdir);
                return 1;
            }
        }

        for (int i = 0; i < NCASES; i++) {
            const case_def_t *c = &ALL_CASES[i];

            if (nfilter > 0) {
                int match = 0;
                for (int j = 0; j < nfilter; j++)
                    if (strcmp(c->name, filter[j]) == 0) { match = 1; break; }
                if (!match) continue;
            }

            if (generate_pair_outputs(c, pair_outdir) != 0)
                return 1;
        }
    }

    fprintf(stderr, "\nDone. %d case(s) processed.\n", processed);
    return 0;
}
