/* Self-test for compare_files: three cases. */
#include <stdio.h>
#include <stdlib.h>
#include "util/compare.h"

#define TMP_A "_cmp_a.txt"
#define TMP_B "_cmp_b.txt"

static void write_file(const char *path, const char *content)
{
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); exit(1); }
    fputs(content, f);
    fclose(f);
}

static int run(const char *label, const char *a, const char *b,
               double atol, double rtol, int expect_mismatch)
{
    write_file(TMP_A, a);
    write_file(TMP_B, b);
    int got = compare_files(TMP_A, TMP_B, atol, rtol, 5);
    int ok = expect_mismatch ? (got > 0) : (got == 0);
    printf("  %-45s %s\n", label, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

int main(void)
{
    int failures = 0;
    printf("=== compare_files self-test ===\n");

    /* (a) Byte-equal files -> expect 0 mismatches */
    failures += run(
        "(a) byte-equal",
        " CSF List \n[\"e\", \"u\", \"d\", \"f\"]\n"
        "       1             0     1.2500000000000000e+00 \n",
        " CSF List \n[\"e\", \"u\", \"d\", \"f\"]\n"
        "       1             0     1.2500000000000000e+00 \n",
        1e-10, 1e-10, 0);

    /* (b) Numeric noise within tolerance -> expect 0 mismatches */
    failures += run(
        "(b) numeric noise within tol=1e-6",
        "  1.0000000000000000e+00\n",
        "  1.0000003000000000e+00\n",
        1e-6, 1e-6, 0);

    /* (c) Real mismatch outside tolerance -> expect >0 mismatches */
    failures += run(
        "(c) real mismatch outside tol=1e-10",
        "  1.0000000000000000e+00\n",
        "  2.0000000000000000e+00\n",
        1e-10, 1e-10, 1);

    /* (d) String token mismatch -> expect >0 mismatches */
    failures += run(
        "(d) string mismatch",
        "SOCI\n",
        "FOCI\n",
        1e-10, 1e-10, 1);

    /* (e) '-' (DRT nil marker) vs '-' -> expect 0 mismatches */
    failures += run(
        "(e) '-' nil markers match",
        "      -      3\n",
        "      -      3\n",
        1e-10, 1e-10, 0);

    /* (f) '-' vs '0' -> mismatch (different semantics in DRT) */
    failures += run(
        "(f) '-' vs '0' is a mismatch",
        "      -\n",
        "      0\n",
        1e-10, 1e-10, 1);

    remove(TMP_A);
    remove(TMP_B);

    printf("%s (%d failure(s))\n",
           failures == 0 ? "PASS" : "FAIL", failures);
    return failures ? 1 : 0;
}
