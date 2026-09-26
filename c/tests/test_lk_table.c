/* test_lk_table.c -- unit test for the two-phase L(k) table POD machinery
 * (lk_pack / lk_get / lk_put / lk_clear). Self-contained: exercises lk_table.h
 * directly, no libhexpr. abort() paths are checked in forked children. */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include "lk_table.h"

static int failures = 0;
static void check(const char *label, int ok)
{
    printf("  %-58s %s\n", label, ok ? "PASS" : "FAIL");
    if (!ok) failures++;
}

/* Run fn in a forked child; return 1 iff it terminated via SIGABRT (abort()). */
static int aborts(void (*fn)(void))
{
    fflush(NULL);
    pid_t pid = fork();
    if (pid == 0) {
        FILE *n = freopen("/dev/null", "w", stderr);  /* silence the diagnostic */
        (void)n;
        fn();
        _exit(0);                                     /* did NOT abort */
    }
    int st = 0;
    waitpid(pid, &st, 0);
    return WIFSIGNALED(st) && WTERMSIG(st) == SIGABRT;
}

static void trigger_out_of_range(void) { (void)lk_pack(LK_R_GF, 0, 0, 0, 0, 0, 0); }

static void trigger_overfill(void)
{
    static lk_table_t t;
    lk_clear(&t);
    for (uint32_t i = 0; i < LK_CAP; i++) lk_put(&t, i, (double)i);  /* passes 3/4 */
}

static void trigger_integrity(void)
{
    static lk_table_t t;
    lk_clear(&t);
    uint32_t k = lk_pack(2, 2, 5, 1, 1, 1, 1);
    lk_put(&t, k, 1.0);
    lk_put(&t, k, 2.0);   /* same key, different value -> abort */
}

int main(void)
{
    printf("=== lk_table unit test ===\n");

    /* pack vs the Ruby mixed-radix formula (emol_expression.rb:268) */
    check("pack(4,4,5,1,1,1,1) == 9520201", lk_pack(4,4,5,1,1,1,1) == 9520201u);
    check("pack(1,1,1,1,1,1,1) == 2376777", lk_pack(1,1,1,1,1,1,1) == 2376777u);
    {
        uint32_t ref = (uint32_t)(((((((unsigned)2*LK_R_FF+3)*LK_R_HF+7)*LK_R_PSI+4)
                                    *LK_R_PHI+2)*LK_R_TH+3)*LK_R_THP+1);
        check("pack(2,3,7,4,2,3,1) matches mixed-radix formula",
              lk_pack(2,3,7,4,2,3,1) == ref);
    }
    check("max pack(7,7,15,15,15,7,7) == 2^24-1",
          lk_pack(7,7,15,15,15,7,7) == (uint32_t)(LK_KEY_SPACE - 1));

    /* get / put / dedup */
    static lk_table_t t;
    lk_clear(&t);
    double v;
    uint32_t k1 = lk_pack(4,4,5,1,1,1,1);
    check("get miss on empty table", lk_get(&t, k1, &v) == 0);
    lk_put(&t, k1, -1.4142135623730951);
    check("get hit after put (value exact)",
          lk_get(&t, k1, &v) == 1 && v == -1.4142135623730951);
    check("size == 1", t.size == 1);
    lk_put(&t, k1, -1.4142135623730951);   /* same value -> no-op */
    check("dedup: re-put same value keeps size 1", t.size == 1);

    /* fill many, all retrievable through collisions */
    lk_clear(&t);
    int ok = 1;
    for (uint32_t i = 0; i < 1000; i++) lk_put(&t, i * 4099u + 7u, (double)i);
    for (uint32_t i = 0; i < 1000; i++) {
        double vv;
        if (!lk_get(&t, i * 4099u + 7u, &vv) || vv != (double)i) ok = 0;
    }
    check("1000-key fill+lookup all correct", ok && t.size == 1000);

    /* abort paths, each in a forked child */
    check("out-of-range component aborts", aborts(trigger_out_of_range));
    check("overfill past 3/4 aborts", aborts(trigger_overfill));
    check("integrity: same key, different value aborts", aborts(trigger_integrity));

    printf(failures ? "%d LK_TABLE CHECK(S) FAILED\n" : "ALL LK_TABLE TESTS PASSED\n",
           failures);
    return failures ? 1 : 0;
}
