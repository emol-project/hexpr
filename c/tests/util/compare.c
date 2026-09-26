#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compare.h"

#define LINE_MAX_LEN 8192

/* Advance *p past whitespace; then consume and NUL-terminate the next
 * non-whitespace run.  Returns pointer to the token, or NULL if none. */
static char *next_tok(char **p)
{
    while (**p == ' ' || **p == '\t' || **p == '\r' || **p == '\n')
        (*p)++;
    if (**p == '\0') return NULL;
    char *start = *p;
    while (**p && **p != ' ' && **p != '\t' && **p != '\r' && **p != '\n')
        (*p)++;
    if (**p) *(*p)++ = '\0';
    return start;
}

/* Returns 1 and sets *val if tok is a valid finite double, 0 otherwise. */
static int parse_double(const char *tok, double *val)
{
    char *end;
    if (!tok || tok[0] == '\0') return 0;
    *val = strtod(tok, &end);
    return (end != tok && *end == '\0' && isfinite(*val));
}

static int nums_match(double a, double b, double abs_tol, double rel_tol)
{
    double diff = fabs(a - b);
    if (diff <= abs_tol) return 1;
    double ref = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
    return (ref > 0.0 && diff <= rel_tol * ref);
}

int compare_files(const char *path_a, const char *path_b,
                  double abs_tol, double rel_tol, int max_report)
{
    FILE *fa = fopen(path_a, "r");
    FILE *fb = fopen(path_b, "r");
    if (!fa || !fb) {
        if (fa) fclose(fa);
        if (fb) fclose(fb);
        fprintf(stderr, "compare_files: cannot open '%s' or '%s'\n",
                path_a, path_b);
        return -1;
    }

    char la[LINE_MAX_LEN], lb[LINE_MAX_LEN];
    int  lineno = 0, mismatches = 0;

    while (1) {
        char *ra = fgets(la, sizeof(la), fa);
        char *rb = fgets(lb, sizeof(lb), fb);
        if (!ra && !rb) break;
        if (!ra || !rb) {
            if (mismatches < max_report)
                fprintf(stderr, "%s:%d: file length mismatch\n",
                        !ra ? path_a : path_b, lineno + 1);
            mismatches++;
            break;
        }
        lineno++;

        char *pa = la, *pb = lb;
        while (1) {
            char *ta = next_tok(&pa);
            char *tb = next_tok(&pb);
            if (!ta && !tb) break;
            if (!ta || !tb) {
                if (mismatches < max_report)
                    fprintf(stderr, "%s:%d: token count mismatch\n",
                            path_a, lineno);
                mismatches++;
                break;
            }
            double va, vb;
            if (parse_double(ta, &va) && parse_double(tb, &vb)) {
                if (!nums_match(va, vb, abs_tol, rel_tol)) {
                    if (mismatches < max_report)
                        fprintf(stderr,
                                "%s:%d: numeric mismatch: '%s' vs '%s'\n",
                                path_a, lineno, ta, tb);
                    mismatches++;
                }
            } else {
                if (strcmp(ta, tb) != 0) {
                    if (mismatches < max_report)
                        fprintf(stderr,
                                "%s:%d: string mismatch: '%s' vs '%s'\n",
                                path_a, lineno, ta, tb);
                    mismatches++;
                }
            }
        }
    }

    fclose(fa);
    fclose(fb);
    return mismatches;
}
