/* test_hexpr/util/compare.h -- numeric file comparator for regression tests */
#ifndef HEXPR_COMPARE_H_
#define HEXPR_COMPARE_H_

/* Compare two text files token by token.
 *
 * Tokens are split on whitespace.  For each pair of corresponding tokens:
 *   - If both parse as finite doubles via strtod -> compare numerically:
 *       |a - b| <= abs_tol  OR  |a - b| <= rel_tol * max(|a|, |b|)
 *   - Otherwise -> exact string comparison.
 *
 * Reports up to max_report mismatches to stderr with file:line info.
 * Returns total mismatch count (0 = equal within tolerance).
 * Returns -1 on I/O error. */
int compare_files(const char *path_a, const char *path_b,
                  double abs_tol, double rel_tol, int max_report);

#endif /* HEXPR_COMPARE_H_ */
