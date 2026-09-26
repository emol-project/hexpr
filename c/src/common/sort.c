#include <stdlib.h>
#include <string.h>
#include "sort.h"

/**
 * @brief Merge two adjacent sorted runs @c base[l..m] and @c base[m+1..r] into
 *        order, using @p tmp as scratch; stable on ties.
 * @param base byte base of the element array.
 * @param l    first index of the left run.
 * @param m    last index of the left run (right run starts at @c m+1).
 * @param r    last index of the right run.
 * @param sz   element size in bytes.
 * @param cmp  qsort-style comparator (negative/zero/positive).
 * @param tmp  scratch buffer of the SAME element count as @p base (see note).
 * @note Internal helper for @ref msort. Stability comes from taking the left
 *       element when @c cmp <= 0.
 */
static void merge_run(char *base, size_t l, size_t m, size_t r,
                      size_t sz, int (*cmp)(const void *, const void *),
                      char *tmp)
{
    /* tmp is indexed by ABSOLUTE element position (tmp + k*sz, k starts at l),
     * not from 0 -- so tmp must span the whole array, and only [l..r] is
     * written/copied back here, leaving the rest of tmp untouched. */
    size_t i = l, j = m + 1, k = l;
    while (i <= m && j <= r) {
        /* cmp <= 0: left element wins ties, preserving stability */
        if (cmp(base + i * sz, base + j * sz) <= 0)
            memcpy(tmp + k++ * sz, base + i++ * sz, sz);
        else
            memcpy(tmp + k++ * sz, base + j++ * sz, sz);
    }
    while (i <= m) memcpy(tmp + k++ * sz, base + i++ * sz, sz);
    while (j <= r) memcpy(tmp + k++ * sz, base + j++ * sz, sz);
    memcpy(base + l * sz, tmp + l * sz, (r - l + 1) * sz);
}

/**
 * @brief Recursively merge-sort the inclusive range @c base[l..r].
 * @param base byte base of the element array.
 * @param l    first index of the range (inclusive).
 * @param r    last index of the range (inclusive).
 * @param sz   element size in bytes.
 * @param cmp  qsort-style comparator.
 * @param tmp  scratch buffer spanning the whole array (see @ref merge_run).
 * @note Internal helper for @ref hexpr_stable_sort. @p l and @p r are an
 *       INCLUSIVE range, so callers pass @c nmemb-1 as the top index.
 */
static void msort(char *base, size_t l, size_t r, size_t sz,
                  int (*cmp)(const void *, const void *), char *tmp)
{
    if (l >= r) return;
    /* Overflow-safe midpoint: l + (r-l)/2 instead of (l+r)/2, which could
     * overflow size_t for very large ranges. */
    size_t m = l + (r - l) / 2;
    msort(base, l, m, sz, cmp, tmp);
    msort(base, m + 1, r, sz, cmp, tmp);
    merge_run(base, l, m, r, sz, cmp, tmp);
}

/**
 * @brief Stable in-place sort with a qsort-compatible signature.
 * @param base   array to sort in place.
 * @param nmemb  number of elements.
 * @param size   element size in bytes.
 * @param compar qsort-style comparator (negative/zero/positive).
 * @note Stable (equal elements keep their input order) -- unlike @c qsort(3),
 *       whose order on ties is unspecified. Uses O(nmemb*size) scratch via a
 *       bottom-up merge sort; if that allocation fails it falls back to
 *       @c qsort, which sorts correctly but is NOT stable. O(n log n) time.
 */
void hexpr_stable_sort(void *base, size_t nmemb, size_t size,
                       int (*compar)(const void *, const void *))
{
    /* Guard is `< 2`, not `== 0`: it both short-circuits the trivially-sorted
     * 0/1-element cases AND prevents the `nmemb - 1` below from wrapping around
     * (size_t is unsigned) when nmemb == 0. */
    if (nmemb < 2) return;
    char *tmp = malloc(nmemb * size);
    if (!tmp) {
        qsort(base, nmemb, size, compar); /* unstable fallback */
        return;
    }
    msort((char *)base, 0, nmemb - 1, size, compar, tmp);
    free(tmp);
}
