#include <stdlib.h>
#include "darr.h"

/**
 * @brief Initialize a ::darr_t to the empty state (no allocation).
 * @param a array to initialize (caller-owned storage, e.g. a stack variable).
 * @note Must be called before any other @c darr_* call on @p a. Pairs with
 *       @ref darr_free. Does not allocate; the first @ref darr_push does.
 */
void darr_init(darr_t *a)
{
    a->data = NULL;
    a->size = 0;
    a->cap  = 0;
}

/**
 * @brief Grow the backing store: capacity 0 -> 4, otherwise doubled.
 * @param a array whose capacity is increased.
 * @warning Calls @c realloc and may move the buffer, invalidating any pointer
 *          previously obtained from @ref darr_get. Internal helper for
 *          @ref darr_push.
 */
static void darr_grow(darr_t *a)
{
    int new_cap = (a->cap == 0) ? 4 : a->cap * 2;
    a->data = realloc(a->data, (size_t)new_cap * sizeof(void *));
    a->cap  = new_cap;
}

/**
 * @brief Append one pointer to the array, growing if full.
 * @param a    target array (must have been @ref darr_init "initialized").
 * @param item pointer to store; stored BY VALUE (the pointer is copied, the
 *             pointee is not). The array does not take ownership of @p item.
 * @note Amortized O(1); a growth step is O(size) and, by reallocating, may
 *       invalidate pointers previously returned by @ref darr_get (indices stay
 *       valid). @p item may be NULL -- this container does not use NULL as a
 *       sentinel.
 */
void darr_push(darr_t *a, void *item)
{
    if (a->size == a->cap) darr_grow(a);
    a->data[a->size++] = item;
}

/**
 * @brief Return the pointer stored at index @p i.
 * @param a array to read.
 * @param i index; the caller must ensure @c 0 <= i < @ref darr_size.
 * @return the stored pointer (as passed to @ref darr_push).
 * @warning No bounds checking is performed; an out-of-range @p i is undefined
 *          behaviour.
 */
void *darr_get(const darr_t *a, int i)
{
    return a->data[i];
}

/**
 * @brief Number of elements currently stored.
 * @param a array to query. @return element count (@c a->size).
 */
int darr_size(const darr_t *a)
{
    return a->size;
}

/**
 * @brief Logically empty the array, keeping the allocation for reuse.
 * @param a array to clear.
 * @note Sets size to 0 but retains capacity (no @c free); stored pointers are
 *       dropped without being freed -- free any owned pointees first. Contrast
 *       @ref darr_free, which releases the backing store.
 */
void darr_clear(darr_t *a)
{
    a->size = 0;
}

/**
 * @brief Release the backing store and reset to the empty state.
 * @param a array to free.
 * @note Frees only the pointer array, NOT the pointees -- if @p a owned its
 *       elements, free them before calling. Safe to @ref darr_push again
 *       afterwards (a is left as if freshly @ref darr_init "initialized").
 */
void darr_free(darr_t *a)
{
    free(a->data);
    a->data = NULL;
    a->size = 0;
    a->cap  = 0;
}
