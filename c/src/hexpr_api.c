/* hexpr_api.c -- public API: version, last_error, init/shutdown */
#define _POSIX_C_SOURCE 200112L
#include "../include/hexpr.h"
#include "wig/wigner.h"
#include "common/hexpr_error.h"
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Version                                                             */
/* ------------------------------------------------------------------ */

/**
 * @brief Return the library's semantic version as a string.
 * @return pointer to a static string literal (e.g. "1.0.0"); the caller must
 *         NOT free or modify it.
 * @note Thread-safe: returns a constant, holds no state.
 */
const char *hexpr_version_string(void) { return "1.0.0"; }

/** @brief Major version component (from @c HEXPR_VERSION_MAJOR). @return int. */
int         hexpr_version_major(void)  { return HEXPR_VERSION_MAJOR; }
/** @brief Minor version component (from @c HEXPR_VERSION_MINOR). @return int. */
int         hexpr_version_minor(void)  { return HEXPR_VERSION_MINOR; }
/** @brief Patch version component (from @c HEXPR_VERSION_PATCH). @return int. */
int         hexpr_version_patch(void)  { return HEXPR_VERSION_PATCH; }

/* ------------------------------------------------------------------ */
/* Thread-local error buffer                                           */
/* ------------------------------------------------------------------ */

#define HEXPR_ERR_BUF_SIZE 512

static _Thread_local char tl_errbuf[HEXPR_ERR_BUF_SIZE] = {0};

/**
 * @brief Return the calling thread's most recent error message.
 * @return pointer to a thread-local buffer owned by the library; the caller must
 *         NOT free it. Valid until the next error is set on the SAME thread.
 * @note Thread-local: a thread observes only errors set on its own buffer.
 *       Returns an empty string ("") when no error has been set on this thread.
 */
const char *hexpr_last_error(void)
{
    return tl_errbuf;
}

/**
 * @brief printf-style setter for the calling thread's error buffer (the message
 *        @ref hexpr_last_error subsequently returns).
 * @param fmt printf-style format string.
 * @param ... arguments for @p fmt.
 * @note NOT part of the public API: declared in @c common/hexpr_error.h and
 *       absent from the public @c include/hexpr.h; intended for library-internal
 *       callers only.
 * @note The message is truncated to @c HEXPR_ERR_BUF_SIZE-1 (511) characters by
 *       @c vsnprintf, and only the calling thread's buffer is written.
 */
void hexpr_set_error(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(tl_errbuf, HEXPR_ERR_BUF_SIZE, fmt, ap);
    va_end(ap);
}

/* ------------------------------------------------------------------ */
/* Library lifecycle                                                   */
/* ------------------------------------------------------------------ */

static int g_inited = 0;

/**
 * @brief Initialize the library (currently: build the Wigner tables via
 *        @c wigner_init()).
 * @return ::hexpr_status_t; @c HEXPR_OK in the current implementation.
 * @note Idempotent: a second call while already initialised returns @c HEXPR_OK
 *       without re-running @c wigner_init() (guarded by the static @c g_inited
 *       flag).
 * @note NOT thread-safe: @c g_inited is an unsynchronised static @c int. Call
 *       once from a single thread before any concurrent use of the library.
 */
hexpr_status_t hexpr_init(void)
{
    if (g_inited) return HEXPR_OK;
    wigner_init();
    g_inited = 1;
    return HEXPR_OK;
}

/**
 * @brief Mark the library uninitialised.
 * @note Asymmetric with @ref hexpr_init: it clears the @c g_inited guard only
 *       (sets it to 0) and invokes no Wigner teardown, so a subsequent
 *       @ref hexpr_init will call @c wigner_init() again. (Control flow only;
 *       the Wigner internals in @c wig/ are out of scope here.)
 * @note NOT thread-safe (same @c g_inited rationale as @ref hexpr_init).
 */
void hexpr_shutdown(void)
{
    g_inited = 0;
}
