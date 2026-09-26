/* loopgen/expr.h -- Energy expression internal header */
#ifndef HEXPR_EXPR_H_
#define HEXPR_EXPR_H_

#include <stdio.h>
#include <stdint.h>
#include "drt.h"    /* which includes hexpr.h */

/* hexpr_expr_t and the public expression API come from hexpr.h, reached
 * through drt.h. They used to be duplicated here behind #ifndef HEXPR_H_;
 * see the comment in drt.h for what that copy cost and why it is gone. */

/* ------------------------------------------------------------------ */
/* Internal only: not declared in hexpr.h, hidden in libhexpr.so       */
/* ------------------------------------------------------------------ */

/* hexpr_expr_free is internal (hexpr_expr_destroy is the public alias) */
void hexpr_expr_free(hexpr_expr_t *expr);

#endif /* HEXPR_EXPR_H_ */
