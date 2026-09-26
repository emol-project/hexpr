#include "fmt.h"

/**
 * @brief Write a double to a stream in the project's fixed
 *        @c "%24.16e" format (matches Ruby's @c printf).
 * @param f destination stream (not closed; caller owns it).
 * @param x value to format.
 * @note No leading label, separator, or trailing newline is emitted -- the
 *       caller composes surrounding layout. Output width is exactly the
 *       @c %24.16e field; no NUL is involved (stream output, not a buffer).
 */
void hexpr_print_double(FILE *f, double x)
{
    fprintf(f, "%24.16e", x);
}

/**
 * @brief Write the first @p norb step characters of a CSF as a Ruby
 *        @c Array#inspect-style JSON array, followed by a newline.
 * @param f     destination stream (not closed; caller owns it).
 * @param steps step-code characters; must hold at least @p norb readable bytes.
 * @param norb  number of step characters to emit.
 * @note Reads exactly @p steps[0..norb-1]; a terminating NUL is neither required
 *       nor consulted (the count, not NUL, bounds the read). Output form is
 *       @c ["e", "u", ...]\\n. With @p norb <= 0 the result is @c []\\n.
 */
void hexpr_print_csf(FILE *f, const char *steps, int norb)
{
    fprintf(f, "[");
    for (int i = 0; i < norb; i++) {
        if (i > 0) fprintf(f, ", ");
        fprintf(f, "\"%c\"", steps[i]);
    }
    fprintf(f, "]\n");
}
