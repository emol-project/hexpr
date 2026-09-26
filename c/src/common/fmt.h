/* common/fmt.h -- output formatters matching Ruby's printf / p conventions */
#ifndef HEXPR_FMT_H_
#define HEXPR_FMT_H_

#include <stdio.h>

/* Print x in Ruby's printf("%24.16e", x) format (glibc %e matches). */
void hexpr_print_double(FILE *f, double x);

/* Print a CSF as Ruby's p-style Array#inspect + newline.
 * steps[0..norb-1] are characters from {'e','u','d','f'}.
 * Example (norb=4, steps="eudf"):  ["e", "u", "d", "f"]\n  */
void hexpr_print_csf(FILE *f, const char *steps, int norb);

#endif /* HEXPR_FMT_H_ */
