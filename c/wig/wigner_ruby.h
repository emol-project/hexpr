/* ====================================================================
 *  wigner_ruby.h -- Ruby C-extension entry points.
 *
 *  These are the Init_* and wrapper functions exposed to MRI's
 *  module loader.  The pure C interface lives in wigner.h; this
 *  header is only for code participating in the Ruby build.
 * ====================================================================
 */
#ifndef WIGNER_RUBY_H_
#define WIGNER_RUBY_H_

#include "ruby.h"

VALUE wigcoef_rb_init  (VALUE self);
VALUE wigcoef_rb_delta (VALUE self, VALUE ia, VALUE ib, VALUE ic);
VALUE wigcoef_rb_f3l   (VALUE self, VALUE ia, VALUE ib, VALUE ic);
VALUE wigcoef_rb_wig6j (VALUE self,
                        VALUE i1, VALUE i2, VALUE i3,
                        VALUE i4, VALUE i5, VALUE i6);
VALUE wigcoef_rb_wig9j (VALUE self,
                        VALUE i1, VALUE i2, VALUE i3,
                        VALUE i4, VALUE i5, VALUE i6,
                        VALUE i7, VALUE i8, VALUE i9);

void  Init_Wigcoef(void);

#endif /* WIGNER_RUBY_H_ */
