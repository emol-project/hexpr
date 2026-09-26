/* ====================================================================
 *  wigner_ruby.c -- thin Ruby C-extension layer over the pure-C
 *                   Wigner library.
 *
 *  The library no longer carries Fortran-style pointer arguments,
 *  so the wrapper just unpacks INT values and forwards them.
 *  Module name and arity match the original (Wigcoef.init, .delta,
 *  .f3l, .wig6j, .wig9x).  The 9j entry is also exposed under the
 *  more conventional name "wig9j"; "wig9x" is kept as an alias for
 *  backwards compatibility.
 * ====================================================================
 */
#include "ruby.h"
#include "wigner.h"
#include "wigner_ruby.h"

VALUE wigcoef_rb_init(VALUE self)
{
    (void)self;
    wigner_init();
    return Qnil;
}

VALUE wigcoef_rb_delta(VALUE self, VALUE ia, VALUE ib, VALUE ic)
{
    (void)self;
    double d = wigner_delta(NUM2INT(ia), NUM2INT(ib), NUM2INT(ic));
    return rb_float_new(d);
}

VALUE wigcoef_rb_f3l(VALUE self, VALUE ia, VALUE ib, VALUE ic)
{
    (void)self;
    double d = wigner_f3l(NUM2INT(ia), NUM2INT(ib), NUM2INT(ic));
    return rb_float_new(d);
}

VALUE wigcoef_rb_wig6j(VALUE self,
                       VALUE i1, VALUE i2, VALUE i3,
                       VALUE i4, VALUE i5, VALUE i6)
{
    (void)self;
    double d = wigner_6j(NUM2INT(i1), NUM2INT(i2), NUM2INT(i3),
                         NUM2INT(i4), NUM2INT(i5), NUM2INT(i6));
    return rb_float_new(d);
}

VALUE wigcoef_rb_wig9j(VALUE self,
                       VALUE i1, VALUE i2, VALUE i3,
                       VALUE i4, VALUE i5, VALUE i6,
                       VALUE i7, VALUE i8, VALUE i9)
{
    (void)self;
    double d = wigner_9j(NUM2INT(i1), NUM2INT(i2), NUM2INT(i3),
                         NUM2INT(i4), NUM2INT(i5), NUM2INT(i6),
                         NUM2INT(i7), NUM2INT(i8), NUM2INT(i9));
    return rb_float_new(d);
}

void Init_Wigcoef(void)
{
    VALUE module = rb_define_module("Wigcoef");

    rb_define_module_function(module, "init",  wigcoef_rb_init,  0);
    rb_define_module_function(module, "delta", wigcoef_rb_delta, 3);
    rb_define_module_function(module, "f3l",   wigcoef_rb_f3l,   3);
    rb_define_module_function(module, "wig6j", wigcoef_rb_wig6j, 6);
    rb_define_module_function(module, "wig9j", wigcoef_rb_wig9j, 9);
    /* Backwards-compatible alias. */
    rb_define_module_function(module, "wig9x", wigcoef_rb_wig9j, 9);
}
