# =====================================================================
#  extconf.rb -- builds the Wigcoef Ruby C extension.
#
#  Usage:
#      ruby extconf.rb
#      make
#      # => produces Wigcoef.so (or .bundle) in the current directory
# =====================================================================
require 'mkmf'

# Library sources live one level up under src/ and include/.
$srcs = %w[
  wigner_ruby.c
  ../src/wigner_init.c
  ../src/wigner_3j.c
  ../src/wigner_6j.c
  ../src/wigner_9j.c
]
$INCFLAGS << " -I../include -I../src -I."

# Math library (sqrt, atan).
have_library('m', 'sqrt')

create_makefile('Wigcoef')
