#!/bin/sh
# check_fortran_coverage.sh -- three derived coverage checks over the
# Fortran seam and the exported symbol table, run from c/tests.
#
#   A. include/hexpr.h  vs  src/hexpr_fortran.h
#      every public C function must have an F77 wrapper, and every wrapper
#      must correspond to a public C function.
#
#   B. src/hexpr_fortran.h  vs  test_fortran_symbols.c
#      every declared wrapper must appear in the test's check list, and
#      nothing else may.
#
#   C. both headers  vs  nm -D on the built library
#      the library must export exactly what the two headers declare, and
#      nothing else.  Skipped when there is no library to look at.
#
# WHY THIS EXISTS
#
# test_fortran_symbols.c asks the dynamic linker whether a list of symbols
# exists in libhexpr.so.  The list is written by hand, so a wrapper added
# to hexpr_fortran.h and never added to that list is simply not checked --
# and the suite still reports PASS, which is worse than having no suite,
# because the green is read as coverage.
#
# That is not hypothetical.  When check B was written the list had drifted
# by eleven symbols: hexpr_drt_walk_ncsfs_ and hexpr_drt_walk_csfs_ from
# the CSFSET work, and the whole pair/block/subspace evaluation family
# -- nine functions -- from the evaluation API.  The list was a prefix of
# the header, correct on the day it was written and silently stale after.
#
# A hand-maintained inventory of anything drifts.  These two are derived.
#
# WHY CHECK A WAS ADDED LATER
#
# Check B compares the F77 header against the test list, so it cannot see
# the upstream gap: a public C function that never got an F77 wrapper at
# all.  csf2det needed exactly that addition by hand.  Check A was left
# out at first because it appeared to need a list of public functions that
# legitimately have no Fortran counterpart, and an unreviewed exemption
# list would just be another hand-maintained inventory.
#
# The difference was then printed and looked at, and it was empty: at the
# time of writing every one of the 45 public functions had a wrapper and
# every wrapper had a public function.  EXEMPT below is therefore empty,
# and an entry in it is a deliberate act with a reason attached rather
# than a list inherited unreviewed.
#
# WHY CHECK C WAS ADDED LATER STILL
#
# libhexpr.so used to export 142 functions.  The two headers declare 91 of
# them.  The other 51 -- internal hexpr_* entry points plus darr_*,
# dhash_*, wigner_*, calc_TensorOp_average and pair_eval_internal -- were
# visible by accident, and the documentation had been pointing readers at
# two of them.  The library is now built with -fvisibility=hidden and the
# headers mark their declarations with a visibility pragma.
#
# Check C was considered at the time and rejected: it would have needed a
# 31-entry exemption list, which is the hand-maintained inventory this
# script exists to abolish.  Once the surface was closed the expected set
# became derivable from the same two headers checks A and B already parse,
# and the exemption list became empty.  That is the only reason it can be
# written this way.
#
# It is the check that keeps the visibility work from rotting.  A new
# internal function declared outside a pragma region, or a public one that
# stops being exported because some internal header declares it first, is
# otherwise completely silent: the build succeeds and every suite passes.
#
# PARSING
#
# Both headers are stripped of comments before any name is extracted.
# Without that, a function named in prose -- "release with hexpr_drt_free()"
# -- is indistinguishable from a declaration, and an internal name
# mentioned in a public header's comment would be reported as a missing
# wrapper.  Stripping did not change the result on the day this was
# written; it removes a way for the checker to be wrong later.
#
# The name pattern requires the identifier to be followed by '(' on the
# same line.  A declaration that breaks the line between the name and its
# parameter list would be missed.  None do today; if that changes, this
# check reports it as a missing wrapper rather than passing silently.
#
# hexpr_fortran_stdout() does not end in '_' and is skipped by the wrapper
# pattern: it is a helper for the bindings, not a wrapper.

set -e

PUB=../include/hexpr.h
HDR=../src/hexpr_fortran.h
SRC=test_fortran_symbols.c

for f in "$PUB" "$HDR" "$SRC"; do
    if [ ! -f "$f" ]; then
        echo "=== Fortran seam coverage ==="
        echo "  FAIL: cannot read $f"
        echo "FAIL (1 failure(s))"
        exit 1
    fi
done

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Public C functions that are exempt from needing an F77 wrapper.
# One name per line, each with a reason. Empty by design -- see the header
# comment. Adding a name here silences a real gap, so say why.
cat > "$tmp/exempt" <<'EOF'
EOF
sed 's/#.*//' "$tmp/exempt" \
    | sed 's/[ 	]*$//' \
    | grep -v '^$' \
    | sort -u > "$tmp/exempt_clean" || :

# --- comment stripper -------------------------------------------------
# Removes /* */ (including multi-line) and //. Does not track string
# literals; neither header contains one that embeds a comment opener, and
# a wrong strip would drop a declaration and show up as a reported gap,
# not as a silent pass.
cat > "$tmp/strip.awk" <<'EOF'
BEGIN { inc = 0 }
{
    line = $0; out = ""; i = 1; n = length(line)
    while (i <= n) {
        c2 = substr(line, i, 2)
        if (inc) {
            if (c2 == "*/") { inc = 0; i += 2 } else { i++ }
        } else if (c2 == "/*") { inc = 1; i += 2 }
        else if (c2 == "//") { break }
        else { out = out substr(line, i, 1); i++ }
    }
    print out
}
EOF

# Public: any identifier hexpr_... immediately followed by '(', outside
# comments. Names ending in '_' are excluded so that this stays a list of
# C entry points even if a wrapper is ever declared in the public header.
awk -f "$tmp/strip.awk" "$PUB" \
    | grep -o 'hexpr_[a-z0-9_]* *(' \
    | sed 's/ *(.*//' \
    | grep -v '_$' \
    | sort -u > "$tmp/public_all"

comm -23 "$tmp/public_all" "$tmp/exempt_clean" > "$tmp/public"

# Declared wrappers: hexpr_..._ immediately followed by '(', outside
# comments.
awk -f "$tmp/strip.awk" "$HDR" \
    | grep -o 'hexpr_[a-z0-9_]*_ *(' \
    | sed 's/ *(.*//' \
    | sort -u > "$tmp/declared"

# The C name each wrapper corresponds to (trailing '_' removed).
sed 's/_$//' "$tmp/declared" | sort -u > "$tmp/wrapped"

# Checked: the string literal inside each check_sym() call, underscore ones.
grep -o 'check_sym(h, *"hexpr_[a-z0-9_]*_"' "$SRC" \
    | sed 's/.*"\([^"]*\)".*/\1/' \
    | sort -u > "$tmp/checked"

npub=$(wc -l < "$tmp/public"  | tr -d ' ')
nwrp=$(wc -l < "$tmp/wrapped" | tr -d ' ')
nexm=$(wc -l < "$tmp/exempt_clean" | tr -d ' ')
ndecl=$(wc -l < "$tmp/declared" | tr -d ' ')
nchk=$(wc -l < "$tmp/checked"  | tr -d ' ')

failures=0

# --- A. public header vs F77 header -----------------------------------

echo "=== F77 wrapper coverage (public header vs. F77 header) ==="

comm -23 "$tmp/public" "$tmp/wrapped" > "$tmp/unwrapped"
if [ -s "$tmp/unwrapped" ]; then
    echo "  FAIL: declared in $PUB but has no F77 wrapper:"
    sed 's/^/    /' "$tmp/unwrapped"
    echo "    -> add a wrapper in $HDR and hexpr_fortran.c,"
    echo "       or add the name to EXEMPT in this script with a reason"
    failures=$((failures + 1))
fi

comm -13 "$tmp/public" "$tmp/wrapped" > "$tmp/orphan"
if [ -s "$tmp/orphan" ]; then
    echo "  FAIL: wrapped in $HDR but not public in $PUB:"
    sed 's/^/    /' "$tmp/orphan"
    echo "    -> a wrapper over an internal function, a renamed public"
    echo "       function, or a name this script parsed wrongly"
    failures=$((failures + 1))
fi

if [ ! -s "$tmp/unwrapped" ] && [ ! -s "$tmp/orphan" ]; then
    echo "  $npub public, $nwrp wrapped, $nexm exempt, no gap        PASS"
fi

# --- B. F77 header vs test list ---------------------------------------

echo "=== F77 symbol coverage (F77 header vs. test list) ==="

comm -23 "$tmp/declared" "$tmp/checked" > "$tmp/unchecked"
if [ -s "$tmp/unchecked" ]; then
    echo "  FAIL: declared in $HDR but never checked:"
    sed 's/^/    /' "$tmp/unchecked"
    echo "    -> add a check_sym() line for each in $SRC"
    failures=$((failures + 1))
fi

comm -13 "$tmp/declared" "$tmp/checked" > "$tmp/undeclared"
if [ -s "$tmp/undeclared" ]; then
    echo "  FAIL: checked in $SRC but not declared in $HDR:"
    sed 's/^/    /' "$tmp/undeclared"
    echo "    -> a renamed or removed wrapper, or a typo in the test"
    failures=$((failures + 1))
fi

if [ ! -s "$tmp/unchecked" ] && [ ! -s "$tmp/undeclared" ]; then
    echo "  $ndecl declared, $nchk checked, no drift             PASS"
fi

# --- C. exported symbol table vs both headers -------------------------
#
# Expected exports: every name declared in the public header, plus every
# hexpr_ name declared in the F77 header.  The latter is 46, not 45: the
# wrappers plus hexpr_fortran_stdout, which fortran/src/hexpr_mod.f90
# binds to by name through the shared library and which therefore has to
# stay visible even though it is not public C API.
#
# public_all is used rather than public: EXEMPT is about which functions
# need an F77 wrapper, not about which are exported.

echo "=== Exported symbols (library vs. both headers) ==="

LIB=
for cand in ../libhexpr.so ../libhexpr.dylib; do
    if [ -f "$cand" ]; then LIB="$cand"; break; fi
done

skipc=
if [ -z "$LIB" ]; then
    skipc="no libhexpr.so or libhexpr.dylib in .. -- run make -C c first"
elif ! command -v nm > /dev/null 2>&1; then
    skipc="nm not found"
elif ! nm -D --defined-only "$LIB" > "$tmp/nm_raw" 2>/dev/null; then
    skipc="nm -D --defined-only is not supported here"
fi

if [ -n "$skipc" ]; then
    # A skip is not a failure: checks A and B read only text and must stay
    # runnable outside a build tree.
    echo "  SKIP: $skipc"
else
    awk '$2 == "T" { print $3 }' "$tmp/nm_raw" | sort -u > "$tmp/exported"

    awk -f "$tmp/strip.awk" "$HDR" \
        | grep -o 'hexpr_[a-z0-9_]* *(' \
        | sed 's/ *(.*//' \
        | sort -u > "$tmp/hdr_all"

    cat "$tmp/public_all" "$tmp/hdr_all" | sort -u > "$tmp/expected"

    nexpect=$(wc -l < "$tmp/expected" | tr -d ' ')
    nexport=$(wc -l < "$tmp/exported" | tr -d ' ')

    comm -23 "$tmp/expected" "$tmp/exported" > "$tmp/not_exported"
    comm -13 "$tmp/expected" "$tmp/exported" > "$tmp/not_declared"

    if [ -s "$tmp/not_exported" ]; then
        echo "  FAIL: declared in a header but not exported by $LIB:"
        sed 's/^/    /' "$tmp/not_exported"
        echo "    -> the declaration is outside the visibility pragma"
        echo "       region, or an internal header duplicates it and the"
        echo "       defining TU sees that copy instead (loopgen/drt.h and"
        echo "       loopgen/expr.h both do this behind #ifndef HEXPR_H_)"
        failures=$((failures + 1))
    fi

    if [ -s "$tmp/not_declared" ]; then
        echo "  FAIL: exported by $LIB but declared in no header:"
        sed 's/^/    /' "$tmp/not_declared"
        echo "    -> every exported symbol is a compatibility surface."
        echo "       Either declare it in a header, or let"
        echo "       -fvisibility=hidden hide it: check that no pragma"
        echo "       region covers it, and that any archive it comes from"
        echo "       is named in --exclude-libs BEFORE --whole-archive"
        failures=$((failures + 1))
    fi

    if [ ! -s "$tmp/not_exported" ] && [ ! -s "$tmp/not_declared" ]; then
        echo "  $nexpect expected, $nexport exported, no gap      PASS"
    fi
fi

# --- verdict ----------------------------------------------------------

if [ "$failures" -eq 0 ]; then
    echo "PASS (0 failure(s))"
    exit 0
fi

echo "  ($npub public, $ndecl declared, $nchk checked)"
echo "FAIL ($failures failure(s))"
exit 1
