#!/bin/bash
# Compiler / sanitizer matrix over c/tests `make check` (which compares with
# the frozen references in reference/). Each configuration builds a fresh
# clone of BRANCH under OUT/<name>, so the working tree is never rebuilt.
#
#   bash validation/safety/run_matrix.sh [REPO] [BRANCH] [OUT]
#
# Needs gcc and clang (MSan and TSan are clang only).
SRC=${1:-$(cd "$(dirname "$0")/../.." && pwd)}
BRANCH=${2:-release/v1.2.0}
OUT=${3:-/tmp/hxm}
mkdir -p "$OUT"
BASE="-Wall -Wextra -pedantic -g -fPIC -std=c11 -fvisibility=hidden"
run() {   # name cc cflags ldextra [env]
  local name=$1 cc=$2 cf=$3 ld=$4 envs=$5
  local d=$OUT/$name
  rm -rf "$d" && git clone -q -b "$BRANCH" "$SRC" "$d"
  {
    echo "== $name: $cc $cf"
    make -s -C "$d/c" CC=$cc CFLAGS="$cf" LDLIBS="-lm $ld" 2>&1 | grep -E "warning|error" | sort -u
    make -s -C "$d/c/tests" CC=$cc CFLAGS="$cf $ld" all 2>&1 | grep -E "warning|error" | sort -u
    ( cd "$d/c/tests" && env $envs make -s check 2>&1 | tail -40 )
  } > "$OUT/$name.log" 2>&1
  echo "$name: $(grep -E 'ALL TESTS PASSED|TEST\(S\) FAILED' "$OUT/$name.log" | tail -1)  \
sanitizer reports: $(grep -c -E 'ERROR: (Address|Memory|Thread|Leak)Sanitizer|runtime error:|WARNING: ThreadSanitizer' "$OUT/$name.log")  \
compiler warnings: $(grep -c warning "$OUT/$name.log")"
}
SAN="-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=undefined"
for cc in gcc clang; do
  for o in -O0 -O2 -O3; do run plain-$cc$o $cc "$BASE $o" ""; done
  run asan-ubsan-$cc $cc "$BASE -O1 $SAN" "$SAN" "ASAN_OPTIONS=detect_leaks=1"
done
INT="-fsanitize=signed-integer-overflow,shift,integer-divide-by-zero,float-cast-overflow,bounds -fno-sanitize-recover=all"
run int-clang clang "$BASE -O1 $INT" "$INT"
run msan-clang clang "$BASE -O1 -fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer" \
    "-fsanitize=memory"
run tsan-clang clang "$BASE -O1 -fsanitize=thread" "-fsanitize=thread"

# valgrind memcheck on the plain gcc -O2 build
if command -v valgrind > /dev/null; then
  cd "$OUT/plain-gcc-O2/c/tests" || exit 1
  export LD_LIBRARY_PATH=..
  bad=0
  for t in test_compare test_wigner test_drt test_expr test_full test_api test_csf_lookup \
           test_pair_one_4o4e test_pair_full_4o4e test_block_subspace_4o4e test_max_terms \
           test_drt_from_csfset test_csfset_pair test_ncsfs_consistency test_csf2det \
           test_pair_full_union test_pair_threads test_lk_table test_lk_twophase test_rank1; do
    valgrind -q --error-exitcode=99 --leak-check=full \
        --errors-for-leak-kinds=definite,indirect ./$t > /dev/null 2> "$OUT/vg_$t.err" \
      || { bad=$((bad + 1)); echo "valgrind: $t failed (see $OUT/vg_$t.err)"; }
  done
  echo "valgrind: $bad failing program(s)"
fi
