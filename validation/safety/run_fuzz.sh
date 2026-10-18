#!/bin/bash
# Coverage-guided fuzzing of hexpr_rank1_pair (fuzz_rank1.c) with the library
# sources instrumented, under ASan + UBSan; then the line coverage of
# csf_rank1.c reached by the corpus.
#
#   bash validation/safety/run_fuzz.sh [SECONDS] [OUT]
#
# Needs clang with libFuzzer, llvm-profdata and llvm-cov.
HERE=$(cd "$(dirname "$0")" && pwd)
C=$HERE/../../c
T=${1:-300}
OUT=${2:-/tmp/hxfuzz}
mkdir -p "$OUT" && cd "$OUT" || exit 1
INC="-I$C/include -I$C -I$C/wig -I$C/src/common -I$C/src/loopgen -I$C/src/csf"
SRCS="$C/src/common/fmt.c $C/src/common/darr.c $C/src/common/dhash.c $C/src/common/sort.c
      $C/src/loopgen/drt.c $C/src/loopgen/expr_tree.c $C/src/loopgen/expr.c
      $C/src/csf/csf_lookup.c $C/src/csf/csf_pair_core.c $C/src/csf/csf_pair.c
      $C/src/csf/csf2det.c $C/src/csf/csf_rank1.c $C/src/hexpr_api.c
      $C/wig/wigner_init.c $C/wig/wigner_3j.c $C/wig/wigner_6j.c $C/wig/wigner_9j.c"
clang -std=c11 -g -O1 -fsanitize=fuzzer,address,undefined -fno-sanitize-recover=undefined \
    $INC "$HERE/fuzz_rank1.c" $SRCS -lm -o fuzz_rank1 || exit 1
mkdir -p corpus
./fuzz_rank1 -max_total_time="$T" -max_len=96 corpus 2>&1 | grep -E "DONE|ERROR|SUMMARY|^Done"
clang -std=c11 -O0 -fprofile-instr-generate -fcoverage-mapping -fsanitize=fuzzer \
    $INC "$HERE/fuzz_rank1.c" $SRCS -lm -o fuzz_cov || exit 1
LLVM_PROFILE_FILE=cov.profraw ./fuzz_cov -runs=0 corpus > /dev/null 2>&1
llvm-profdata merge -o cov.profdata cov.profraw
llvm-cov report ./fuzz_cov -instr-profile=cov.profdata "$C/src/csf/csf_rank1.c" | tail -1
echo "lines of csf_rank1.c never reached:"
llvm-cov show ./fuzz_cov -instr-profile=cov.profdata "$C/src/csf/csf_rank1.c" \
    | grep -E '^ +[0-9]+\| +0\|'
