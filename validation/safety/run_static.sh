#!/bin/bash
# Static analysis of the library sources: cppcheck and gcc -fanalyzer.
#
#   bash validation/safety/run_static.sh [REPO]
#
# Run from any directory; analyses REPO/c (default: this checkout).
C=${1:-$(cd "$(dirname "$0")/../.." && pwd)}/c
cd "$C" || exit 1
INC="-Iinclude -I. -Iwig -Isrc/common -Isrc/loopgen -Isrc/csf"
FILES=$(find src wig -name '*.c' ! -name 'wigner_ruby.c' ! -name 'wigner_python.c' \
        ! -name 'run_reference.c')
echo "== cppcheck"
cppcheck --quiet --std=c11 --enable=warning,performance,portability \
    --check-level=exhaustive --inline-suppr $INC \
    --template='{file}:{line}: {severity}: {id}: {message}' $FILES 2>&1 | sort -u
echo "== gcc -fanalyzer"
for f in $FILES; do
  gcc -std=c11 -O1 -fanalyzer -c -o /dev/null -fPIC $INC "$f" 2>&1
done | grep -E "warning:|error:" | sort | uniq -c | sort -rn
