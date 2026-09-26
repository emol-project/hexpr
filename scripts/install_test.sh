#!/bin/bash
# Verify a prefix install works end-to-end.
#
# Builds nothing: expects c/libhexpr.<so|dylib> to exist already.
# Creates install_test_<timestamp>/ in the repo root, installs there,
# then runs a smoke test for each of C, Python, Julia, Ruby.

set -uo pipefail

cd "$(dirname "$0")/.."

# Platform shared-library name and runtime loader env var.
case "$(uname -s)" in
  Darwin) SOEXT=dylib; LIBVAR=DYLD_LIBRARY_PATH ;;
  *)      SOEXT=so;    LIBVAR=LD_LIBRARY_PATH ;;
esac

TS="$(date +%s)"
PREFIX="$(pwd)/install_test_${TS}"

echo "=== Prefix install verification ==="
echo "PREFIX=$PREFIX"
echo

echo "--- Build check ---"
if [ ! -f "c/libhexpr.$SOEXT" ]; then
    echo "[FAIL] c/libhexpr.$SOEXT does not exist; run 'make -C c' first"
    exit 1
fi
echo "[OK] c/libhexpr.$SOEXT exists"
echo

echo "--- Install ---"
mkdir -p "$PREFIX"
INSTALL_LOG="$PREFIX/install_err.log"
if scripts/install.sh "$PREFIX" >/dev/null 2>"$INSTALL_LOG"; then
    echo "[OK] scripts/install.sh succeeded"
else
    echo "[FAIL] scripts/install.sh failed:"
    sed 's/^/    /' "$INSTALL_LOG"
    exit 1
fi
echo

export HEXPR_PREFIX="$PREFIX"
export HEXPR_LIBRARY="$HEXPR_PREFIX/lib/libhexpr.$SOEXT"
export $LIBVAR="$HEXPR_PREFIX/lib:${!LIBVAR:-}"
export PYTHONPATH="$HEXPR_PREFIX/share/hexpr/python:${PYTHONPATH:-}"
export JULIA_LOAD_PATH="$HEXPR_PREFIX/share/hexpr/julia:${JULIA_LOAD_PATH:-}"
export RUBYLIB="$HEXPR_PREFIX/share/hexpr/ruby/lib:${RUBYLIB:-}"

echo "--- Environment (script-local) ---"
echo "HEXPR_PREFIX=$HEXPR_PREFIX"
echo "HEXPR_LIBRARY=$HEXPR_LIBRARY"
echo "$LIBVAR=${!LIBVAR}"
echo "PYTHONPATH=$PYTHONPATH"
echo "JULIA_LOAD_PATH=$JULIA_LOAD_PATH"
echo "RUBYLIB=$RUBYLIB"
echo

echo "--- Smoke tests ---"

PASS=0
FAIL=0
SKIP=0

have() { command -v "$1" >/dev/null 2>&1; }

# --- C ---
C_ERR_FILE="$PREFIX/c_err.log"
C_OUT_FILE="$PREFIX/c_out.log"
cat > "$PREFIX/smoke.c" <<'EOF'
#include <stdio.h>
#include "hexpr.h"
int main(void) {
    if (hexpr_init() != HEXPR_OK) {
        fprintf(stderr, "hexpr_init failed\n");
        return 1;
    }
    printf("C: %s\n", hexpr_version_string());
    hexpr_shutdown();
    return 0;
}
EOF
(
    # Bake an rpath so the smoke binary is self-sufficient on macOS (SIP
    # strips DYLD_LIBRARY_PATH from many contexts); harmless on Linux.
    gcc -I"$PREFIX/include" "$PREFIX/smoke.c" -L"$PREFIX/lib" -lhexpr \
        -Wl,-rpath,"$PREFIX/lib" \
        -o "$PREFIX/smoke" 2>"$C_ERR_FILE" || exit 2
    env "$LIBVAR=$PREFIX/lib:${!LIBVAR:-}" \
        "$PREFIX/smoke" >"$C_OUT_FILE" 2>"$C_ERR_FILE"
)
C_RC=$?
if [ $C_RC -eq 0 ]; then
    C_OUT="$(cat "$C_OUT_FILE")"
    PASS=$((PASS + 1))
    printf "[C]      hexpr %s\n" "${C_OUT#C: }"
else
    FAIL=$((FAIL + 1))
    if [ -s "$C_ERR_FILE" ]; then
        C_ERR="$(head -n 3 "$C_ERR_FILE" | tr '\n' ' ')"
    else
        C_ERR="exit $C_RC"
    fi
    printf "[C]      FAIL: %s\n" "$C_ERR"
fi

# --- Python ---
if ! have python3; then
    SKIP=$((SKIP + 1))
    printf "[Python] SKIP: python3 not found\n"
else
PY_OUT=""
PY_ERR_FILE="$PREFIX/py_err.log"
PY_OUT="$(python3 -c 'import hexpr; print("hexpr", hexpr.version())' 2>"$PY_ERR_FILE")"
PY_RC=$?
if [ $PY_RC -eq 0 ]; then
    PASS=$((PASS + 1))
    printf "[Python] %s\n" "$PY_OUT"
else
    # Fall back to import + init only
    PY_OUT="$(python3 -c '
import hexpr
try:
    hexpr.init()
    print("OK")
    hexpr.shutdown()
except AttributeError:
    print("OK")
' 2>"$PY_ERR_FILE")"
    PY_RC=$?
    if [ $PY_RC -eq 0 ]; then
        PASS=$((PASS + 1))
        printf "[Python] %s\n" "$PY_OUT"
    else
        FAIL=$((FAIL + 1))
        if [ -s "$PY_ERR_FILE" ]; then
            PY_ERR="$(head -n 3 "$PY_ERR_FILE" | tr '\n' ' ')"
        else
            PY_ERR="exit $PY_RC"
        fi
        printf "[Python] FAIL: %s\n" "$PY_ERR"
    fi
fi
fi

# --- Julia ---
if ! have julia; then
    SKIP=$((SKIP + 1))
    printf "[Julia]  SKIP: julia not found\n"
else
JL_ERR_FILE="$PREFIX/jl_err.log"
JL_OUT="$(julia --startup-file=no -e \
    'using HExpr; HExpr.init(); println("hexpr ", HExpr.version_string()); HExpr.shutdown()' \
    2>"$JL_ERR_FILE")"
JL_RC=$?
if [ $JL_RC -eq 0 ]; then
    PASS=$((PASS + 1))
    printf "[Julia]  %s\n" "$JL_OUT"
else
    FAIL=$((FAIL + 1))
    if [ -s "$JL_ERR_FILE" ]; then
        JL_ERR="$(head -n 3 "$JL_ERR_FILE" | tr '\n' ' ')"
    else
        JL_ERR="exit $JL_RC"
    fi
    printf "[Julia]  FAIL: %s\n" "$JL_ERR"
fi
fi

# --- Ruby ---
RB_ERR_FILE="$PREFIX/rb_err.log"
if ! have ruby; then
    SKIP=$((SKIP + 1))
    printf "[Ruby]   SKIP: ruby not found\n"
elif ! ruby -e 'require "ffi"' >/dev/null 2>&1; then
    SKIP=$((SKIP + 1))
    printf "[Ruby]   SKIP: ffi gem not installed; run: gem install ffi\n"
else
    RB_OUT="$(ruby -e \
        'require "hexpr"; HExpr.init; puts "hexpr " + HExpr.version_string; HExpr.shutdown' \
        2>"$RB_ERR_FILE")"
    RB_RC=$?
    if [ $RB_RC -eq 0 ]; then
        PASS=$((PASS + 1))
        printf "[Ruby]   %s\n" "$RB_OUT"
    else
        FAIL=$((FAIL + 1))
        if [ -s "$RB_ERR_FILE" ]; then
            RB_ERR="$(head -n 3 "$RB_ERR_FILE" | tr '\n' ' ')"
        else
            RB_ERR="exit $RB_RC"
        fi
        printf "[Ruby]   FAIL: %s\n" "$RB_ERR"
    fi
fi

echo
echo "--- Summary ---"
TOTAL=4
if [ $FAIL -eq 0 ] && [ $SKIP -eq 0 ]; then
    echo "PASS: $PASS/$TOTAL"
    RC=0
elif [ $FAIL -eq 0 ]; then
    echo "PASS: $PASS/$TOTAL ($SKIP skipped)"
    RC=0
else
    echo "FAIL: $PASS/$TOTAL ($SKIP skipped, $FAIL failed)"
    RC=1
fi
echo
echo "The $(basename "$PREFIX")/ directory is scratch; delete it when done."

exit $RC
