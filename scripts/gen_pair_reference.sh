#!/bin/bash
# Generate per-pair reference data for hexpr cases.
#
# Usage:
#   gen_pair_reference.sh <case1> [<case2> ...]
#
# Output goes to $PAIR_OUTDIR/<case>/{pair,pair_one,pair_summary}.txt
# Log goes to $LOG (also tail-printed at the end).
#
# Designed for use with `at` for scheduled runs:
#   echo "$HOME/src/claude/hexpr/scripts/gen_pair_reference.sh s_foci_small t_foci_small" | at 01:00

set -euo pipefail

# Platform runtime loader env var (Linux LD_LIBRARY_PATH / macOS DYLD_*).
case "$(uname -s)" in
  Darwin) LIBVAR=DYLD_LIBRARY_PATH ;;
  *)      LIBVAR=LD_LIBRARY_PATH ;;
esac

# === Settings ============================================================
REPO_ROOT="$HOME/src/claude/hexpr"
DIR="$REPO_ROOT/c"                                  # working dir for run
LIB_PATH="$REPO_ROOT/c"                             # dir holding libhexpr
RUN_REF="$REPO_ROOT/c/run_reference"                # binary
PAIR_OUTDIR="$REPO_ROOT/reference_pair"             # output target
LOG="/tmp/pair_gen_$(date +%Y%m%d_%H%M%S).log"      # log file
COM=("$RUN_REF" --pair --pair-outdir "$PAIR_OUTDIR")  # command + base args
# =========================================================================

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <case_name> [<case_name> ...]" >&2
    exit 2
fi

if [ ! -x "$RUN_REF" ]; then
    echo "ERROR: $RUN_REF not found or not executable" >&2
    echo "Build with: make -C $DIR" >&2
    exit 1
fi

echo "log: $LOG"

# === Run =================================================================
{
    echo "=== gen_pair_reference.sh ==="
    echo "started:   $(date)"
    echo "repo_root: $REPO_ROOT"
    echo "cases:     $*"
    echo ""

    cd "$DIR"
    env "$LIBVAR=$LIB_PATH" "${COM[@]}" "$@"

    echo ""
    echo "finished:  $(date)"
} >> "$LOG" 2>&1
# =========================================================================

echo "=== Last lines of log ==="
tail -10 "$LOG"
