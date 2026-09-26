#!/bin/bash
# Copy libhexpr.so and language bindings into a prefix directory.
#
# Usage: scripts/install.sh [PREFIX]
#   PREFIX defaults to $HOME/local/hexpr.

set -euo pipefail

cd "$(dirname "$0")/.."

PREFIX="${1:-$HOME/local/hexpr}"

# Platform shared-library name and runtime loader env var.
case "$(uname -s)" in
  Darwin) SOEXT=dylib; LIBVAR=DYLD_LIBRARY_PATH ;;
  *)      SOEXT=so;    LIBVAR=LD_LIBRARY_PATH ;;
esac

if [ ! -f "c/libhexpr.$SOEXT" ]; then
    echo "error: c/libhexpr.$SOEXT not found; build it first with 'make -C c'" >&2
    exit 1
fi

mkdir -p "$PREFIX/lib"
mkdir -p "$PREFIX/include"
mkdir -p "$PREFIX/share/hexpr/python"
mkdir -p "$PREFIX/share/hexpr/julia/HExpr/src"
mkdir -p "$PREFIX/share/hexpr/ruby/lib"
mkdir -p "$PREFIX/share/hexpr/fortran"
mkdir -p "$PREFIX/share/hexpr/example_ci"

# Track every file installed for the manifest.
INSTALLED=()

install_file() {
    # install_file <src> <dst>
    cp "$1" "$2"
    INSTALLED+=("$2")
}

# --- Library + headers ---
install_file "c/libhexpr.$SOEXT" "$PREFIX/lib/libhexpr.$SOEXT"

install_file c/include/hexpr.h "$PREFIX/include/hexpr.h"

# --- Python ---
PY_DEST="$PREFIX/share/hexpr/python/hexpr"
rm -rf "$PY_DEST"
mkdir -p "$PY_DEST"
# Copy *.py files only; skip __pycache__, *.pyc, and *.egg-info.
while IFS= read -r -d '' f; do
    rel="${f#python/hexpr/}"
    dest="$PY_DEST/$rel"
    mkdir -p "$(dirname "$dest")"
    install_file "$f" "$dest"
done < <(
    find python/hexpr -type f \
        ! -path '*/__pycache__/*' \
        ! -name '*.pyc' \
        ! -path '*.egg-info/*' \
        -print0
)

# --- Julia (package layout: HExpr/src/HExpr.jl + HExpr/Project.toml) ---
install_file julia/src/HExpr.jl "$PREFIX/share/hexpr/julia/HExpr/src/HExpr.jl"
install_file julia/Project.toml "$PREFIX/share/hexpr/julia/HExpr/Project.toml"

# --- Ruby ---
install_file ruby/lib/hexpr.rb "$PREFIX/share/hexpr/ruby/lib/hexpr.rb"

# --- Fortran (source module only) ---
install_file fortran/src/hexpr_mod.f90 "$PREFIX/share/hexpr/fortran/hexpr_mod.f90"

# --- example_ci (worked end-to-end examples; same exclusions as Python) ---
EX_DEST="$PREFIX/share/hexpr/example_ci"
rm -rf "$EX_DEST"
mkdir -p "$EX_DEST"
while IFS= read -r -d '' f; do
    rel="${f#example_ci/}"
    dest="$EX_DEST/$rel"
    mkdir -p "$(dirname "$dest")"
    install_file "$f" "$dest"
done < <(
    find example_ci -type f \
        ! -path '*/__pycache__/*' \
        ! -name '*.pyc' \
        ! -path '*.egg-info/*' \
        -print0
)

# --- Write install manifest (lists every installed file, including itself) ---
MANIFEST="$PREFIX/share/hexpr/install_manifest.txt"
{
    for p in "${INSTALLED[@]}"; do
        printf '%s\n' "$p"
    done
    printf '%s\n' "$MANIFEST"
} > "$MANIFEST"

# --- Final env-var hint ---
cat <<EOF
Installed hexpr to $PREFIX

To use hexpr, add the following to your shell rc (e.g. ~/.bashrc):

  export HEXPR_PREFIX=$PREFIX
  export HEXPR_LIBRARY=\$HEXPR_PREFIX/lib/libhexpr.$SOEXT
  export $LIBVAR=\$HEXPR_PREFIX/lib:\$$LIBVAR
  export PYTHONPATH=\$HEXPR_PREFIX/share/hexpr/python:\$PYTHONPATH
  export JULIA_LOAD_PATH=\$HEXPR_PREFIX/share/hexpr/julia:\$JULIA_LOAD_PATH
  export RUBYLIB=\$HEXPR_PREFIX/share/hexpr/ruby/lib:\$RUBYLIB

For Fortran, the source module is at:
  \$HEXPR_PREFIX/share/hexpr/fortran/hexpr_mod.f90
Compile it together with your Fortran project as needed.

To uninstall later:
  bash scripts/uninstall.sh $PREFIX
EOF
