#!/bin/bash
# Remove a prefix install of hexpr by reading its install_manifest.txt.
#
# Usage: scripts/uninstall.sh PREFIX

set -uo pipefail

if [ "$#" -lt 1 ] || [ -z "${1:-}" ]; then
    echo "usage: scripts/uninstall.sh PREFIX" >&2
    exit 1
fi

PREFIX="$1"
MANIFEST="$PREFIX/share/hexpr/install_manifest.txt"

if [ ! -f "$MANIFEST" ]; then
    echo "error: no install manifest found at $MANIFEST; nothing to uninstall" >&2
    exit 1
fi

REMOVED=0
SKIPPED=0

while IFS= read -r path || [ -n "$path" ]; do
    [ -z "$path" ] && continue
    if [ -e "$path" ] || [ -L "$path" ]; then
        rm -f "$path"
        REMOVED=$((REMOVED + 1))
    else
        SKIPPED=$((SKIPPED + 1))
    fi
done < "$MANIFEST"

# Clean up empty directories under $PREFIX/share/hexpr/ only.
# Leave $PREFIX/lib, $PREFIX/include, $PREFIX/share alone (other packages
# may share them).
if [ -d "$PREFIX/share/hexpr" ]; then
    find "$PREFIX/share/hexpr" -depth -mindepth 1 -type d -empty -exec rmdir {} \; 2>/dev/null
    if [ -d "$PREFIX/share/hexpr" ] && [ -z "$(ls -A "$PREFIX/share/hexpr")" ]; then
        rmdir "$PREFIX/share/hexpr"
    fi
fi

echo "Removed $REMOVED files; skipped $SKIPPED (already gone)."
echo
echo "Uninstalled hexpr from $PREFIX."
echo "Remember to remove the HEXPR_* / LD_LIBRARY_PATH (DYLD_LIBRARY_PATH on"
echo "macOS) / PYTHONPATH / JULIA_LOAD_PATH / RUBYLIB export lines from"
echo "~/.bashrc (or wherever you put them)."
