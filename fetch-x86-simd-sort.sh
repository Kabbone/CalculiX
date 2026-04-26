#!/bin/sh
# Clone x86-simd-sort next to this repository (where src/Makefile expects it).
# Run from the CalculiX repo root.
set -e

DEST="$(dirname "$0")/../x86-simd-sort"

if [ -d "$DEST" ]; then
    echo "x86-simd-sort already present at $DEST"
    exit 0
fi

git clone https://github.com/intel/x86-simd-sort "$DEST"
