#!/usr/bin/env bash
# Verify INDOTERM-Julia input data provenance, and unpack the derived CSV if needed.
#
# Run from the repository root:  bash scripts/verify_data.sh
set -euo pipefail

cd "$(dirname "$0")/.."
MANIFEST="data/checksums.sha256"
[ -f "$MANIFEST" ] || { echo "FAIL: $MANIFEST not found"; exit 1; }

echo "== committed inputs =="
( cd data && sha256sum -c checksums.sha256 )

# The pipeline reads data/national_data.csv, which is gitignored and regenerated
# from data/national_data.zip. Unpack it if absent, then check it against the
# 'derived' line recorded in the manifest.
if [ ! -f data/national_data.csv ]; then
    echo "== data/national_data.csv absent; unpacking from national_data.zip =="
    unzip -oq data/national_data.zip -d data/
fi

EXPECTED=$(awk '/^# derived[[:space:]]+national_data\.csv/ {print $NF}' "$MANIFEST")
if [ -z "$EXPECTED" ]; then
    echo "FAIL: no 'derived national_data.csv' line in $MANIFEST"; exit 1
fi
ACTUAL=$(sha256sum data/national_data.csv | cut -d' ' -f1)

echo "== derived artifact =="
if [ "$EXPECTED" = "$ACTUAL" ]; then
    echo "national_data.csv: OK"
else
    echo "national_data.csv: FAILED"
    echo "  expected $EXPECTED"
    echo "  actual   $ACTUAL"
    exit 1
fi
echo "All data provenance checks passed."
