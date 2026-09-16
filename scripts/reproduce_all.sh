#!/usr/bin/env bash
# Reproduce every published INDOTERM-Julia result, from a clean checkout.
#
#   bash scripts/reproduce_all.sh
#
# Stages, in order — each must pass before the next runs:
#   1. input data provenance (sha256)
#   2. environment instantiation from the committed Manifest.toml
#   3. fast gates        (seconds)
#   4. integration gates (minutes — benchmark solve must replicate)
#   5. full V&V gates    (hours  — V1/V2/V3/V5/V9 + homogeneity + arclength)
#   6. known-open gates  (reported, never fatal: V4, V6)
#
# Expect several hours end to end. Logs land in logs/gates/.
set -euo pipefail

cd "$(dirname "$0")/.."
JULIA="${JULIA:-julia}"

echo "════ 1/6  input data provenance ════"
bash scripts/verify_data.sh

echo
echo "════ 2/6  environment ════"
if [ ! -f Manifest.toml ]; then
    echo "WARNING: Manifest.toml is not committed — this run is NOT reproducible."
    echo "         Generate and commit it, then re-run:"
    echo "           $JULIA --project=. -e 'using Pkg; Pkg.instantiate()'"
    echo "           git add -f Manifest.toml && git commit -m 'Pin dependency versions'"
fi
"$JULIA" --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.status()'

for tier in fast integration full nongating; do
    case $tier in
        fast)       n="3/6" ;;
        integration) n="4/6" ;;
        full)       n="5/6" ;;
        nongating)  n="6/6" ;;
    esac
    echo
    echo "════ $n  gates: $tier ════"
    bash scripts/run_gates.sh "$tier"
done

echo
echo "All reproduction stages complete. Gate logs: logs/gates/"
