#!/usr/bin/env bash
# Run INDOTERM-Julia gates for one tier, per test/gates.tsv.
#
#   bash scripts/run_gates.sh fast          # every tier: fast | integration | full | nongating
#   JULIA=julia-1.10 bash scripts/run_gates.sh full
#
# Exit status is 0 only if every gate in the tier passed. `nongating` always
# reports but never fails the run (those gates are known-open; see VV_PLAN.md).
#
# WHY THE DUAL CRITERION: the V-gate scripts print their verdict and exit 0 even
# when they FAIL. Keying CI on the exit code alone would produce a green build
# that certifies nothing. Gates marked `marker` in the manifest must therefore
# also emit a ✅ and must not emit ❌/FAILED. Converting these scripts to real
# @testset gates is tracked in README "Known gap: gate verdict contracts".
set -uo pipefail

cd "$(dirname "$0")/.."
TIER="${1:-}"
JULIA="${JULIA:-julia}"
MANIFEST="test/gates.tsv"

case "$TIER" in
  fast|integration|full|nongating) ;;
  *) echo "usage: $0 {fast|integration|full|nongating}" >&2; exit 2 ;;
esac
[ -f "$MANIFEST" ] || { echo "FAIL: $MANIFEST not found" >&2; exit 2; }

LOGDIR="${GATE_LOGDIR:-logs/gates}"
mkdir -p "$LOGDIR"

pass=0; fail=0; failed_names=""
while IFS=$'\t' read -r tier script verdict gate desc; do
    case "$tier" in ''|\#*) continue ;; esac
    [ "$tier" = "$TIER" ] || continue

    name="$(basename "$script" .jl)"
    log="$LOGDIR/${name}.log"
    printf '\n=== [%s] %s (gate %s) — %s\n' "$tier" "$name" "$gate" "$desc"

    if [ ! -f "$script" ]; then
        echo "  ✗ MISSING: $script"; fail=$((fail+1)); failed_names="$failed_names $name"; continue
    fi

    "$JULIA" --project=. "$script" >"$log" 2>&1
    rc=$?
    ok=1
    [ $rc -ne 0 ] && { ok=0; echo "  ✗ exit status $rc"; }

    if [ "$verdict" = "marker" ]; then
        if grep -qE '❌|FAILED|FAILS' "$log"; then
            ok=0; echo "  ✗ failure marker in output:"
            grep -nE '❌|FAILED|FAILS' "$log" | head -3 | sed 's/^/      /'
        fi
        if ! grep -q '✅' "$log"; then
            ok=0; echo "  ✗ no ✅ success marker in output (gate did not assert a pass)"
        fi
    fi

    if [ $ok -eq 1 ]; then
        pass=$((pass+1)); echo "  ✓ PASS  (log: $log)"
    else
        fail=$((fail+1)); failed_names="$failed_names $name"; echo "  ✗ FAIL  (log: $log)"
    fi
done < "$MANIFEST"

printf '\n──────────────────────────────────────\n'
printf 'tier=%s  passed=%d  failed=%d\n' "$TIER" "$pass" "$fail"
[ -n "$failed_names" ] && printf 'failed:%s\n' "$failed_names"

if [ "$TIER" = "nongating" ]; then
    echo "(nongating tier — reported only, never fails the build)"; exit 0
fi
[ $fail -eq 0 ] || exit 1
exit 0
