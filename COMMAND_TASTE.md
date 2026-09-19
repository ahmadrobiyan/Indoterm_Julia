# INDOTERM-Julia — Command Taste (6-region locked model)

> Copy-paste runbook. All paths relative to `IndotermJulia/`.
> Julia stdout is block-buffered when redirected — tail the log file, don't wait on console.

## 0. One-time setup
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 1. Fast health check (~1 min)
```bash
bash scripts/run_gates.sh fast
# fast tier (test/gates.tsv): test_syntax, load_check, verify_scenario_complete
```

## 2. Integration (~5-10 min)
```bash
bash scripts/run_gates.sh integration
# check_residual_6reg, solve_benchmark_6reg (‖F‖∞ ≈ 1.4e-9), gdp_report_6reg
```

## 3. Full V-gates (hours — run overnight)
```bash
bash scripts/run_gates.sh full
# V9 coalprice, V2 numeraire, V1 walras, V3 path, V5 multiplier,
# homogeneity, arclength fallback, swapped closure
# Verdict rule: exit 0 AND ✅ AND no ❌/FAILED (scripts exit 0 even on fail)
```

## 4. Single gate
```bash
julia --project=. test/verify_coalprice_reference.jl
julia --project=. test/verify_numeraire.jl
julia --project=. test/verify_path_independence.jl
julia --project=. test/verify_multiplier.jl
julia --project=. test/verify_arclength_fallback.jl
```

## 5. Interactive research (Julia REPL)
```julia
include("test/pipeline_cache.jl")
using IndotermJulia
agg, params = cached_pipeline(6)  # HIT in ~1s, MISS ~250s once

# externally-validated boom (V9): +50% coal export price
r = run_model!(agg, params, COALPRICE_REFERENCE; h0=0.25, tol=1e-8, verbose=true)

# read national headlines (bmk=1 indices → % change)
using IndotermJulia: MAINMACROS
nm = r.values["NatMacro"]
pct(n) = 100*(nm[findfirst(==(n), MAINMACROS)]-1)
@show pct("RealGDP"), pct("RealGNE"), pct("CPI"), pct("AggEmploy")

# regional screen (mandatory before quoting regional signs — V4 limit)
print_regional_confidence_report(r.values, benchmark_levels(params), 6)
```

## 6. Custom shock (negative coal price — contraction)
```julia
sc = Scenario(
  name="research — world coal export price -20%",
  source="constructed; mirrors COALPRICE_REFERENCE with logpct(-20)",
  swaps=COALPRICE_SWAPS,                    # coalprice closure, capital frozen
  shocks=[("fpexp_d", 5) => logpct(-20)],   # coal=sector 5, LOG-additive!
  numeraire=:exrate,                        # phi exogenous (coalprice.CMF:72 commented out)
)
r = run_model!(agg, params, sc; h0=0.25, tol=1e-8, verbose=true)
# fpexp_d is log-space: use logpct(p), NEVER pct(p) (+50% = log(1.5), not 1.5)
```

## 7. Tight stall repro (minutes, not hours)
```bash
julia --project=. test/scratch/_repro_coaldrop_stall.jl 2>&1 | tee logs/repro_$(date +%F).log
# stages t=0.015625/0.03125/0.046875 (must converge), then fires t=0.0625
# RED = STALL-ABORT / :no_progress at ~4e-7 ; GREEN = 3-5 exact-LU its
```

## 8. Row-level X-ray (names the guilty equations)
```julia
# inside any solve_newton! call:
r = solve_newton!(m, vars; verbose=true, diagnose_worst_row=true)
# prints top-5 scaled-residual rows at exit: row id, JuMP name, |F|,
# touching cols with colnrm/value/near_floor flags
```

## 9. Logs & heartbeat
```bash
ls -t logs/*2026-*.log | head -5      # newest first
tail -n 40 logs/<run>.log             # solver speaks per-iteration (verbose=ON, flushed)
tail -n 5 logs/heartbeat.jsonl        # phase heartbeats for long probes
```

## Rules that bite
- 6 regions × 25 sectors is LOCKED. 12/34 = test instruments only (needs explicit authorization).
- Levels form only. No Johansen linearisation rewrite.
- `origin/` is read-only reference. TAB/CMF win over code.
- `fpexp_d`/`fqexp_d` are LOG-additive → `logpct`. `xcap` is a FLOW → `pct_shocks`.
- Closure + numeraire belong to the Scenario. Never override at the call site.
- `u` order: hou=na+1, inv=na+2, gov=na+3, exp=na+4. Coal sector = 5.
