# INDOTERM-Julia

A Julia translation of **INDOTERM** — a TERM-family, ORANI-G-derived, multi-region
computable general equilibrium (CGE) model of Indonesia — from its original
GEMPACK/TABLO source. The port runs the full model without a GEMPACK licence.

Base data: 2016 Indonesia national input-output table plus regional supplementary
data — **185 commodities/industries × 34 provinces**, aggregated for solving to
**25 sectors × 6 island groups**.

| | |
|---|---|
| Model source (authoritative) | `origin/TERM.TAB` (equations), `origin/TERM.CMF` (closure) |
| Julia implementation | `src/` |
| Gates / verification | `test/`, manifest at `test/gates.tsv` |
| Project state | `PLAN.md` (gate table at top), `VV_PLAN.md` (V&V) |
| Working guide | `AGENTS.md` |

> **When Julia behaviour and `origin/*.TAB`/`*.CMF` disagree, the TAB/CMF files win.**
> Several defects found so far were the port silently departing from TABLO semantics.

---

## Quick start

```bash
git clone https://github.com/ahmadrobiyan/Indoterm_Julia
cd Indoterm_Julia

bash scripts/verify_data.sh                                  # provenance + unpack data
julia --project=. -e 'using Pkg; Pkg.instantiate()'          # environment
bash scripts/run_gates.sh fast                               # seconds
```

## Windows (PowerShell)

The helper scripts in `scripts/` are bash. On Windows either use **Git Bash** (ships with Git
for Windows) or **WSL** to run them:

```bash
bash scripts/verify_data.sh
bash scripts/run_gates.sh fast
```

**Individual gates need no bash at all** — run them directly in PowerShell:

```powershell
Expand-Archive -Force data\national_data.zip -DestinationPath data\
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test\sensitivity.jl        # V8 elasticity sweep
julia --project=. test\runtests.jl           # V7 regression suite
```

Gate definitions, including what each one is, are in `test/gates.tsv`.

## Reproducing published results

```bash
bash scripts/reproduce_all.sh        # all stages; expect several hours
```

Or one tier at a time:

```bash
bash scripts/run_gates.sh fast          # seconds  — parses, loads, scenario provenance
bash scripts/run_gates.sh integration   # minutes  — benchmark solve replicates
bash scripts/run_gates.sh full          # hours    — V1/V2/V3/V5/V9, homogeneity, arclength
bash scripts/run_gates.sh nongating     # known-open gates; reported, never fatal
```

Gate definitions and tiers live in **`test/gates.tsv`**. Logs land in `logs/gates/`.

CI runs `fast` on every push and `integration` + `full` + `nongating` nightly.

## Data provenance

`data/checksums.sha256` pins the committed inputs. `data/national_data.csv` is a
derived artifact (gitignored, ~66 MB) unpacked from `data/national_data.zip`;
`scripts/verify_data.sh` unpacks it and checks it against the recorded hash.

## Environment reproducibility

**Julia version.** Development and CI both target **Julia 1.12** (`.github/workflows/ci.yml`).
`Project.toml` allows >= 1.10, but the committed `Manifest.toml` records the version it was
generated on — generate it on the same version CI runs, or the lock is not actually pinning
the environment.

`Manifest.toml` **must be committed**. It is the lock file that makes a result
reconstructible. This is not a formality — an unpinned dependency update once took
`ras_balance!` from under 1 s to ~124 s (`AGENTS.md`). CI fails if it is absent.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.status()'
git add Manifest.toml && git commit -m "Pin dependency versions"
```

---

## What has and has not been established

Verification is strong; validation is partial. Both matter when reading any number
out of this model. Full detail in `VV_PLAN.md`.

| Check | State |
|---|---|
| Benchmark replication (zero shock → zero change) | ✅ ‖F‖∞ ≈ 1.4e-9 |
| Square system, full-rank Jacobian | ✅ deficiency 637 → 0 |
| Price homogeneity, degree 0 | ✅ worst cell 2.4e-15 |
| Numeraire invariance (`:gdppi` vs `:cpi`) — V2 | ✅ 0 unexplained |
| Walras identity — V1 | ✅ nationally scoped |
| Path independence — V3 | ✅ Legs A/B agree to 2.1e-7 |
| Income-side GDP decomposition — V5 | ✅ residual 0.7% |
| **External validation vs published GEMPACK result — V9** | ✅ all 8 Table 2 columns, sign + within 0.5 pp |
| Aggregation consistency — V4 | 🟤 closed as a **disclosed limitation** |
| Closure ordering (short-run leg) — V6 | 🔴 **unresolved**, time-boxed |
| Elasticity sensitivity — V8 | 🟡 implemented (`test/sensitivity.jl`), **not yet run** |

### Known limitations — read before citing

1. **Results are 25 sectors × 6 island groups, not 34 provinces.** Full 34-region
   solving needs Excerpt 49 (`Substitute`/`Backsolve`) condensation: Jacobian LU
   fill-in already exceeds 16 GB at 20 regions. Do not describe results as
   provincial.
2. **No elasticity sensitivity analysis has been run (V8).** The sweep is implemented
   (`test/sensitivity.jl`, `bash scripts/run_gates.sh full`) but has never been executed. Benchmark replication
   is structurally insensitive to exactly the parameters that drive counterfactual
   results — a mistranslated substitution elasticity changes nothing at the
   benchmark and everything off it. Magnitudes are not yet defensible against a
   sensitivity challenge.
3. **V6 (short-run closure ordering) is unresolved.** Report as "ordering not
   established for this scenario," not as a pass.
4. **V4 aggregation consistency is a disclosed limitation.** 29 of 534 region-level
   totals disagree in sign under a shocked 12→6 comparison, concentrated near
   zero-crossings. National GDP is unaffected (0.04–0.05%).
5. **The shipped base data has a regional GDP identity imbalance of up to 36.1%.**
   This is a defect in the source data, not the translation.

### Known gap: gate verdict contracts

The `verify_*` gates **print** their verdict and exit 0 even when they fail. Keying
CI on exit status alone would therefore produce a green build that certifies
nothing. As a stopgap, `scripts/run_gates.sh` applies a dual criterion to gates
marked `marker` in `test/gates.tsv`: the run must exit 0, must emit a `✅`, and must
not emit `❌`/`FAILED`/`FAILS`.

This is a workaround, not the fix. **The real fix is converting these scripts to
`Test.jl` `@testset`s that throw on failure** — tracked as remaining Phase-1 work.
Until then, treat a green `full` tier as "no gate printed a failure," which is
weaker than "every gate asserted a pass."

---

## Repository layout

```
src/                  model implementation
test/                 gates (see test/gates.tsv)
test/scratch/         one-off diagnostics — NOT gates, not run by CI
origin/               original GEMPACK/TABLO source — provenance, never modified
data/                 base-year data + checksums
scripts/              provenance, gate runner, reproduction entrypoint
analysis/             standalone analyses (e.g. backward/forward linkages)
logs/                 gate logs and session history
```

## Licence and third-party provenance

The Julia implementation is MIT-licensed — see `LICENSE`.

**`origin/` is not covered by that licence.** It contains the original INDOTERM /
TERM GEMPACK source (TERM designed by Mark Horridge). Its redistribution terms are
the original authors', and must be confirmed before this repository is made public
or archived. See `NOTICE`.

## Citation

See `CITATION.cff`.
