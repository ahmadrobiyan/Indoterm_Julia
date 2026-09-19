# Full gate tier driver (PowerShell equivalent of scripts/run_gates.sh full).
# Runs the 9 full-tier gates sequentially; each gate appends to the shared log.
# Verdict contract per test/gates.tsv: exit code + (for marker gates) require ✅, forbid ❌/FAILED/FAILS.
# Launch detached: Start-Process powershell -ArgumentList "-NoProfile","-File","test/scratch/_run_full_gates.ps1"
$ErrorActionPreference = "Continue"
$root = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent  # test/scratch -> repo root
Set-Location $root
$log = "logs/gates_full_2026-09-14.log"
"===== FULL TIER START $(Get-Date -Format o) =====" | Tee-Object $log
$gates = @(
  "test/runtests.jl",
  "test/verify_coalprice_reference.jl",
  "test/verify_numeraire.jl",
  "test/verify_walras_gdp.jl",
  "test/verify_path_independence.jl",
  "test/verify_multiplier.jl",
  "test/verify_homogeneity_tol.jl",
  "test/verify_arclength_fallback.jl",
  "test/verify_swapped_closure.jl"
)
$summary = @()
foreach ($g in $gates) {
  "----- GATE $g START $(Get-Date -Format o) -----" | Tee-Object $log -Append
  $t = Measure-Command { julia --project=. $g >> $log 2>&1 }
  $code = $LASTEXITCODE
  "----- GATE $g EXIT $code  elapsed $([int]$t.TotalSeconds)s -----" | Tee-Object $log -Append
  $summary += [pscustomobject]@{gate=$g; exit=$code; secs=[int]$t.TotalSeconds}
}
"===== SUMMARY =====" | Tee-Object $log -Append
$summary | Format-Table -AutoSize | Out-String | Tee-Object $log -Append
"===== MARKER SCAN (fail words) =====" | Tee-Object $log -Append
Select-String -Path $log -Pattern "❌|FAILED|FAILS" | Select-Object -First 20 | Tee-Object $log -Append
"===== FULL TIER END $(Get-Date -Format o) =====" | Tee-Object $log -Append
