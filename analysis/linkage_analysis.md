# Backward/forward linkage analysis — Metals, Machinery & Electrical, Transport Equipment

Data: 2016 Indonesia national IO table underlying INDOTERM (`data/national_data.csv`, i.e.
`1BAS` domestic intermediate flows + `MAKE` output), aggregated 185 → 25 sectors using
`AGGCOM` / `SEC_MAP_185_to_25` from `src/aggregation_data.jl`. **`origin/sec.agg` (GEMPACK)
is treated as legacy and not used** — the Julia port's own aggregation mapping is the source
of truth. Method and full numbers: `analysis/linkage_analysis.jl`,
`analysis/output/linkage_25sector.csv`, `analysis/output/triad_technical_coeffs.csv`.

Rasmussen–Hirschman indices (25-sector average = 1.0): **BL** (backward linkage / "power of
dispersion", from the Leontief inverse) measures how much a sector pulls on the rest of the
economy when its output expands; **FL** (forward linkage / "sensitivity of dispersion")
measures how much the rest of the economy depends on that sector as a supplier — reported
both from the Leontief inverse (demand-side proxy) and the Ghosh inverse (supply-side,
the theoretically preferred measure for forward linkage).

## Headline result for the triad

| Sector | BL (rank/25) | FL-Ghosh (rank/25) | Classification |
|---|---|---|---|
| **BasMetals** | 1.140 (**3rd**) | 1.026 (12th) | **Key sector** (above-average both ways) |
| **MetalMach** | 0.956 (16th) | 0.863 (17th) | Weakly linked (below average both ways) |
| **TranspEquip** | 0.914 (18th) | 0.835 (19th) | Weakly linked (below average both ways) |

Only **BasMetals** clears both thresholds economy-wide. `MetalMach` and `TranspEquip` are
below the national average on *both* dimensions — despite being visually "central" to
industrial policy narratives, they are not high-multiplier sectors in the aggregate 2016
structure. This mirrors a common finding in Indonesian IO studies: metal/machinery/transport
manufacturing is dominated by assembly-type value-added-light activity with heavy reliance
on imported inputs (captured here as leakage, since only *domestic* flows are counted) and
narrow domestic onward-selling relative to sectors like Utilities, FoodProc, or Chemicals.

## Direct interconnectedness of the triad itself

Technical-coefficient sub-matrix `A[i,j]` = share of sector *j*'s total input cost sourced
directly from sector *i* (domestic flows only):

| from ＼ to | BasMetals | MetalMach | TranspEquip |
|---|---|---|---|
| **BasMetals** | 3.7% | 5.1% | 1.5% |
| **MetalMach** | 1.1% | 7.2% | 5.5% |
| **TranspEquip** | 0.1% | 0.3% | 12.1% |

Reading down the columns: **BasMetals is the metallurgical base of the chain** — it supplies
5.1% of MetalMach's total inputs (its single largest sector-to-sector metal linkage) but
almost nothing to TranspEquip directly (0.1%). **MetalMach is the true hinge of the triad**:
it draws inputs from BasMetals (5.1% of BasMetals' output goes there — see supplier list
below) *and* feeds TranspEquip directly (5.5% of TranspEquip's inputs, its largest
manufacturing input link after its own sector). **TranspEquip is highly self-contained**
(12.1% own-sector coefficient, the largest diagonal in the triad) and barely buys from
BasMetals — it is downstream of MetalMach, not of BasMetals, in this data.

So the linkage chain runs **BasMetals → MetalMach → TranspEquip**, mediated almost entirely
through MetalMach; a direct BasMetals↔TranspEquip channel is negligible.

- Share of the triad's combined intermediate **inputs** sourced from *within* the triad: **41.1%**
- Share of the triad's combined intermediate **sales** going *to* the triad itself: **35.9%**

i.e. the three sectors already form a fairly tight cluster relative to the rest of the
economy — more than a third of what they buy and sell stays inside the cluster — but the
bulk of that internal trade is not diffuse: it concentrates on the BasMetals→MetalMach and
MetalMach→TranspEquip legs.

## Who else the triad depends on / sells to

- **BasMetals'** dominant input is **OthMining** (28.7% of its cost structure — ores/basic
  minerals), well above any intra-triad supplier. Its output is overwhelmingly absorbed by
  **Construction** (58.6% of its domestic intermediate sales) and MetalMach (25.8%) — so
  BasMetals' economy-wide "key sector" status is driven mostly by construction demand, not
  the machinery/transport chain.
- **MetalMach's** own-sector input share (7.2%) exceeds any external one; Chemicals (5.7%)
  and BasMetals (5.1%) follow. Its output splits across Construction (25.1%), itself (21.1%),
  TranspEquip (10.1%) and Transport (7.8%) — a genuinely intermediate-goods sector selling
  into several downstream chains at once.
- **TranspEquip** is largely self-referential on the input side (12.1% own-sector, MetalMach
  5.5% a distant second) and sells mostly to itself (40.2%) and to Transport services (14.9%)
  — consistent with vehicle assembly/parts trade and vehicle-fleet replacement rather than
  broad economy-wide diffusion.

## Bottom line

Backward/forward linkage analysis on INDOTERM's national base data does **not** support
treating "Metals + Machinery&Electrical + Transport Equipment" as a single high-multiplier
industrial cluster in the aggregate: only BasMetals is an above-average key sector, and its
key-sector status is anchored in construction demand and mining inputs, not the
machinery/transport chain. What *is* real is a directional, moderately concentrated internal
supply chain — BasMetals feeds MetalMach, which in turn feeds TranspEquip — with MetalMach
acting as the connective sector between upstream metals and downstream transport-equipment
assembly. A policy or investment case for "clustering" these three would need to target that
BasMetals→MetalMach→TranspEquip corridor specifically (e.g. via MetalMach capacity/quality
upgrading) rather than assume broad, symmetric interdependence across all three.
