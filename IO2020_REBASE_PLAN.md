# Rebasing INDOTERM Julia onto the 2020 Indonesian IO table

**Status: 🟡 Step 1 (gap audit) COMPLETE and passing, 2026-08-12. Step 2 blocked on data sourcing.**

**⏸ PARKED 2026-09-12 (user decision).** Do not start Step 2 until BPS publishes the
regional-side sources that §3.D lists as the hard limit — **IRIO 2020** and **SUT 2020**.
Without them the result is a 2020-national / 2016-regional hybrid, and the decision is not to
build that. The national-side downloads (§3.B/C, Transaksi Total and purchaser-price tables)
are also on hold; nothing here is a prerequisite for the coal paper. Re-audit with
`test/audit_io2020.jl` when a new BPS release appears, then revisit.

Companion to `VV_PLAN.md`. The audit gate is [`test/audit_io2020.jl`](test/audit_io2020.jl) —
read-only, writes nothing, ~40 s. Re-run it any time:

```bash
julia --project=IndotermJulia IndotermJulia/test/audit_io2020.jl
```

---

## 1. Why

The model's entire database is 2016 vintage, inherited from the authors' shipped `national.har`
/ `regsupp.har`. BPS published the official **Tabel Input-Output Indonesia 2020** on 2025-01-27.
A 2020 baseline matters because the 2016 structure predates the nickel/bauxite downstreaming
build-out, the coal-price cycle that `coalprice.CMF` — the model's *only* external validation
scenario — is about, and COVID-year demand composition. All three sit directly on the
intermediate quadrant this model solves.

**"Using the original author approach" does not cover this step, and the plan says so.**
`origin/readme.txt` states the shipped suite generates the TERM database *from a national CGE
database already in ORANI-G format*; `origin/init.tab:1-45` merely *checks* an "initial
ORANIG03-style database NATIONAL.HAR". Everything in `origin/` — `reg0.tab` … `preagg.tab`,
driven by `mkdata.bat` — is the **regionalization** method (national ORANI-G → bottom-up regional
TERM), not the IO-table → ORANI-G construction. That upstream build is not in `origin/` at all,
and replicating the authors' route byte-for-byte needs a source-code GEMPACK licence.

The authors' approach is therefore honoured by **reusing the downstream half unchanged**: the
Julia port of `build_reg0!` → `aggregate_regions!` (~2,790 LOC) *is* their method. This project
adds only the missing upstream stage.

---

## 2. Audit results (Step 1, all four sections pass)

### §1 Classification — ✅ the single biggest risk is retired

The 2020 table's 185 products use the **same classification and the same order** as the model's
hardcoded `COM` array in `src/prepare_sets.jl:13-42`. BPS codes run 1–185 in sequence; all 10
positional anchors match:

| pos | `COM` | 2020 product name |
|---|---|---|
| 1 | Rice | Padi |
| 26 | Livestock | Ternak dan Hasil-hasilnya kecuali Susu Segar |
| 37 | CoalLignite | Batubara dan lignit |
| 38 | CrudeOil | Minyak Bumi |
| 50 | CoarseSalt | Garam Kasar |
| 96 | BasChemicals | Kimia Dasar Kecuali Pupuk |
| 111 | Glass | Kaca dan Barang-barang dari Kaca |
| 147 | WaterSupply | Pengadaan Air |
| 178 | GovEducation | Jasa Pendidikan Pemerintah |
| 185 | OthSvc | Jasa Lainnya |

This matters more than it looks. `SEC_MAP_185_to_25` (`src/aggregation_data.jl:66-134`) and
`REG_MAP_34_to_6` (`:158-`) are **positional integer arrays** and nothing in the codebase has ever
validated that position N still names the same product. A classification change would have made a
rebase hopeless. It isn't one — `COM`, `IND`, `SEC_MAP_185_to_25`, `REG_MAP_34_to_6` all carry
over **unchanged**. The audit now hard-`error`s if a future BPS release reorders anything.

### §4 Internal balance — ✅ exact

All five identities close at **0.000e+00** relative:

| check | result |
|---|---|
| intermediate column sums vs row `1900` | ✅ 0.000e+00 |
| `1900`+`1950`+`2000`+`2090` vs row `2100` | ✅ 0.000e+00 |
| intermediate row sums vs col `1800` | ✅ 0.000e+00 |
| `1800`+`3090` vs col `3100` | ✅ 0.000e+00 |
| col `8000` (supply) vs col `3100` (demand) | ✅ 0.000e+00 |

The published table is clean. Safe to build on.

### §2 Coverage — 1 supplied, 9 partial, 20 absent, of 30 headers

`src/build_reg0!.jl` consumes 30 `nat["..."]` headers. Against the 2020 domestic/basic release:

| verdict | count | what it means |
|---|---|---|
| ✅ supplied | 1 | `1PTX` (row `2030`) |
| 🟡 partial | 9 | a **control total** exists but not the interior split → these are RAS targets |
| ❌ absent | 20 | of which ~6 are elasticities (2016 carry-over **by design**, not a gap) and ~14 are genuine missing data |

Full table is printed by §2 of the audit script.

### §3 Structural change — ✅ the rebase is justified, after one correction

> **⚠ Correction recorded deliberately.** The first run reported Trade's supplier share jumping
> **+9.79 pp**, which would have been the headline. It was an artifact. ORANI-G strips
> trade/transport margins *out* of `1BAS` into a separate `1MAR` matrix; the BPS basic-price use
> table does **not** — each cell is at basic price and the margin appears as a direct purchase of
> trade/transport services by the using industry (standard SNA convention). `1MAR` is **9.3% of
> `1BAS`**. Folding it back in at `[margin commodity, using industry]` collapses the move to
> **+2.52 pp** and replaces the entire top-cell list (which had been five `Trade → *` cells) with
> economically readable ones. **Any comparison against a BPS IO table must fold `1MAR` in first.**

Largest sector-share moves on the corrected basis (shares of total domestic intermediate flow —
the two vintages are in different units and are **never** compared in levels):

| as supplier | Δpp | | as user | Δpp |
|---|---|---|---|---|
| Trade | +2.52 | | Transport | −1.27 |
| Crops | +1.15 | | OthServices | +1.07 |
| OthMining | +1.01 | | Construction | −0.82 |
| Chemicals | −1.01 | | Utilities | −0.60 |
| Transport | −0.92 | | InfoComm | +0.56 |

Largest cell-level shifts: FoodProc→HotelsRest **+0.64**, NonMetalPrd→Construction **−0.62**,
FoodProc→FoodProc **−0.50**, Chemicals→OthServices **+0.48**, Crops→FoodProc **+0.46**.

**Verdict: structure has moved materially (max 2.52 pp). The rebase is justified.**

> **⚠ Do not quote the factor shares.** 2016 labour 51.2% vs 2020 38.9% is **not** a factor-share
> change. BPS row `2020 Surplus Usaha Bruto` contains **mixed income of the self-employed**, which
> is large in Indonesia and which ORANI-G's `1LAB`/`1CAP` split allocates differently. The gap is a
> definitional difference of unknown size.

---

## 3. Files needed — with availability, checked 2026-08-12

### A. In hand

| # | File | Supplies |
|---|---|---|
| 1 | **IO 2020 Transaksi Domestik, Harga Dasar, 185 Produk** — local CSV, [BPS table 2270](https://www.bps.go.id/en/statistics-table/1/MjI3MCMx/tabel-input-output-indonesia-transaksi-domestik-atas-dasar-harga-dasar--185-produk---2020--juta-rupiah-.html) | domestic intermediate quadrant + final-demand/supply controls |
| 2 | `data/national_data.csv` (71 MB) | the 2016 ORANI-G database — **structural donor** for everything 2020 lacks |
| 3 | `data/regsupp_data.csv`, `data/distgone_data.csv` | 2016 regional trade/distance/margin structure |

### B. To download — ✅ available now

| # | File | Unblocks | Link |
|---|---|---|---|
| 4 | **IO 2020 Transaksi Total, Harga Dasar, 185 Produk** — *the critical one* | imports by subtraction (file 4 − file 1) → `src=imp` halves of `1BAS`,`2BAS`,`3BAS`,`5BAS`,`6BAS`; without it there is no Armington structure | [BPS table 2272](https://www.bps.go.id/id/statistics-table/1/MjI3MiMx/tabel-input-output-indonesia-transaksi-total-atas-dasar-harga-dasar--185-produk---2020--juta-rupiah-.html) |
| 5 | **Tabel Input-Output Indonesia 2020** (publication, 2025-01-27) | methodology, concordances, definitions — read before trusting any subtraction | [BPS publication](https://www.bps.go.id/id/publication/2025/01/27/ea74e60a8662a70a734fa431/tabel-input-output-indonesia-2020.html) |

### C. To download — 🟡 very likely available, **verify by hand**

| # | File | Unblocks | Note |
|---|---|---|---|
| 6 | **IO 2020 Transaksi Total, Harga Pembeli, 185 Produk** | purchaser − basic = margin+tax wedge per cell → `1MAR`…`5MAR`, `1TAX`…`5TAX`, `0TAR` | BPS states the 2020 IO is presented at **both** basic and purchaser prices. The 17-product purchaser version is confirmed at [table 2273](https://www.bps.go.id/id/statistics-table/1/MjI3MyMx/tabel-input-output-indonesia-transaksi-total-atas-dasar-harga-pembeli--17-produk---2020--juta-rupiah-.html) and the 185-product purchaser version exists for 2016 at [table 2118](https://www.bps.go.id/en/statistics-table/1/MjExOCMx/tabel-input-output-indonesia-transaksi-total-atas-dasar-harga-pembeli--185-produk---2016---juta-rupiah-.html). The 2020/185 variant is almost certainly table **2274** (`MjI3NCMx`) but automated checking is blocked (BPS returns HTTP 403), so **confirm in a browser**. |

### D. ❌ NOT published for base-year 2020 — this is the hard limit

| # | File | Needed for | Status |
|---|---|---|---|
| 7 | **IRIO 2020** (inter-regional IO, 34 provinces) | `R001`–`R005` (185×34) — the **entire regional dimension** | **Not published.** Only IRIO **2016** exists: [34 prov × 52 industries](https://www.bps.go.id/id/statistics-table/1/MjEyOCMx/tabel-inter-regional-input-output-indonesia-transaksi-domestik-atas-dasar-harga-produsen-menurut-34-provinsi-dan-52-industri--2016--juta-rupiah-.html), [publication](https://www.bps.go.id/en/publication/2021/12/29/3ea49c0d856eceaba836792d/tabel-interregional-input-output-indonesia-tahun-2016-tahun-anggaran-2021.html). IRIO 2016 took a 3-year compilation cycle; a 2020 edition is not in the public domain. |
| 8 | **SUT Indonesia 2020**, supply side | `MAKE` (185×185) | **Not published.** The latest confirmed SUT is **SUT 2016**. Without it `MAKE` stays 2016, which effectively freezes the multi-product/CET block (`SCET`). |

### E. Partial substitutes for D — obtainable, lower fidelity

| # | File | Substitutes for |
|---|---|---|
| 9 | **PDRB Provinsi menurut Lapangan Usaha 2020** + **PDRB Pengeluaran Provinsi 2020** | rebuild `R001`–`R005` by concordance from provincial GRDP instead of IRIO. 17-sector granularity vs the model's 185 — needs the 2016 within-sector shares as the splitter. |
| 10 | **Sakernas 2020** industry × occupation cross-tab | `1LAB`'s 4-way `OCC` split (row `2010` is a single total) |
| 11 | **PMTB menurut jenis barang modal dan lapangan usaha 2020** | `2BAS` — splitting the single `3030` column by investing industry |
| 12 | **National-accounts detail 2020** (mixed-income decomposition) | splitting row `2020` into `1CAP` / `1LND` / mixed income |

### F. 2016 carry-over by design — **not a gap**

All elasticities and behavioural parameters: `SLAB`, `P028`, `P015`(`1ARM`), `SMAR`, `SCET`,
`P018`, `ALFA`, `QRAT`, `TARG`, `DPRC`, `TFRO`, `RADJ`, `REXP`, `ENGL`, `ITEX`, `TAU`, `XPEL`.
These are estimated/assumed, never observed in an IO table. The authors' values stand.

---

## 4. What is actually achievable

| Block | Achievable vintage | Why |
|---|---|---|
| Domestic intermediate quadrant | **2020, observed** | file 1 |
| Import split (`src=imp`) | **2020, observed** | file 4, by subtraction |
| Final demand column totals | **2020, observed** | file 1 |
| Margins + taxes per cell | **2020, estimated** | file 6 gives the wedge; the `1MAR`/`1TAX` split uses 2016 proportions, RAS'd to 2020 totals |
| `1LAB` occupation split, `1CAP`/`1LND` | **2020 totals, 2016 proportions** | files 10/12 improve this |
| `MAKE` | **2016** | SUT 2020 not published |
| `R001`–`R005`, regional trade | **2016** | IRIO 2020 not published |
| Elasticities | **2016** | by design |

**The honest conclusion: the national block can go to 2020; the regional block cannot.** For a
model whose entire purpose is regional detail, that must be disclosed, not glossed. This is a
**2020-national / 2016-regional hybrid**, and every result quoted from it names that — the same
convention as always naming the closure.

---

## 5. Code files to create / modify

| Path | Action | Note |
|---|---|---|
| `test/audit_io2020.jl` | ✅ **created, run, passing** | the Step-1 gate |
| `data/io2020/` | create | 2020 raw inputs, kept out of `data/` so the validated 2016 database is never shadowed |
| `src/read_io2020.jl` | create | parse the BPS CSVs (UTF-8 BOM, thousands separators, `-` = **zero not minus**) into the `Dict{String,Array}` contract `read_data.jl` returns. The quoted-CSV splitter and `bps_number` in `test/audit_io2020.jl` are already written and tested — **lift them, don't rewrite** |
| `src/build_national2020!.jl` | create | the new upstream stage: assemble an ORANI-G-shaped `nat` dict from B/C, using A as controls and 2016 as the structural donor |
| `src/read_data.jl` | modify | add `vintage` argument (`:y2016` default, `:y2020`). 117 LOC, nothing downstream changes |
| `test/pipeline_cache.jl` | modify | thread `vintage`; key the cache on it (`data/cache_6reg_2020.jls`) so both vintages coexist |
| `VV_PLAN.md`, `AGENTS.md` | modify | record the vintage composition of every block |

**Do not touch** `build_reg0!` → `aggregate_regions!` (Steps 0–4), `prepare_sets.jl`,
`aggregation_data.jl`, or anything in `origin/`. That the authors' pipeline runs **unmodified**
on the new database *is* the fidelity test.

**Reuse, do not reimplement:** `src/ras_balance!.jl` (270 LOC) and `src/build_pstras!.jl` (230
LOC) already implement proportional fill + balancing. Every estimated block is RAS'd to a
published 2020 control total, never left free. Do not write a second RAS.

---

## 6. Verification

1. **Step 1 (done)** — classification assertion passes; coverage accounts for all 30 headers;
   the 2020 table balances at 0.000e+00.
2. **Post-build, pre-solve** — the `1BAS` domestic half reproduces the raw CSV **cell-for-cell**
   (it is observed, not estimated; any deviation is a loader bug). Every RAS-filled block
   reproduces its published 2020 control total to ~1e-8.
3. **Benchmark solve** — zero-shock `run_model!` reproduces the 2020 benchmark to the residual the
   2016 database achieves (project record: worst 1.8e-4). Materially worse ⇒ the assembled
   database is not internally consistent, *not* a solver regression.
4. **Re-gate on 2020** — `test/verify_walras.jl` (V1), `verify_homogeneity.jl` (V2),
   `verify_path_independence.jl` (V3), `verify_numeraire.jl` (V5). Same scripts, new database, no
   new test code. **V2 is the sharpest test** of a hand-assembled database: it fails loudly if any
   price/quantity split was fabricated inconsistently.
5. **Both vintages coexist** — `cached_pipeline(6)` (2016) still passes every gate after the
   `vintage` refactor. Regression-protects the validated 2016 result.

## 7. Not in scope

- Changing the 6-region or 25-sector locks — both carry over unchanged.
- Re-estimating elasticities (§3.F stays 2016).
- Reaching 34 regions — still gated on the Excerpt 49 condensation work, unrelated to vintage.
