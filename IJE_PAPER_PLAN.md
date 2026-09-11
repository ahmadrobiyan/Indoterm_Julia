# IJE Paper Plan — regional incidence of a coal-price shock in INDOTERM-Julia

Target: **Indonesian Journal of Energy** (Purnomo Yusgiantoro Center; ISSN
2549-760X e / 2549-1016 print; DOI 10.33116/ije; 2 issues/year, Feb & Aug;
double peer-reviewed; **no submission fee**; English; `.docx` template;
APA references). Scope explicitly covers *energy economics, energy analysis,
energy modeling and prediction* — a CGE incidence paper fits without stretching.
Contact: ije@pycenter.org.

## 1. Proposed paper

**Working title:** *Who gains from a coal boom? Regional incidence of a world
coal-price shock across Indonesian islands in a validated TERM-family CGE model*

**Research question:** a +50% world coal export price moves national GDP barely
(−0.2%) — but which islands gain, which lose, and how robust are those signs?

**Core result (already in hand, validated):**
- National: Real GDP −0.205%, CPI +1.42%, exports −3.88%, coal output +2.35%
  (V9: all 8 macro columns match published GEMPACK results in sign, ≤0.5pp).
- Regional: **Kalimantan GDP +5.2% at both 6- and 12-region resolutions**;
  every other island <0.4%; Sumatra/Kalimantan GDP signs resolution-robust (V4R6).
- Every headline number carries an elasticity-driven interval (V8), e.g. coal
  employment +14.5%…+23.5% — quoted as ranges, never points.

**Novelty beyond replication:** the published reference (Horridge et al.
draftreport) reports national + 18-region rows; this paper adds (a) island-level
incidence with a cross-resolution robustness test almost no regional CGE paper
runs, (b) sensitivity intervals on every headline, (c) a Julia/levels
re-implementation with a public verification record (VV_PLAN.md).

## 2. New work required (in order)

1. **Emissions/subsidy check — DONE 2026-09-11, negative.** Zero hits for
   `emission|co2|carbon|subsid` across `src/` and `origin/TERM.TAB`. Tax-revenue
   `(change)` variables exist (`delTAXhou/inv/gov/exp` by commodity×source×
   region, TERM.TAB:273-299) but no source `.CMF` ever shocks them and none is
   validated as an instrument. **Decision: incidence paper, not subsidy/carbon.**
   A fuel-subsidy or carbon extension would need new-instrument development plus
   its own validation — weeks, and a second paper, not this one.
2. **One new scenario** for novelty beyond replication (replication alone will
   not clear review). Candidate: combined coal+oil/gas price shock, or a
   counterfactual export-restriction sketch. Must run under
   `COALPRICE_REFERENCE`-class closure; each scenario ≈ 2–5 min at 6 regions.
3. **Regional confidence screen on every reported scenario**
   (`print_regional_confidence_report`) — near-zero regional signs excluded
   from directional claims per V4.
4. **Write to the IJE `.docx` template**, APA references, English.

## 3. Mandatory disclosures (reviewer-proofing section)

- 2016 database (pre-nickel-downstreaming, pre-coal-cycle) — retrospective
  validation study, not a current-market forecast.
- 6 island groups, not 34 provinces — provincial claims out of scope.
- V4: near-zero regional flow signs (notably BaliNusa) are
  resolution-sensitive; excluded from directional claims.
- V8: all ranges are lower bounds (1 of 13 points unattainable); SCET inert.
- V6 closure ordering not established for this scenario class.

## 4. What NOT to promise

- Provincial (34-region) results — needs condensation work, a separate project.
- Carbon/emissions analysis — unless step 1 finds the variables.
- Point estimates — every number ships with its V8 interval.

## 5. Schedule sketch

- Steps 1–2 (model checks + new scenario): ~1–2 weeks elapsed (minutes of
  compute, hours of interpretation).
- Writing + template: 2–3 weeks.
- Target an August or February issue per the biannual cycle; submission itself
  is free.
