# 09 — RQR REBUILD: PRE-DECLARATION (frozen 2026-08-19, before any data was generated)

**Frozen.** Reported against verbatim. Implements the finding of `08_UNIVERSAL_VERDICT.md` §10.

---

## 1. What is being tested, and why

The features carry grouped **Pearson** residuals, which are only asymptotically N(0,1). Measured
directly (`code/rqr_pivotality.R`, 12 designs × 2000 nulls, no battery, no learning): their scale
collapses to sd **0.763** at n = 100 with 5% events, while **Dunn & Smyth (1996) randomized quantile
residuals hold sd 0.951 in every design**. Under one fixed cutoff, the realized-size spread was
Pearson **.0189** versus RQR **.0090**.

**Hypothesis:** the non-pivotality that defeated Routes A, B and calibrate-then-combine enters
through the *features*, and replacing the residual coding repairs it downstream.

## 2. The design — PAIRED, so the comparison cannot be confounded

Both residual codings are computed **on the same datasets, in the same pass**. Every row carries:

- 11 battery z-scores (unchanged)
- 22 **Pearson** codings: `g1..g10`, `s1..s10`, `hl`, `mabs` (unchanged, the incumbent)
- 22 **RQR** codings: `qg1..qg10`, `qs1..qs10`, `qhl`, `qmabs` — identical construction, Dunn–Smyth
  residuals substituted
- the 10 design descriptors

Two rules are then trained **on identical data**, differing only in which 22 residual columns they
receive. Any difference is attributable to the residual definition alone.

**RQR construction** (binary y, fitted p̂): `u ~ U(0, 1−p̂)` if y = 0; `u ~ U(1−p̂, 1)` if y = 1;
`r = Φ⁻¹(u)`. Grouped by fitted-risk decile as `Σr / √(group size)`. **One draw per dataset**, seed
tied to the dataset seed so the corpus is reproducible.

## 3. Grid, seeds, volume

Identical to `07_UNIVERSAL_PREDECLARATION.md` §1–2 — same training grid, same ten held-out classes
H1a–H8 — so results are directly comparable with the run already completed.
**Seed base 41000** (never used). 200k null + 200k alternative training rows, 2000 null + 2400
alternative per held-out class. Every row deposited to `results/rqr/`.

## 4. Success criteria (judged ONLY on held-out classes)

| # | Criterion | Threshold |
|---|---|---|
| **R1** | RQR size-spread across H1a–H8 (non-extrapolation) beats the Pearson rule trained on the same data | strictly smaller |
| **R2** | RQR size-spread | ≤ **.020** (Route B's was .0240; the pre-declared invariance bar) |
| **R3** | Realized size in every non-extrapolation class | ∈ [.035, .065] |
| **R4** | Mean power vs the pre-declared comparator (CCT-full, analytic null) | ≥ **−.01** |
| **R5** | Both rules obey **power ≥ size** in every alternative cell | asserted in code; a violation halts the run as a bug, not a finding |

R5 exists because a tail inversion in `calibrate_then_combine.py` produced power .028 against size
.05–.075 and was caught only because that is arithmetically impossible. The check is now automatic.

## 5. Declared fallbacks

- **R1–R4 all pass** → the residual definition was the root cause; rebuild the programme on RQR
  features and re-open the "no calibration needed" goal.
- **R1 passes, R2 fails** → RQR is a genuine improvement but not sufficient; report it as a component
  finding and keep the deployment-gap thesis as the paper.
- **R1 fails** → the pivotality gain measured at feature level does **not** survive into the trained
  rule. That is itself worth reporting: it would show the battery z-scores, not the residuals, are
  the binding constraint.
- **R4 fails while R1–R3 pass** → exactness bought at the cost of power; report the trade honestly
  and let the reader choose.

## 6. Forbidden after this document

Changing the residual construction, the grouping, the number of draws, or the criteria after seeing
held-out results. Reporting the Pearson arm only if RQR wins, or vice versa — **both arms are
reported regardless of outcome**, since they are paired by construction.

## 7. Execution

`code/rqr_features.R` (feature builder) → `code/rqr_corpus.R` (paired corpus) →
`code/rqr_train.py` (both arms, Routes A and B) → `code/universal_eval.R` (unchanged criteria).
