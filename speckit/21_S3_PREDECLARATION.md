# 21 — GATE S3: THE CONFIRMATORY GRID (pre-declaration, frozen 2026-08-21 before any data)

**Frozen.** This is the paper's main table. Criteria come from `13_SMALLN_PAPER_SPECKIT.md` S3 and
are restated here unchanged; what this document ADDS is the three things S3 cannot run without and
which spec 13 left open: the p>2 projection, the seed count, and the rule-fitting choice.

## 1. Configuration inherited from closed gates (not re-opened here)

- **K = 6** (36 cells) — gate S4, `20_S4_PREDECLARATION.md`. NOT the historical K = 12.
- **Deployment-matched training, mix pinned at 1.0** — gate PB/EM.
- **Randomized matched-level comparison at alpha = .05**, NA = a failed test that cannot reject —
  `probit_common.py`, shared by every gate.
- **Separation policy**: Firth refit when `any(|beta_j|>8)` or non-convergence, applied to the
  observed fit AND every bootstrap replicate. Measured rate 0.03% at n=50, 0.00% beyond.
- **B = 199**, exact rank p-value.

## 2. ⚠️ THE p>2 PROJECTION — a first-time construction, declared before use

Every map result in this project (phases 11-21, gates PB/EM/S4) is **p = 2**. The map is a K x K
grid over two covariate axes and has never been run with more. S3's design grid requires p = 3, so a
projection must be fixed. **Declaring it here, before any S3 data exists.**

**RULE: grid on the two covariates with the largest |standardised coefficient| in the fitted model,
by rank.** Standardised means `|beta_j| * sd(x_j)`, computed from the observed fit and **held fixed
for that dataset's entire bootstrap block** so the geometry cannot shift between the observed row
and its replicates — shifting it would break the exchangeability Proposition 1 rests on. A
3-level factor contributes its dummies; the factor's score is the largest |beta*sd| among them, and
if selected the grid axis is the factor's own levels (so a 3-level factor gives 3 columns, not 6).

**Why this rule.** (i) It reduces EXACTLY to the current construction when p = 2, so S3 is a strict
generalisation of everything already measured rather than a different method. (ii) It is
deterministic given the fit, scale-free through ranks, and needs no tuning. (iii) It does not use
the response beyond the fit itself. (iv) Gate EM already ruled out the obvious alternative of
projecting onto eta-hat: that geometry did not reliably help and its effect was smaller than the
corpus noise.

**What is NOT claimed:** that this is the best projection. It is a declared, defensible default.
**A p = 2 cell is included at every n specifically so the projection's cost is measurable** — if
p = 3 cells underperform their p = 2 counterparts, that is reported as the price of the projection,
not hidden.

## 3. Seeds — THREE, and this is not negotiable

Gate EM measured a learned rule swinging **.0219 - .1786 across three training corpora on one
cell** — a range LARGER than any design effect in this project. **A single-corpus S3 would be
uninterpretable.** Every cell runs 3 corpora (101000/102000/103000) and every number in the main
table is a mean across them WITH its range. No single-seed number may be quoted anywhere.

## 4. The grid

- **n in {50, 75, 100, 150, 200}**
- **Designs:** D1 textbook p=2, 50% events · D2 p=3 with a 3-level factor · D3 p=3 with AR(1)
  rho=.7 · D4 p=2 with 15% events
- Per cell: 2000 nulls (1000 train / 1000 deployed size), 1000 random departures (training alt
  class), 2000 HELD-OUT random departures with disjoint family seeds (evaluation)
- **Amplitude `2*sqrt(100/n)`** — gate S2's n=100 mid-power convention carried by spec 13's
  1/sqrt(n) rule, identical to gate S4b. Fixed here so no cell can be re-amplitude-tuned.
- 5 n x 4 designs x 3 seeds = **60 blocks**

## 5. One rule per n, pooled across designs — decided by the pre-declared 2-cell pilot

Spec 13 left this open. **Declared: one deployment-matched rule per (n, seed), trained on the D1
corpus and applied UNCHANGED to all four designs at that n.** This is the honest test of the
paper's actual claim — a *pretrained* test is one an analyst downloads and applies to their own
design without retraining. Fitting a fresh rule per design would measure something the method does
not offer. The 2-cell pilot (D1 and D3 at n=100) reports the per-design-rule numbers alongside, so
the cost of pooling is measured and reported, not assumed.

## 6. Criteria (from spec 13 S3, unchanged)

| # | Criterion | Threshold |
|---|---|---|
| **S3a** | deployed size | in [.035, .065] in EVERY cell |
| **S3b** | matched-level gain over the classical comparator | > 0 in >= 80% of cells |
| **S3c** | mean matched-level gain | >= +.03 |
| **S3d** | trend: gain at n=50 > gain at n=200 | the small-n thesis |

Comparator, fixed here: the **Cauchy combination of the eleven** (NA-robust), plus `max` over the
eleven reported alongside. Both at randomized matched level.

## 7. Fallbacks (from spec 13, restated)

- All pass -> the paper as thesised.
- Size holds but mean gain in (0,.03) -> the paper leads with CALIBRATION (T1+T2) and reports power
  parity honestly.
- Size fails anywhere -> that cell is reported and the scope shrinks to where it holds.
- **S3d fails** -> the small-n framing is retired and the paper is re-aimed at whatever n the
  evidence supports. The pivot is reported, not buried.

## 8. Forbidden

Re-running any cell with a different K, amplitude, or seed set after seeing it. Quoting a
single-seed number. Fitting a per-design rule and reporting it as the headline (it is a reported
diagnostic only). Dropping a design because it underperforms. Revising the comparator after seeing
results.

## 9. Execution

`code/s3_grid.R` -> `code/s3_score.py` -> `code/s3_verdict.py`. Deposit `results/s3/`.
2-cell pilot FIRST (D1, D3 at n=100, seed 101000), reported here, before the full grid launches.
