# 07 — UNIVERSAL TEST: PRE-DECLARATION (frozen 2026-08-19, before any data was generated)

**This document is frozen.** Results are reported against it verbatim. Written before the corpus
generator was run. Implements `06_UNIVERSAL_SPECKIT.md`.

---

## 1. The training grid (what the model may see)

Seed base **31000** (never used in this project). Every cell deposits its feature row, its design
descriptor vector, and its label.

| Axis | Training values |
|---|---|
| n | 100, 150, 250, 400, 600, 900, 1400 |
| p (covariates) | 2, 3, 5, 8 |
| Covariate law | `indep` (mixed normal/uniform), `ar1` with ρ = 0.5 |
| Categorical | none, or one 3-level factor |
| Target event rate | 0.50, 0.15 |
| Departures | random families from `make_family_p()` (the p-generalized generator), regenerated per cell |

## 2. The held-out grid (never trained on; the only grid the claims are judged on)

| Class | Definition | Why |
|---|---|---|
| **H1 n-interp** | n ∈ {200, 500, 1000}, otherwise training values | keeps C1 / G-BOOT comparable |
| **H2 n-extrap** | n = 2000 | outside the training range |
| **H3 p-interp** | p = 4 | between trained values |
| **H4 p-extrap** | p = 12 | outside the trained range |
| **H5 dependence** | AR(1) ρ = 0.8 | stronger correlation than trained |
| **H6 skew** | one log-normal covariate | covariate law never trained |
| **H7 categorical** | one 5-level factor | more levels than trained |
| **H8 sparse** | target event rate 0.05 | the sparse regime this family cares about |

Each held-out class is evaluated with **2000 null datasets** (size) and **300 datasets per
confirmatory family** (power), at n = 400 unless the class fixes n.

## 3. The rules under test — fixed now, no changes after this document

**Shared score features (33):** 11 battery z-scores (`FEATS11`) + grouped residuals dual-coded
(`g1..g10`, `s1..s10`, `hl`, `mabs`). G = 10 regardless of n.

**Design descriptors d (10), computable from the user's own fit, never from truth:**
`n`, `p`, `mean(π̂)`, `sd(π̂)`, `min(π̂)`, `max(π̂)`, observed event rate, `log cond(X'WX)`,
`max leverage`, share of binary columns in X.

| Route | Rule |
|---|---|
| **A** | 20 nnets, size = 8, decay = .05, seeds 5001–5020, trained on the wide corpus. Threshold = the pooled-null 95th percentile, shipped as a single number. |
| **B** | Route A's score network unchanged + a calibration MLP (3 hidden layers, 128-128-64, ReLU) mapping d → score quantiles on the τ-grid {.80,.85,.90,.925,.95,.96,.97,.975,.98,.99}, pinball loss weighted ×3 for τ ≥ .95, monotone in τ enforced by cumulative-softplus output. Trained on **null datasets only**. p-value by interpolation of F̂₀(·|d). |
| **C** | Same score network trained jointly against an adversary (MLP 64-64) predicting d from the score on **null data only**, gradient reversal, λ ∈ {0, 0.1, 0.3, 1, 3, 10}. λ = 0 is Route A by construction. |

**The score network never receives d.** This is a hard constraint, not a preference — the leakage
autopsy (`MANIFEST.md` 1.2) showed these models read design fingerprints when allowed to.

## 4. Success criteria — judged ONLY on the held-out classes of §2

| # | Criterion | Threshold |
|---|---|---|
| **P1** | Realized size at the shipped rule, in **every** held-out class | ∈ [.035, .065], Clopper–Pearson covering .05 |
| **P2** | Pooled held-out null p-values | KS vs U(0,1), p > .01 |
| **P3** | **Invariance**: max − min realized size across H1–H8 | ≤ **.02** (current state: .033) |
| **P4** | Mean power over the 12 confirmatory families vs the pre-declared comparator (CCT-full), each at its own realized size | ≥ **−.01** (a tie is a win — the deliverable is *no calibration needed*) |
| **P5** | Extrapolation classes H2 and H4 | size ∈ [.03, .07], **reported separately, never pooled into P1** |

## 5. Declared fallbacks

- **All pass (any route)** → the tool as originally envisioned; venue question reopens.
- **P1–P3 pass, P4 fails** → a calibration-free test that ties classical power. Ship it; state plainly
  that the power advantage did not survive.
- **B passes, C fails** → "learning the threshold works; forcing invariance does not."
- **All routes fail P1–P3** → the deployment gap is structural. `05_GBOOT_VERDICT.md`'s warning
  becomes the paper, strengthened by a serious failed attempt to repair it.
  **This is a legitimate outcome and will not be described as failure.**

## 6. Forbidden after this document

Tuning λ, the descriptor set, the corpus grid, or the τ-grid after seeing held-out results.
Reporting in-grid coverage as if it were the claim. Dropping a held-out class that fails.
Re-running with new seeds because the first answer was unwelcome.

## 7. Execution

`code/universal_corpus.R` (generator) → `code/universal_train.py` (Routes A/B/C) →
`code/universal_eval.R` (held-out evaluation). Everything deposited under
`results/universal/`. Uses `code/fast_battery.R`, whose output is verified identical to the
original battery (80 datasets, 0 mismatches, max |diff| = 0.000e+00).
