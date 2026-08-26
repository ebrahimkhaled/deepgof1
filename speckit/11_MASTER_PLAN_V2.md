# 11 — MASTER PLAN v2: the mechanism, the fix, and what remains

**Written 2026-08-20**, after the residual map beat the classical rule under oracle calibration,
collapsed under deployment, and was then substantially repaired by adversarial pivoting.
Supersedes `04_MASTER_PLAN.md` wherever they disagree.

---

## §1. What we now know that we did not know a week ago

### 1.1 The mechanism (measured, phase 12)

A classifier separates **bootstrap-generated null data** from **real null data** at **AUC .62**,
while their marginal statistics differ by **under 1%** (mean |cell| .3139 vs .3111; sd .3967 vs
.3948).

> **A smooth classical statistic averages that discrepancy away. A trained network is built to find
> whatever separates two distributions, so it finds this one too — and its observed score lands at
> the .006 quantile of its own bootstrap distribution.**

This explains every collapse in the project, including why the *best* statistic collapsed *worst*:
power and bootstrap-validity both come from the same sensitivity.

| rule | oracle | deployed | change |
|---|---|---|---|
| 20-net ensemble (battery z + grouped residuals) | .4093 | .2954 | −.114 |
| **12×12 residual-map CNN** | **.6315** | **.2830** | **−.349** |
| classical Cauchy combination | .2912 | .3103 | +.019 |

### 1.2 The diagnostic (this may outlast the test itself)

**One number predicts the failure**: the AUC at which a rule's own score separates bootstrap-null
from true-null data. Drive it toward .50 and level *and* power recover together.

| score-AUC | deployed size | deployed power |
|---|---|---|
| .622 | .006 | .296 |
| .580 | .012 | .350 |
| .550 | .032 | .390 |
| .529 | .046 | .434 |
| **.518** | **.048** | **.442** |

It requires no knowledge of the truth and is not specific to our method. **Any learned test
statistic can be audited this way before deployment.**

### 1.3 The fix that works

Adversarial pivoting (Louppe–Kagan–Cranmer, arXiv:1611.01046) aimed at **the bootstrap-vs-truth
discrepancy** — not at sample size or dimension, which was the wrong target.
`code/adversarial_bootstrap.py`. Result at λ=50: **size .0060 → .0480, power .2830 → .4420**,
gap to classical **−.207 → −.050**.

### 1.4 A defect the author's question exposed

Map cells were `Σr/√(count)`. Under H0, Var(Σr) = **Σπ̂(1−π̂)**, not the count. Measured
consequence: per-cell SD ranged **.152–.494 (3.26× spread)**; with the correct denominator,
**.845–1.072 (1.27×)**. The scale pattern is a function of π̂ — **exactly what differs between
bootstrap and real data** — so part of the sensitivity was built into the feature.
Corrected in `map_corpus.R` / `map_deployed.R`; the v2 rebuild tests it.

---

## §2. The construction, stated for the paper

**Offline, once:** simulate null and misspecified datasets across the design space; for each, form
the **residual map** — a K×K grid over the *ranks* of two covariates, cell value
`Σ(y−π̂)/√(Σπ̂(1−π̂))`, which is ~N(0,1) per cell and scale-free by construction. Train a CNN to
separate null from misspecified, **with an adversarial penalty forcing the score to be blind to the
bootstrap-vs-truth distinction**, and optionally a distillation head predicting what the classical
battery would say (teacher at training time, never an input).

**At use time:** compute the map, score it with the frozen network, calibrate by per-dataset
parametric bootstrap, report the exact rank p-value.

**Why each piece is there, with the measurement that justifies it:**

| piece | why | evidence |
|---|---|---|
| map over covariate ranks, not fitted risk | misfit lives in covariate space; a risk-ordered summary smears a local bump across every bin | phase 1.3 (−.146); map .632 vs summary .392 |
| ranks, not raw values | scale- and distribution-free input | the axis that broke transfer at H5/H6/H7 |
| `√(Σπ̂(1−π̂))` denominator | makes cells ~N(0,1) and removes a π̂-dependent scale pattern | 3.26× → 1.27× SD spread |
| convolution | finds a local departure anywhere in the frame | the property a 1-D summary cannot have |
| adversarial penalty | the score must not see the bootstrap artifact, or its calibration breaks | AUC .62 → .518, size .006 → .048 |
| distillation | the battery's knowledge without its miscalibration | +.013 |

---

## §3. Gates remaining, with pre-declared criteria

### M3 — v2 (RUNNING): does correct standardisation reduce the price of validity?
**Criterion.** At the λ giving size ∈ [.035,.065], v2 power > v1's .4420.
**Fallback.** If not, the standardisation defect was not the binding constraint; report both and
keep v1's numbers as the headline.

### M4 — multi-design transfer (REQUIRED before any claim)
Everything above is **one design** (n=400, p=2). The residual channel's calibration advantage was
measured across ten designs; the map's power was not.
**Criterion.** Size ∈ [.035,.065] in every held-out design class, and positive deployed gain in at
least 5 of 7. **Fallback:** if the map works only at p=2, the contribution narrows to
*a construction for low-dimensional GOF* — still publishable, stated plainly.
⚠️ **Open design problem:** the map is defined for two covariates. For general p the candidates are
2-D grids over leading principal components, p channels of (x_j, fitted risk), or a
permutation-invariant set network. **Decide only after M3.**

### M5 — sample-size transfer
n ∈ {100, 200, 1000, 2000} with the rule trained at n=400 only.

### M6 — real data
Pima2 and LBW through the full pipeline, reported whatever they show.

### M7 — the honest comparison table
Every rule ever built, oracle vs deployed, on identical replicates. **This table is the paper**,
regardless of which row wins.

---

## §4. The three papers this can become

1. **Positive** (if M3–M5 pass): *a pretrained goodness-of-fit test that is correctly calibrated and
   more powerful than the classical benchmark under deployment* — with the adversarial construction
   as the enabling device.
2. **Quantitative** (if the gap plateaus): *the price of validity for a learned test statistic* —
   the curve of power surrendered against calibration recovered, plus the AUC diagnostic.
3. **Cautionary** (if the map fails to transfer): *learned test statistics collapse under the only
   calibration analysts can run, and here is the mechanism, the diagnostic, and the partial repair.*

**All three are publishable. Do not force outcome 1.** The evidence base — 1.5M+ simulated
datasets, every p-value deposited, criteria frozen before each run — supports whichever is true.

**Venue:** JCGS remains the floor and the likely home for 2 or 3. Outcome 1 with M4–M5 passing
would justify reconsidering a reach.

---

## §5. Standing rules (violations have already cost this project twice)

1. **Never report an oracle-calibrated number as a result.** It manufactured +.118 (retracted) and
   +.159 (M1, now known to collapse to −.207).
2. **Both arms get the same calibration device**, on the same replicates.
3. **Assert power ≥ size** in code; a violation is a bug, not a finding. (Caught a tail inversion.)
4. **Report size per design class, never pooled** — pooling passed a rule that was wrong everywhere.
5. **Judge spread against the noise floor** (~.014 at 2000 nulls/cell); differences below ~.015 are
   not interpretable. State the floor beside every spread.
6. **Deposit every per-replicate p-value.** It has already saved an 11-hour rerun.
7. **Every trial enters the manifest**, including failures, with its numbers and its verdict.

---

## §6. Figures (deferred deliberately, now unblocked)

Held back because every headline number moved at least once — `fig_money.png` still shows retracted
values. The findings now stable enough to draw:

1. **The AUC diagnostic curve** — score-AUC against deployed size and power. The night's best result
   and a tool others can use.
2. **Oracle vs deployed, every rule** — the collapse, with the classical rule flat.
3. **The residual map itself** — a null map beside a misspecified one. Shows in one glance why
   covariate space beats risk ordering.
4. **The price-of-validity frontier** — power against λ, with the valid-size region shaded.
