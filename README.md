# DeepGOF-1

A pretrained goodness-of-fit test for logistic regression. The test statistic is a small
convolutional network (18,273 frozen parameters) that reads the model's *residual map* — a 6×6
grid of standardized residual sums over the ranks of the two strongest covariates — and the
p-value is a rank inside the analyst's own parametric bootstrap, so **the level is a property of
the calibration, not of what the network learned**.

This is a **versioned first model**. Everything needed to train a better one — the corpus
generator, the training code, the benchmark harness, and every per-replicate p-value behind the
paper's tables — is in this tree, together with the frozen pre-declared criteria we were judged
by, so a **DeepGOF-2 is judged the same way**.

---

## 0. Two ways in

**Just want to run the test?** Install the R package -- the deployable test ships as
`deepgof1()` in `ebrahim.gof` (>= 2.6.0), with the frozen weights inside:

```r
install.packages("ebrahim.gof")     # or remotes::install_github("ebrahimkhaled/ebrahim.gof")
library(ebrahim.gof)
deepgof1(glm(y ~ x1 + x2, data = d, family = binomial()))
```

**Want to check the paper, or train a successor?** That is what this archive is for: the
training corpus and its generator, the benchmark harness, every per-replicate p-value
behind every table, and the frozen pre-declared criteria. Start at Section 3.

## 1. How do I use the test from this archive?

Pure base R. No Python, no GPU, no training, no tuning. One test ≈ 2–5 s.

```r
source("R/deepgof1.R")
fit <- glm(y ~ x1 + x2, data = d, family = binomial())
deepgof1(fit)
#>  DeepGOF-1: pretrained residual-map goodness-of-fit test
#>  S = 8.2852, B = 199, p-value = 0.0050
#>  grid: 6 x 6 over ranks of x1 and x2
```

The p-value is `(1 + #{bootstrap scores >= observed}) / (B + 1)` with `B = 199`; a bootstrap
refit that fails scores `+Inf`, which counts *against* rejection (the conservative direction).
Run `Rscript release/smoke_test.R` (one level up) for a working example, and
`Rscript training/VERIFY.R` to confirm the base-R forward pass matches PyTorch on the shipped
weights (**agreement 4.6e-11**).

## 2. How was it trained?

Once, offline, on simulated data. Two decisions carry the whole design:

- **Deployment-matched nulls.** The null class is *bootstrap replicates* — simulate from a
  fitted model, refit, map — because that is exactly the object the deployed test compares
  against. Training against plain correct-model data broke calibration in our sweeps: the score
  learns to see the bootstrap itself (paper, Fig. 4; the score-AUC diagnostic).
- **A versioned departure prior.** Alternatives come from ~100 random departure *families*
  (quadratic, cubic, interaction, sine, bump, step, |x|, and index-type terms), amplitudes
  scaled on the linear predictor. The prior is a first-class object: the same architecture
  scored **+.090** under one prior and **−.071** under another. Ours is deposited verbatim as
  `training/ship_corpus.csv` (8,400 labeled maps) and defined by `training/01_corpus.R`.

Architecture: `1×6×6 → Conv3×3(16) → Conv3×3(32) → MaxPool2 → Conv3×3(32) →
GlobalMax⊕GlobalMean(64) → FC64 → FC1`, 25 epochs, Adam, BCE (`training/02_train.py`).
Export to the R deployment format: `training/03_export_to_R.py` (read its port-trap note).

## 3. How do I make DeepGOF-2?

1. Edit **`training/01_corpus.R`** — the one file a successor edits (change the departure
   prior, the null class, the sample size, or the resolution `K`).
2. `python training/02_train.py` on the new corpus; export with `training/03_export_to_R.py`.
3. `Rscript training/VERIFY.R` — must print max |difference| < 1e-6, or do not deploy.
4. Run the benchmark (`benchmark/`, see `benchmark/PORTING.md`) and **compare against the same
   frozen criteria we did**: the pre-declarations in `speckit/` ship with this release exactly
   so a successor cannot move the goalposts.

The five measured directions most likely to pay (each with the deposit behind it, none
guaranteed — one already failed once, as recorded):

1. **Representation, not learner:** p = 3 designs lose 3–4× what p = 2 does because the 2-D
   grid discards a covariate. A permutation-invariant set network over (xᵢ, π̂ᵢ, rᵢ) removes the
   projection entirely.
2. **Corpus-bagging:** training-corpus variance dominated design comparisons (one cell moved
   .0219→.1786 across three corpora). Bagging over M corpora shrinks it like √M.
3. **Selection over the classical battery:** the Cauchy combination loses ~.057 to its own best
   member. Our gate SEL attempt failed (−.0354) for a diagnosed reason — the selector predicted
   which arm fires most extremely, not which has most matched-level power. The objective is
   fixable; the gap is still there.
4. **Conditioning on the sufficient statistic:** our probe was ruled an *invalid measurement*
   (the AUC diagnostic had not converged: .509→.607 from 1.5k to 30k maps). Re-run at N ≈ 30,000
   before drawing any conclusion.
5. **Tune the departure prior alone** — no architecture change needed (see the ±.09/−.07 fact
   above).

## 4. What are the known weaknesses? (verbatim, with numbers)

- **Power rank 7th of 13** on the reproduced Liu et al. (2024) benchmark at matched level
  (mean .637). The test **beats every partition-based test at every n** (+.13/+.15/+.06 vs
  Hosmer–Lemeshow at n = 50/100/500) but loses to the covariate-space and smoothing tests
  (proj .769, le-Cessie .756). It is a **small-sample instrument**: its rank degrades
  5th → 5th → 8th as n grows 50 → 100 → 500.
- **The 2-D map discards covariates:** against the best PUBLISHED test per cell, p = 2 designs
  are at parity (+.009, +.003) while p = 3 designs pay −.060 and −.082 (gate S3) — a
  representation limit, formalized in the paper's Corollary S.6.
- **A thin blind cone exists:** directions with negative consistency certificate
  (min found −0.281), all sub-cell-resolution oscillations — the boundary every binned
  statistic shares. Spherical measure < 1.2e-5 (95% one-sided, 250k directions).
- **Probit-type alternatives:** the map signal μ is tiny (symmetric, tails-only), so finite-n
  power is near α for every grouped test including this one; consistency still holds but with a
  poor constant.
- Size is the strong suit: **58/60 pre-declared design cells in the nominal band**, and on the
  external benchmark the level is the most STABLE of the 13 arms — mean .0506 at nominal .05,
  with a block-to-block SD of .004 against .012 (le-Cessie) and .011 (proj). Those three arms
  tie on the mean; they differ in dispersion, and that is the claim that survives.

## Layout

```
R/            deepgof1() + the base-R forward pass (no Python at deploy time)
inst/extdata/ frozen weights, scaler, and the cross-language verify bench
training/     01_corpus.R · 02_train.py · 03_export_to_R.py · VERIFY.R · ship_corpus.csv
benchmark/    liu/ (Liu et al. 2024 settings, 13 arms) · scenarios/ (9×5×3 lattice)
              · grid/ (60-cell design grid) · theory/ (certificate scripts) · harness/
              · PORTING.md
results/      every per-replicate p-value behind the paper's tables, as deposited
speckit/      the frozen pre-declarations for every gate (the goalposts, shipped)
```

## Licence

Code: MIT. Deposited results and corpora: CC-BY 4.0. See `LICENSE.md`.

## Citing this archive

Archived at Zenodo. Two DOIs, and which one you want depends on why you are citing:

- **`10.5281/zenodo.22113219`** — the **concept** DOI. Always resolves to the newest version.
  Use this one by default.
- **`10.5281/zenodo.22113220`** — the **version** DOI for **v1.0.0**. Use this when you need to
  pin DeepGOF-1 exactly, for example when reproducing a number from the paper.

> Ebrahim, E. K. (2026). *DeepGOF-1: pretrained goodness-of-fit test for logistic regression —
> weights, training corpus, benchmark harness and per-replicate results* (v1.0.0) [Software].
> Zenodo. https://doi.org/10.5281/zenodo.22113219
