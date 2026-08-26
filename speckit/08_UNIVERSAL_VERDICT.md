# 08 — UNIVERSAL TEST VERDICT (2026-08-19)

Corpus: 250k null + 250k alternative training datasets spanning n ∈ {100…1400}, p ∈ {2,3,5,8},
two covariate laws, a 3-level factor, two event rates. Judged **only** on ten held-out classes
containing sample sizes, dimensions, covariate laws, category counts and event rates that never
appeared in training. Criteria frozen in `07_UNIVERSAL_PREDECLARATION.md` before any data existed.
Raw output: `results/universal/UNIVERSAL_VERDICT.txt`; every score deposited.

---

## 1. The verdict against the frozen criteria

| | P1 size | P2 uniform | P3 invariance | P4 power |
|---|---|---|---|---|
| **Route A** (augmentation only, one fixed threshold) | FAIL | PASS | FAIL (.0348) | FAIL (−.037) |
| **Route B** (amortized calibration) | FAIL | FAIL | FAIL (.0240) | **PASS (−.0012)** |

**Neither route met the pre-declared bar.** Reported as it fell; no criterion adjusted, no re-run.

**Route A answers the author's question directly: breadth alone is not enough.** Training across
seven sample sizes and four dimensions left a single fixed threshold running conservative almost
everywhere (.027–.048) and costing .037 of power. Augmentation without a mechanism that *rewards*
invariance does not produce invariance.

---

## 2. ⚠️ THE RESULT THAT OVERTURNS THE PREVIOUS VERDICT

`05_GBOOT_VERDICT.md` concluded: *"the learned score is not pivotal; the classical Cauchy
combination is."* **That conclusion was wrong, and this run refutes it.** It rested on two designs
(2-covariate and 4-covariate, n = 200, 50% events) — too narrow to detect the classical rule's own
miscalibration.

Across ten genuinely different designs, using its **analytic Cauchy null** (the way it is actually
used in practice — no calibration at all):

| | Route B | **CCT-full (classical)** |
|---|---|---|
| spread of realized size (max − min) | **.0285** | .0482 |
| worst case | **.0715** | **.0977** (≈2× nominal) |
| classes outside [.035, .065] | **2 of 10** | **4 of 10** |
| mean \|size − .05\| | **.0089** | .0158 |
| KS vs U(0,1), pooled | p ≈ 0 | **p = 1.2e-40** |

Per class, the classical rule reaches **.0810 at p = 12** and **.0977 in the 5%-event regime**. It
is not pivotal either. **Neither method holds its level across a wide design space, and the learned
rule with amortized calibration is the better-calibrated of the two — by roughly a factor of two on
every summary — while matching its power (−.0012).**

## 3. Power

Mean over the non-extrapolation classes: Route B −.0012 versus CCT-full, i.e. **a tie**, which the
pre-declaration counts as a pass because the deliverable is that no per-dataset calibration is
needed. Route B wins where the design is unusual (H5 correlated .485 vs .432; H7 categorical .405
vs .365) and loses on the skewed covariate (H6 .323 vs .356) and the sparse regime (H8 .132 vs
.184). Route A is uniformly weaker (−.037).

## 4. A methodological lesson worth reporting

**Route A passed the pooled uniformity test (P2, p = .024) while being badly miscalibrated within
every class.** Pooling averages conservative classes against liberal ones and hides both. P3 —
the spread of realized size *across* design classes — is the criterion that actually detects the
failure. Any paper checking a learned test's calibration on pooled null p-values can pass that
check while shipping a rule that is wrong everywhere. Report per-design size, always.

## 5. What the paper is now

The thesis is stronger than either previous version, and it is no longer a claim about our method
being better:

> **Goodness-of-fit tests for logistic regression — classical and learned alike — do not hold their
> nominal level across the range of designs analysts actually encounter.** A pretrained score with
> amortized calibration is the better-behaved of the two, but neither is exact. The commonly
> reported evidence hides this twice over: by evaluating against a null reference drawn from the
> true model, and by pooling null p-values across designs.

That is a claim about the field's standard practice, supported by 544,000 datasets, ten held-out
design classes, and a deposited per-replicate record. It is more useful, and much harder to
dismiss, than "our test is more powerful."

## 6. Honest status of every claim now on record

| Claim | Status |
|---|---|
| Learned test beats classical on power | ⛔ dead (G-BOOT: −.015 deployed) |
| Learned score is not pivotal | ✔ holds |
| **Classical CCT *is* pivotal** | ⛔ **RETRACTED here — .0977 at 5% events, .0810 at p = 12** |
| Amortized calibration removes the bootstrap | ✔ works, and beats the analytic classical route |
| One fixed threshold works everywhere | ⛔ dead (Route A) |
| No GOF rule tested here holds level across designs | ✔ **the new headline** |

## 7. Next

1. **Route C (adversarial pivoting)** — now clearly worth running: it is the only route that
   directly optimises the quantity (P3) that everything failed on. Needs the CUDA torch build.
2. **Widen the calibration descriptors.** Route B's two failures (p = 12, ρ = .8) are both
   design-complexity directions; the descriptor set may simply lack the information. Cheap to test
   — the corpus is already on disk, only the calibration network is retrained.
3. **Sparse regime (H8) deserves its own study** — every rule degrades there, and it is the regime
   this research family cares most about.
4. Re-check whether the classical rule's miscalibration is already known in the literature. If it
   is not, that alone is a publishable note.

---

## 9. Appendix (2026-08-19, same day): "calibrate-then-combine" was tried and is WORSE

**The hypothesis.** The 11 classical inputs are themselves miscalibrated (A.poly .052-.084 liberal
everywhere; A.spline .046-.076; Xie .020-.058 conservative), and a Cauchy combination is dominated
by the smallest p-value, so it *amplifies* whichever component is most liberal - which is why the
combination's spread (.048) exceeds any single component's (worst .037). Dropping the offenders is
not an option: it costs **.161 of power** (.444 -> .283), because the miscalibrated tests are
precisely the powerful ones. So the proposal was to fix each component's level conditionally on
design (11 easy low-dimensional problems), then combine p-values that are already uniform - at
which point the Cauchy null is analytic and pivotal by construction.

**The measurement refutes it.** `code/calibrate_then_combine.py`, same held-out classes:

| | calibrate-then-combine | classical CCT | Route B (combine-then-calibrate) |
|---|---|---|---|
| spread of size | **.1125** | .0482 | **.0285** |
| worst case | .1640 | .0977 | .0715 |
| out of [.035,.065] | 8/10 | 4/10 | 2/10 |
| mean power | .4487 | .4436 | (tie with CCT) |

Power was preserved (indeed marginally best), but **calibration collapsed**.

**Why - two mechanisms, separated by a uniformity diagnostic.**

1. *Near the training grid* (H1a-c, H3, H6) the recalibrated combination has median ~.50 and
   P(p<=.5) ~ .50 - the **body** of the conditional distribution is fitted correctly - while
   P(p<=.05) runs .052-.072. Only the **tail** is wrong. The probability integral transform needs
   the *whole* conditional CDF accurate; the pinball loss spread over 30 quantile levels fits the
   body at the tail's expense. Route B needs only **one** conditional quantile (tau = .95) and can
   weight the loss there - a far easier estimation problem.
2. *Outside the training grid* (H2 n=2000, H4 p=12, H7 five-level factor) the median falls to
   .381 / .305 / .379 - the whole distribution is displaced. The quantile model does not
   extrapolate, and these are exactly the three axes whose held-out values exceed the trained range.

**The structural lesson, which is the transferable part.** Per-component calibration multiplies the
number of models that must extrapolate correctly - eleven instead of one - and the Cauchy rule then
**amplifies whichever one extrapolated worst**. Concentrating all the estimation into a single
calibration model is both an easier statistical problem and immune to that amplification.
**Combine-then-calibrate is the better architecture, for a reason that was not obvious in advance
and is worth stating in the paper.**

⚠️ **Process note.** The first run of this experiment reported power .028 against a size of .05-.075.
Power below size is arithmetically impossible against genuine alternatives, which is the only reason
a tail inversion in `to_uniform` (returning F(z|d) where a p-value needs 1-F(z|d)) was caught in
minutes. A subtler sign error - one producing merely disappointing numbers rather than impossible
ones - would have read as a real negative finding. **The evaluator should assert
power >= size on every alternative cell**, rather than relying on someone noticing.

---

## 10. Appendix (2026-08-19): THE AUTHOR'S RESIDUAL IDEA WORKS -- the features themselves were the leak

**The question:** "can't we use the residuals better?" It turned out to name the root cause.

The 33 features carry grouped **Pearson** residuals, (O-E)/sqrt(V), which are only ASYMPTOTICALLY
N(0,1). The approximation fails exactly where the test fails. Measured directly on 12 design cells
x 2000 nulls (`code/rqr_pivotality.R`, no battery, no learning -- the FEATURE alone):

| design | sd of grouped Pearson residual | sd of grouped RQR |
|---|---|---|
| n=400, 50% events | 0.920 | 0.947 |
| n=400, 5% events | 0.900 | 0.954 |
| **n=100, 5% events** | **0.763** | **0.951** |
| n=2000 | ~0.92 | ~0.95 |

**Pearson residuals lose a quarter of their scale in the sparse small-n cell; Dunn & Smyth (1996)
randomized quantile residuals hold 0.95 everywhere.** For binary y: u ~ U(0,1-phat) if y=0,
u ~ U(1-phat,1) if y=1, r = qnorm(u) -- EXACTLY N(0,1) under correct specification at any n and any
event rate, not asymptotically.

**Consequence for pivotality of the derived statistic (sum of squares over 10 groups):**

| | Pearson | **RQR** |
|---|---|---|
| spread of the .95 null quantile across 12 designs | 2.81 | **0.84** |
| relative spread (max/min) | 1.211 | **1.051** |
| **realized-size spread under ONE fixed cutoff** | **.0189** | **.0090** |

**Halving the size spread by changing one feature definition exceeds every gain any architecture
change produced.** The non-pivotality was entering through the inputs, and no downstream calibration
layer can fully undo an input whose scale moves with the design.

**Honest caveats.** (a) RQR is still not perfectly pivotal (.0090 > 0), because beta-hat estimation
remains; part of the residual .0090 is Monte-Carlo noise at 2000 reps per cell. (b) RQR introduces
**randomization** -- the p-value becomes random given the data. That is a genuine cost and must be
stated: either report the randomized test honestly, or average several draws and lose exactness.
Dunn & Smyth's own framing accepts the randomization as the price of exact normality.

**Status: EXPLORATORY.** The next step is to rebuild the feature vector with RQR-based groupings and
re-run Routes A/B on a fresh corpus under a pre-declaration, exactly as C1 was handled. Nothing here
enters a manuscript until that happens.

**Prior-art note:** randomized quantile residuals are standard GLM diagnostics (Dunn & Smyth 1996,
JCGS). Using them as *features for a learned goodness-of-fit test* is, as far as the C3 sweep and
this one show, not established -- but it must be checked properly before any novelty claim.
