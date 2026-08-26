# 10 — THE RESIDUAL-MAP IDEA (author's, 2026-08-20): extract power from the clean channel

## The measurement that makes this the right move

`decompose_features.py` split the rule into its two input channels:

| channel | power | size-spread across held-out designs |
|---|---|---|
| 11 classical z-scores | **.4149** | **.0365** (noise floor .0139) |
| 22 Pearson residual codings | .1575 | .0190 |
| 22 RQR residual codings | .0824 | **.0160** |
| both (the full rule) | .4066 | .0380 |

**All the power sits in the miscalibrated channel; all the good calibration sits in the powerless
one.** Every attempt so far has tried to repair the broken channel (recalibrate the components,
amortize, adversarially pivot). The author's proposal inverts it: **extract power from the channel
that is already well calibrated**, and the level problem disappears rather than being corrected.

The reason the residual channel is powerless is not that residuals lack information — it is that we
summarise them terribly: sort by fitted risk, average into ten bins, take order statistics. Ten
numbers from a single one-dimensional projection.

## The author's framing, and what it corrects

> "like a photo we make a feature map to detect the object — for the raw data we have, what to
> detect the GOF"

The correction it supplies: a convolutional feature map finds an object **anywhere in the frame**.
Our residual "image" has always been **one-dimensional** — ordered by fitted risk. But misfit lives
in **covariate space**. A local bump in a corner of the (x1, x2) plane is smeared across every risk
bin. Phase 1.3 measured exactly this (raw grouped residuals scored −.146 because the bump lands in a
random risk bin each replicate) and we worked around it by *sorting* the bins — which discards
location rather than handling it.

**The representation should be the residual surface over covariate space**, where locality is real
and convolution is the right tool.

## The second half: the tests as teacher, not as input

> "give him the opinion of some tests so he learns the opinion of tests in those patterns"

Correct instinct, with one trap. Feeding the z-scores as **inputs** is what currently imports their
miscalibration. Instead use them as **auxiliary training targets**: a second head predicts what each
classical test would report for this pattern. The representation is then forced to encode what the
battery knows, while **inference uses only the residual map**, so the calibration stays clean.
Distillation, not concatenation.

## Why this is not a rerun of C2

C2 (manifest phase 3) put a CNN on 50-bin residual **curves** and lost twice (.7546 / .7479 vs
.7955). The objection is fair and must be answered in the paper. What differs:

| | C2 | this |
|---|---|---|
| geometry | 1-D, ordered by fitted risk | **2-D over covariate space** — location is meaningful |
| scale | 12k training rows | 400k+ available |
| supervision | label only | **+ test distillation** |
| bins | 50 bins × ~20 obs = noise | coarse grid, many obs per cell |

C2's lesson was "the right resolution beats the deeper reader" — a *resolution* finding along one
projection. This changes the projection itself.

## Open design problem: variable p

A 2-D image is natural at p = 2. For general p the candidates are (a) 2-D grids over the leading
principal components of X, (b) p channels of (x_j, fitted risk) grids, (c) a permutation-invariant
set network over the point cloud {(x_i, r_i)}. **Decide only after the p = 2 proof of concept** —
if a residual map cannot beat the 10-group summary where the geometry is easiest, the idea is dead
and no amount of generalisation machinery rescues it.

## GATE M1 — the decisive first experiment (declared before running)

**Question.** At p = 2, where a residual map is unambiguous, does a 2-D residual map extract
materially more power than the 10-group summary — from residuals **alone**, with no z-scores?

**Fixed:** n = 400, p = 2, the standard departure library, 60k null + 60k alternative training rows,
2000 null + 2400 alternative held-out. Map = 12×12 grid over the **ranks** of (x1, x2) (rank grid
makes the input scale-free by construction), cell value = sum of residuals / sqrt(cell count).

**Arms:** (1) 10-group summary (the incumbent, residuals only); (2) CNN on the 12×12 map;
(3) CNN + test-distillation auxiliary head. All residuals-only at inference.

**Criterion M1.** Arm 2 or 3 reaches power ≥ **.25** while keeping size-spread ≤ **.025**.
Rationale: the residual channel currently gives .16 (Pearson) / .08 (RQR); the z-score channel gives
.41. Reaching .25 from the clean channel alone would be a genuine step and would justify the larger
build. Below that, the idea does not pay for its complexity.

**Fallback if M1 fails.** The residual channel is information-poor, not merely badly summarised;
revert to attacking the three miscalibrated battery tests (A.poly, A.spline, Xie) directly, and
record M1 as a measured negative — which also retires the "we never tried a proper feature map"
objection permanently.

**Cost.** ~120k datasets at n = 400 with the fast battery ≈ 25 min; CNN training on GPU (now
available, torch 2.6.0+cu124, GTX 1660 Ti verified) ≈ minutes.
