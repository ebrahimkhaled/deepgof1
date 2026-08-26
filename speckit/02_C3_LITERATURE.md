# 02 — C3 LITERATURE SWEEP: what is new here, stated so a referee can check it

**Run 2026-08-17.** Purpose: find the paper that scoops us, before a referee does. It does not
exist — but four literatures come close, and the manuscript must position against all four
*explicitly*, because a JASA/Biometrika referee will know at least two of them.

Verification column: **F** = abstract fetched and read; **S** = search-result summary only
(must be fetched before the manuscript cites it).

---

## 1. The four neighbouring literatures

### A. Neural networks as the ALTERNATIVE model (sieve GOF)

| Work | What it does | Why it is not us |
|---|---|---|
| Shen, Jiang, Kong & Yang (2021), *Statistics & Probability Letters* — "A goodness-of-fit test based on neural network sieve estimators" (S) | Fits a neural network as a nonparametric sieve alternative; the test statistic compares the parametric fit to that sieve; asymptotic null distribution derived | The net is **fitted to the user's data** as an alternative model. Nothing is learned offline; nothing transfers between datasets. Ours never fits the user's data — it *reads* a fixed 33-number summary of it |
| ReLU-network GOF statistics for genetic association testing (2025, PMC) (S) | Deep ReLU nets to test input-feature significance; more powerful than the t-test under nonlinear signal | Feature-significance testing, not model-adequacy testing; again fitted per dataset |

**Position in the paper:** these are the "neural network as a flexible alternative" school. State
in one sentence that our net is not an alternative model at all — it is a **learned decision rule
over classical diagnostics**, which is why its cost at use time is one forward pass.

### B. Learn a distinguisher on the user's own data (the adaptive/split school) — the closest rivals

| Work | What it does | Why it is not us |
|---|---|---|
| **GRASP** — Javanmard & Mehrabi (2024), **JRSS-B 86(1) 215** (S; arXiv 2209.02064) | Randomization test for a black-box binary classifier; f-divergence tolerance hypothesis; **finite-sample** type-I control, distribution-free | Per-dataset; targets a black-box classifier's conditional law, not the adequacy of a *fitted parametric regression*; no offline training, no transfer |
| **BAGofT** — Zhang, Ding & Yang (2021) (S; CRAN `BAGofT`) | Learns an **adaptive partition** on a held-out split, then a chi-square style statistic on the other split | The learning happens **inside each application** and costs sample. Ours spends zero of the user's sample on learning. ⚠️ Already a comparator in the family — the `nsplits=1` handicap trap is recorded in [[reference_simulation_registry]] |
| **Chen & Lei (2026)**, arXiv:2604.05470, "Evaluating Black-Box Classifiers via Stable Adaptive Two-Sample Inference" (**F**, 7 Apr 2026) | Draws a synthetic sample Y′ ~ Multinom(η̂(X)), trains a **distinguisher** to separate real from synthetic, rank-sum statistic, sample splitting + cross-validation CLT under stability | The newest and most dangerous neighbour rhetorically — it is *also* "train a classifier to detect misfit". But it trains **on the user's data with sample splitting**; validity comes from a stability condition + CV-CLT. Ours trains **once, offline, on simulation**, and validity comes from Monte Carlo calibration, which needs no stability condition at all |

**This is the row that decides the paper's framing.** The distinction is not "we use a network and
they don't" — it is **where the training sample comes from**. Theirs is the user's; ours is the
simulator's. Everything the paper claims (no sample cost, fixed shipped weights, zero marginal
training cost, validity independent of the learner) follows from that one difference. Say it in
the abstract.

### C. Learned test statistics with Monte Carlo calibration (physics)

| Work | What it does | Why it is not us |
|---|---|---|
| **NPLM** — D'Agnolo & Wulzer (2019), "Learning New Physics from a Machine" (arXiv 1806.02350); Grosso et al. (2025), arXiv 2508.02275 (S) | Classifier estimates the log density ratio to a reference sample; statistic = 2·max log-LR; **calibrated by Monte Carlo toys**; explicitly framed as a general-purpose GOF methodology | The calibration philosophy is **identical to ours** (and to LEGofT's) — cite it as the precedent, do not pretend MC calibration is new. But NPLM **trains on the observed dataset in-sample** against a reference sample, per analysis. It also presumes the reference can be sampled at will, which in our setting is the fitted model itself |

**Position:** NPLM is the honest precedent for "a network can be a valid test statistic if you
calibrate it by simulation". Our contribution on top of it is **amortization** — the training is
moved out of the analysis entirely and reused across datasets.

### D. Amortized / simulation-based inference (the framing literature)

| Work | Relevance |
|---|---|
| Cranmer, Brehmer & Louppe (2020), *PNAS* — the SBI frontier (S) | The canonical citation for "train on simulation, apply to data" |
| Zammit-Mangion et al., "Neural Methods for Amortized Inference", *Annual Review of Statistics* (S) | Names the two-stage structure we use: an offline training stage, then a cheap inference stage per dataset |
| Radev et al. (2023), "A Deep Learning Method for Comparing Bayesian Hierarchical Models", arXiv 2301.11873 (S) | **The nearest amortized-testing neighbour**: amortizes *model comparison* over a fixed candidate set of models |

**The gap, stated precisely.** Amortization is established for posteriors, likelihoods, likelihood
ratios and **comparison among an enumerated set of candidate models**. It has not been used for
**omnibus goodness-of-fit of a fitted parametric regression**, where the alternative is not
enumerable — which is exactly why the training prior has to be a *randomized library of
departures* rather than a list of rival models. That randomized library is the paper's technical
core, and it is what makes the kind-holdout result (a shape absent from training still detected)
meaningful rather than trivial.

---

## 2. The novelty claim, in the form it will be defended

> **Existing learned GOF tests spend the analyst's sample to learn. This one spent a simulator's.**

Three claims follow, each falsifiable:

1. **Amortized GOF for a fitted parametric regression is new.** No prior test ships fixed weights
   trained offline on a randomized misspecification library. (Nearest: Radev 2023 amortizes model
   *comparison*; NPLM learns per-analysis; GRASP/BAGofT/Chen–Lei learn per-dataset from the user's
   own sample.)
2. **Zero sample cost.** The user's n is spent entirely on the test, not split between learning
   and testing. This is a *measurable* advantage over BAGofT/GRASP/Chen–Lei, not a rhetorical one
   — and it is the mechanism behind the n=200 result (+.118, the largest gain of the three sizes,
   precisely where a split hurts most).
3. **Validity does not rest on the learner.** Monte Carlo calibration (Barnard 1963) makes size
   exact by construction regardless of what the net learned — verified on untouched test-null
   cells at .059/.040/.065. No stability condition, no CV-CLT, no asymptotics in the learner.

**What we must NOT claim:** not that learned GOF testing is new (it is not — B and C); not that MC
calibration is new (Barnard 1963, and NPLM); not that neural GOF is new (A). The novelty is
**amortization + the randomized departure library**, and nothing wider.

---

## 3. The probit claim — VERIFIED SAFE, and stronger than we thought

The literature independently documents both halves of what our probit result needs:

- **"Low power of most goodness-of-fit tests to detect link misspecification"** is a reported,
  known limitation in the GLM literature (S).
- **Probit-vs-logit is specifically the hardest case** because both links are symmetric; authors
  routinely switch to cloglog *because* probit is too hard to distinguish. Stukel's test is named
  as the most powerful available for link departures (S).

So the C1 probit numbers — ensemble .095/.114/.170 vs CCT-full .032/.064/.059 across
n=200/500/1000, with the equal-weight rule never leaving the noise floor — land on an alternative
the field has already conceded is near-undetectable. **This is the paper's single most citable
sentence, and it is defensible.** ⚠️ Both statements above are from search summaries; before the
manuscript cites them, fetch and quote the primary sources (Stukel 1988 and one GLM link-power
study) so the sentence rests on a citation, not on a search snippet.

---

## 4. What C3 changes about the paper

1. **The abstract's first line changes** from "a pretrained test" to the *where-does-the-training-
   sample-come-from* contrast — that is the sentence a JRSS-B referee who knows GRASP needs.
2. **A comparator is now mandatory, not optional.** Reviewers will ask for BAGofT and, after
   April 2026, plausibly Chen–Lei. BAGofT is already wired into the family's battery (with its
   handicap trap documented). Chen–Lei has no released software found in this sweep — check again
   before submission; if none exists, say so plainly rather than implementing a rival by hand.
3. **A related-work section of four paragraphs**, one per literature A–D, each ending in the one
   sentence that separates us. Do not bury these in a citation list.
4. **Venue reading unchanged** (JASA T&M → Biometrika → JCGS). But note GRASP landed in **JRSS-B**,
   which is evidence the top methodology venues are actively buying work in exactly this corner.

## 5. Open verification debts (before any submission)

- [ ] Fetch and read in full: GRASP (JRSS-B 2024), BAGofT (2021), Shen et al. (2021), Radev et al.
      (2023), Cranmer et al. (2020), D'Agnolo–Wulzer (2019). Six papers marked **S** above.
- [ ] Re-run the Chen–Lei check nearer submission — it is four months old and may have shipped code.
- [ ] Search once more under GOF-specific vocabulary the SBI community does not use
      ("calibration test", "model adequacy", "specification test") — this sweep used the ML framing.
- [ ] Confirm no `ebrahim.gof`-family paper already claims amortization loosely, which would create
      a self-priority problem.
