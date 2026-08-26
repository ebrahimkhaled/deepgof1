# 17 — THE PROBIT SHOWCASE: PRE-DECLARATION (frozen 2026-08-21, before any data was generated)

**Frozen.** Reported against verbatim. This gate exists because the author proposed a dedicated
link-misspecification showcase, and because **a probit claim was RETRACTED from this project on
2026-08-17** (`03_C1_DEFECTS.md` D4). Nothing from that retraction is reused here. This is a fresh,
pre-declared experiment whose entire purpose is to answer honestly a question we previously
answered carelessly.

---

## 1. Why this case, scientifically

Probit-versus-logit is the acknowledged **hardest** link alternative: both links are symmetric and
differ only in tail weight, and the applied literature routinely switches to cloglog precisely
because probit is too hard to distinguish. Low power against link misspecification is a documented,
general weakness of the classical battery, with Stukel's test named as the best of a weak field.
If a pretrained test detects what eleven classical tests individually cannot, that is the single
most persuasive demonstration available — **provided it is measured properly this time.**

## 2. What went wrong before (so it cannot recur)

| failure in the retracted claim | the fix here |
|---|---|
| per-family claim taken from a run whose pre-declaration forbade per-family claims | this gate is *about* one family and says so in advance |
| corpus-seed variance ±.1 on single families, never measured | **three independent training corpus seeds**; the seed spread is a REPORTED quantity, and the headline is the across-seed MEAN with its range |
| oracle calibration | **deployed** per-dataset bootstrap for the learned test; analytic references for the classical tests, exactly as a user would run them |
| ceiling comparison against a single "best individual test" chosen post hoc | all 11 tests reported individually; the comparator is fixed in advance as `max` over the 11 **and** the Cauchy combination |

## 3. The design (fixed)

- **Alternative:** probit link, fitted as logistic `y ~ x1 + x2`, at **two amplitudes**:
  - `mild`  — `y ~ Bernoulli(Phi(0.55*x1 - 0.35*x2))`, the existing generator in
    `learned-def/R/scenarios.R:36`, unchanged.
  - `cloglog` — `y ~ Bernoulli(1 - exp(-exp(0.3 + 0.8*x1 - 0.5*x2)))`, an **asymmetric** link.
    The contrast cell: link misspecification that the classical battery *can* sometimes detect.
  - `probit-extreme` — `y ~ Bernoulli(Phi(1.10*x1 - 0.70*x2))`, doubled slopes. **Reported as a
    CONSERVATISM cell, not a power cell** (see the amendment note).

  **Amendment note (2026-08-21, after a 140-dataset timing pilot on throwaway seed 999000, before
  any confirmatory data).** The pilot was run to measure wall-clock, and it falsified my second
  amplitude. Doubling the slopes does not amplify the link difference — it drives η to ±3.5, so
  fitted probabilities sit near .0002/.9998 and the data is nearly deterministic. The classical CCT
  then rejected **0.000** with a mean p-value of **.758**: the battery is not blind there, it is
  massively *conservative*, which is a different phenomenon and a well-known property of these tests
  under extreme fitted probabilities. Keeping that cell as a "power" comparison would have
  manufactured a flattering result from an artifact. It is retained, relabelled as a conservatism
  measurement, and `cloglog` replaces it as the second genuine link cell. The pilot block is
  deleted and none of its numbers enter any claim.

  **Why more than one link (phase 17.1's lesson).** A comparison run on the **power floor** is
  uninformative: when every arm rejects at ≈ α, no arm can separate from another. Probit is the
  hardest alternative in the battery, so at n = 50 it may sit on that floor; `cloglog` guarantees at
  least one link cell is off the floor at every n. **Criteria apply to each cell separately and
  every cell is reported**, floor or not.

  **Stated expectation, recorded before the result (so it cannot be revised afterwards):** a link
  departure is smooth and monotone along the linear predictor, which is precisely what
  `DEF.stukel` — a test purpose-built for link misspecification — is designed to catch, whereas the
  residual map is a *location* detector over (x1,x2). **PB1 may well fail on the probit cell.** If
  it does, that is the reported result.
- **Sample sizes:** n ∈ {50, 100, 200} (the paper's scope).
- **Covariate design:** identical for null and both alternatives — `x1 ~ U(-3,3)`, `x2 ~ N(0,1)` —
  so the link is the ONLY thing that differs between H0 and H1.
- **Replicates:** 2000 null + 2000 probit-alternative datasets per n, per corpus seed.
- **Seeds:** three training corpora (SEED0 = 101000, 102000, 103000), each producing its own
  deployment-matched rule. **Everything is reported per seed and as the mean.**
- **Learned rule:** the shipped construction — rank-grid residual map, variance-standardised
  cells, small CNN, **deployment-matched training** (null class = bootstrap replicates), B = 199
  exact rank p-value. Trained on the **general random-departure library, with NO probit family in
  training** — this is a held-out-kind test, and that is the point.
- **The mix parameter is pinned to `mix = 1.0`** (pure deployment-matched: the training null class is
  entirely bootstrap replicates). This is the shipping construction of `13_SMALLN_PAPER_SPECKIT.md`
  §2.3 and the phase-15 winner. `mix = 0.5` is computed and deposited **for information only** — no
  criterion may be evaluated on it, and no claim may cite it. Pinning it here removes the
  after-the-fact choice between two arms that was available at phase 15.
- **Classical arms:** all 11 tests at their analytic references, plus the Cauchy combination.
- **Matched-level comparison is RANDOMIZED** so every arm sits at exactly α = .05 despite the
  deployed p-value living on a 1/200 grid: with `t0` the largest attainable threshold whose null
  rejection rate `F(t0) ≤ α` and `t1` the next one up, reject with probability
  `γ = (α − F(t0))/(F(t1) − F(t0))` in between. Applied identically to all thirteen arms. Both the
  raw nominal-.05 numbers and the randomized-matched numbers are reported.
- **Comparison:** at each rule's own realised level (matched-level power), and additionally at raw
  nominal .05 with realised sizes printed beside every number.

## 4. Criteria (declared now)

| # | Criterion | Threshold |
|---|---|---|
| **PB1** | Learned matched-level power exceeds the **best of the 11 individual tests** at every n | mean over seeds, and in ≥ 2 of 3 seeds individually |
| **PB2** | Learned matched-level power exceeds the **Cauchy combination** at every n | same rule |
| **PB3** | Across-seed **range** of the learned probit power | ≤ **.10** (the retracted claim's variance was ~.09 across corpora; if it is this large again, no per-seed number may be quoted — only the mean with its range) |
| **PB4** | Learned deployed size at each n | ∈ [.035, .065] |

## 5. Fallbacks (declared now)

- **All four pass** → the paper gets its showcase figure and one sentence: *on the hardest known
  link alternative, at the sample sizes where it matters, the pretrained test detects what none of
  eleven classical tests can.* Stated with the seed range attached, always.
- **PB1/PB2 pass, PB3 fails** → the effect is real but unstable; report the mean and range, and
  demote the figure from headline to supplement. **No single-seed number is ever quoted.**
- **PB1 or PB2 fails at any n** → report the failure in the manuscript's honest-boundary section.
  The retraction stands and is cited as such. **This is an acceptable outcome and will not be
  re-run with new seeds.**
- **PB4 fails** → the learned arm is not validly calibrated in this cell; report and investigate,
  claim nothing.

## 6. Forbidden

Re-using any number from the retracted claim (.147 / .080 / .170). Quoting a single-seed figure if
PB3 fails. Choosing the "best individual test" after seeing results — it is `max` over the 11,
fixed here. Adding seeds until the answer improves.

## 7. Execution

`code/probit_showcase.R` (corpora + deployed bootstrap, per seed, per n) →
`code/probit_showcase_score.py` (deployment-matched training + scoring + verdict) →
`code/figs/fig5_probit.py` (the figure, built only if PB1–PB2 pass).
Deposit: `results/probit/` — every per-replicate p-value for every arm, every seed.

---

# RESULT (2026-08-21, after all 9 blocks — criteria above were frozen before any data)

**PB1 FAIL · PB2 FAIL · PB3 FAIL on cloglog (PASS on probit) · PB4 PASS.**

| cell | n | learned (mean, 3 seeds) | range | best of the 11 | gap |
|---|---|---|---|---|---|
| probit | 50 | .0415 | .0081 | .0537 (Pigeon-Heyse) | −.0122 |
| probit | 100 | .0512 | .0182 | .0620 (Tsiatis) | −.0108 |
| probit | 200 | .0460 | .0240 | .0507 (Tsiatis) | −.0047 |
| cloglog | 50 | .0244 | .0142 | .0717 (HL-equalwidth) | −.0473 |
| cloglog | 100 | .0844 | .0376 | .1125 (DEF.poly2) | −.0281 |
| cloglog | 200 | .1064 | **.1489** | **.2438** (DEF.poly2) | **−.1375** |

Deployed size .0467 / .0450 / .0533 — PB4 passes at every n. The test is valid here; it is not more
powerful. The winning classical arm is a **directional** test every time (`DEF.poly2`, `Tsiatis`),
which is the mechanism the pre-declaration named in advance.

**The fallback that fires** (§5, third bullet): *"PB1 or PB2 fails at any n → report the failure in
the manuscript's honest-boundary section. The retraction stands and is cited as such. This is an
acceptable outcome and will not be re-run with new seeds."* — plus §5 second bullet, since PB3 also
failed on cloglog: **no single-seed number from this gate may be quoted anywhere**, mean and range
only, and Figure 5 goes to the supplement rather than the body.

**What is now FORBIDDEN by this result**, in addition to §6: claiming any probit or link-family
advantage; presenting Figure 5 as a showcase; using the word "showcase" for this material; adding a
fourth seed, a different amplitude, or an η̂-indexed map geometry *and then reporting the outcome as
if it belonged to this gate*. An η̂-indexed map is a legitimate scientific question but it is a NEW
pre-declared gate with its own frozen criteria — never a rescue of this one.

**What this gate DID establish, and what the manuscript should carry:**
1. A characterized boundary: the residual map detects departures with a **location** in covariate
   space; directional tests own departures **monotone along the linear predictor**. Complementary,
   not competing. Backed by a −.1375 measurement at n=200, ~14 SE.
2. The corpus instability that caused the 2026-08-17 retraction, measured: 7.6× swing (.0226 /
   .1249 / .1716) on cloglog at n=200 while the analytic rival moved .005. Every learned-power
   number in this paper must carry an across-seed range.
3. The Cauchy combination dilutes a lone strong component (.1037 vs .2438); adding the learned arm
   does not rescue it (.1152).
4. score-AUC declines monotonically with n in all three corpora — the bootstrap-vs-truth gap is a
   small-sample artifact of estimating β̂.
5. The battery's small-n fragility (S1): Pigeon-Heyse NA on 2.50% of n=50 datasets, separation flag
   1.14%, both ~0 by n=200; R's `mean()` propagates NA into the Cauchy combination.
