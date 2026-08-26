# 06 — THE UNIVERSAL TEST: spec-kit for a pretrained GOF test that needs no per-dataset calibration

**Written 2026-08-19**, after G-BOOT showed the learned score is not pivotal
(`05_GBOOT_VERDICT.md`). This document specifies the research programme that decides whether the
original vision is reachable — *download the model, feed it any dataset, get a p-value* — or
whether the deployment-gap warning is the paper.

Criteria and fallbacks are declared here, **before any run**, as in every gate of this project.

---

## §1. The goal, stated as the author stated it

> Like learning "Ebrahim" from many photos — any photo comes in, it knows, whatever the resolution.

Formally: a score S(data) whose **null distribution does not depend on the nuisance
characteristics of the dataset** — sample size, number of covariates, covariate law, event rate.
If that holds, one threshold shipped with the software is valid everywhere and the bootstrap is
never needed.

That property has a name in statistics — **pivotality** — and G-BOOT measured that our current
score does not have it: null rejection .0250 on 2-covariate designs versus .0583 on 4-covariate
designs, and a KS test against uniform at p = .0001. The classical Cauchy combination *is* pivotal
(.0500 / .0533, KS p = .236), which is exactly why it needs no bootstrap.

**The author's own analogy contains the diagnosis.** A face model trained only on sharp photos does
degrade on blurry ones; the standard fix is augmentation — train across resolutions. Our corpus was
trained at **one sample size (n = 200)**, **two designs**, **one covariate law**. We never asked it
to be invariant, so there is no reason it would be.

## §2. Prior art — what already exists, and what we may claim

Two literatures solve exactly this, and both must be cited as precedent rather than reinvented.

| Work | What it gives us | What is still ours |
|---|---|---|
| **Louppe, Kagan & Cranmer (2017), NeurIPS — "Learning to Pivot with Adversarial Networks"** ([arXiv:1611.01046](https://arxiv.org/abs/1611.01046)) | The training procedure for **forcing** a network's output to be independent of a nuisance parameter: an adversary tries to predict the nuisance from the score; the score network is trained to defeat it. Tends to a minimax solution that is both optimal and nuisance-independent *when one exists*, with a tunable power/robustness trade-off. | Applying it to make a **GOF statistic** pivotal in *n* and design. The nuisance here is the dataset's own shape, not a physical systematic. |
| **LF2I / ALFFI — Dalmasso, Masserano, Zhao, Izbicki & Lee** ([arXiv:2107.03920](https://arxiv.org/pdf/2107.03920), [arXiv:2306.07769](https://arxiv.org/abs/2306.07769)) | The alternative: don't force invariance — **learn the critical value**. ALFFI approximates the null CDF with a neural network as a function of *both* the parameter point and the observed statistic, amortizing calibration over all future uses. Provably valid confidence sets given a good simulator. | Their target is parameter inference with a simulator. Ours is **model adequacy of a fitted regression**, where the "parameter" is the nuisance design and the simulator is the fitted model itself. |

⚠️ **This narrows the novelty claim again, and it must be narrowed in writing before drafting.**
We may **not** claim amortized calibration, adversarial pivoting, or learned critical values as new.
What remains defensible: *applying amortized calibration / adversarial pivoting to omnibus
goodness-of-fit of a fitted parametric regression, and measuring whether the resulting rule
transfers to designs never seen in training.* If even that is taken, the honest contribution is
the **measurement** — see §8.

Checked and NOT a conflict: Guchenko (arXiv:2607.15015, July 2026) augments GOF tests with
secondary statistics under per-dataset Monte Carlo calibration — no learning, no amortization.

## §3. The three routes — run all three, they share one corpus

The score network and the calibration network are kept **strictly separate**, and this separation
is the design's main safeguard:

- **Score network sees only misfit features.** It is *design-blind by construction*. This is not a
  stylistic choice — the leakage autopsy (`MANIFEST.md` phase 1.2) proved that giving these models
  design descriptors makes them read the design fingerprint instead of the misfit (`lev_sd`
  |coef| 462; +.509 on 4-cov versus +.005 on 2-cov). Design information must never enter the score.
- **Calibration network sees only design descriptors and the score.** It never sees which family
  generated the data, so it cannot learn to detect misfit; its only job is "what does this score
  look like when nothing is wrong, for a dataset shaped like this?"

### Route A — Pivot by augmentation (the cheap baseline; run first)

Train the existing architecture on the wide corpus of §4 with no other change. Ask only: does
breadth alone buy invariance? This is the author's analogy taken literally, and it costs one
training run. It may partly work; it is unlikely to fully work, because nothing in the loss
*rewards* invariance.

### Route B — Amortized calibration (primary; most likely to succeed)

Keep the score network. Add a second network **C** estimating the null distribution of the score
conditional on design descriptors **d**:

- **d** = (n, p, mean/sd of fitted risk, event rate, min/max fitted risk, condition number of X'WX,
  number of distinct covariate patterns, max leverage, share of binary covariates) — all computable
  from the user's own fitted model, none requiring knowledge of truth.
- **C** is trained *only on null datasets* to output the conditional quantile function
  Q(τ | d) of the score, via **pinball (quantile) loss** on a grid of τ, with a monotonicity
  constraint across τ.
- **Use time:** p-value = 1 − F̂₀(S_obs | d_obs), one forward pass. **No bootstrap.**

This is ALFFI's device, applied to GOF. Its great advantage: C is solving a smooth, low-dimensional
regression problem, which is far easier than forcing a 33-input network to be invariant.

### Route C — Pivot by adversarial training (the ambitious arm)

Louppe et al.'s procedure. Score network f, adversary g predicting **d** from f's output on **null
data only**; train f to minimize misfit-detection loss − λ · (adversary's success). Sweep
λ ∈ {0, 0.1, 0.3, 1, 3, 10} and report the whole power/invariance frontier — the trade-off curve
*is* a result, whichever point we ship.

**If Route C succeeds, it is the better paper** (one threshold, no auxiliary network, closest to the
author's vision). **If only Route B succeeds, that is still the tool the author asked for** — the
user still does nothing but supply data.

## §4. The corpus — this is the actual contribution, and where the work is

The current corpus is one point in nuisance space. The new one spans it.

| Axis | Training values | **Held-out values (never trained on)** |
|---|---|---|
| Sample size n | 100, 150, 250, 400, 600, 900, 1400 | **200, 500, 1000** (so C1/G-BOOT stay comparable), **2000** (extrapolation) |
| Covariates p | 2, 3, 5, 8 | **4, 12** |
| Covariate law | independent normal/uniform; AR(1) ρ = 0.5 | **AR(1) ρ = 0.8**, one skewed (log-normal) covariate |
| Categorical | none, or one 3-level factor | **one 5-level factor** |
| Event rate | ~50%, ~15% | **~5%** (the sparse regime this family cares about) |
| Departure families | the 100 random families, regenerated per cell | the 12 confirmatory families |

**Held-out design classes are the whole point.** "Works on any data" is only demonstrated by
coverage on designs absent from training. Training and testing on the same grid proves nothing.

**Volume.** ~250,000 null datasets (calibration network + threshold estimation) and ~250,000
alternative datasets (score network). Sampling of n is weighted toward small n, where the cost is
low and the statistical interest is high.

**Deposit rule.** Every dataset's feature row, design descriptor vector, and score are written to
disk. No repeat of the C1 mistake.

## §5. Architecture — deliberately modest, and why

**Score network.** Keep the 33 features and the small-MLP ensemble. C2 already settled the "does a
deeper reader win" question: a CNN on 50-bin residual curves **lost** twice (.7546 / .7479 versus
.7955) because at these training volumes the fine-resolution input is mostly noise. Do not reopen
it. The one change worth testing is **DeepSets / attention over the 10 grouped residuals**
(permutation-equivariant by construction rather than by the sorted-copy trick) — a clean ablation,
not a headline.

**Calibration network (Route B).** MLP, ~3 hidden layers of 64–128 units, inputs = d (≈10 numbers),
outputs = the score's quantiles at a τ-grid, pinball loss, monotone in τ. Small and easy.

**Adversary (Route C).** MLP predicting each component of d from the score; gradient reversal.

**Why not something bigger.** The binding constraint is not model capacity — it is corpus
generation, which is CPU-bound R. A larger network trains in minutes; the corpus takes hours. Spend
the budget on breadth of data, not depth of model. This is also the honest lesson of C2.

## §6. Compute and the GPU

| Stage | Where it runs | Cost |
|---|---|---|
| Corpus generation (~500k datasets) | **CPU, R, 12 workers** — unavoidable, it is the GOF battery | ~4–6 h with `code/fast_battery.R` (4.4× verified, identical output). Without it this stage would be ~20 h. |
| Score-network training | CPU or GPU | minutes |
| **Route C λ-sweep + adversary** | **GPU** | this is where a GPU actually pays: ~6 configurations × repeated adversarial epochs. CPU ≈ 6+ h; GPU ≈ under 1 h |
| Calibration network | GPU or CPU | minutes |
| Held-out evaluation | CPU (battery again) | ~1–2 h |

**GPU status: the installed torch is `2.11.0+cpu` — CUDA is not available.** The GTX 1660 Ti
(Turing, compute capability 7.5, 6 GB) is fully supported by current CUDA builds. To enable:

```bash
pip install --force-reinstall torch --index-url https://download.pytorch.org/whl/cu124
```

~2.5 GB download, and it replaces the CPU build (only `code/torch_curve.py` from the closed C2 gate
uses it; low risk). **Worth doing for Route C, not needed for Routes A and B.**

**The honest note:** the GPU cannot touch the dominant cost. Corpus generation is 11 classical
statistical tests in R on small matrices — sequential, interpreter-bound, the opposite of what a GPU
accelerates. The 4.4× battery fix bought far more wall-clock than the GPU will.

## §7. The decisive experiment and the pre-declared criteria

**One fixed rule is chosen per route before evaluation. Then, on the held-out design classes only:**

| # | Criterion | Threshold |
|---|---|---|
| P1 | **Coverage on unseen designs** — realized size at the shipped threshold | within [.035, .065] in **every** held-out class, Clopper–Pearson covering .05 |
| P2 | **Uniformity** — pooled null p-values | KS against U(0,1), p > .01 |
| P3 | **Invariance** — spread of realized size across held-out classes | max − min ≤ **.02** (currently .033: .025 → .058) |
| P4 | **Power is not sacrificed** — mean power over the 12 families versus the pre-declared comparator, both at their realized sizes | ≥ −.01 (a tie is a win, because the deliverable is *no calibration needed*) |
| P5 | **Extrapolation** — n = 2000 and p = 12, both outside the training range | size in [.03, .07]; reported separately and never merged into P1 |

**Fallbacks, declared now:**

- **All five pass (any route).** The paper is the tool the author wanted: a pretrained test with a
  fixed threshold, validated on designs it never saw. Strongest outcome. Venue decision reopens —
  JCGS remains the floor, and a clean pivotality result plus theory could justify a reach.
- **P1–P3 pass but P4 fails** (invariant but weaker). Still a real contribution — a *calibration-free*
  test that ties classical power. Ship it, and say plainly that the power advantage did not survive.
- **Route B passes, Route C fails.** Report both. "Learning the threshold works; forcing invariance
  does not" is a clean, useful finding.
- **All routes fail P1–P3.** Then the deployment gap is structural, not a training gap. The paper is
  `05_GBOOT_VERDICT.md`'s warning, now backed by a serious attempt to fix it — which makes the
  warning much stronger, not weaker. **This is a good outcome too, and it must not be spun as failure.**

**Forbidden:** tuning λ, the corpus, or the descriptor set *after* seeing held-out results; reporting
in-training-grid coverage as if it were the claim; quietly dropping a held-out class that fails.

## §8. Risks

| Risk | Why it is real | Mitigation |
|---|---|---|
| The score learns to infer n or p from the features anyway | It already did once — the leakage autopsy | Score network never receives design descriptors; report size by design class always, never pooled |
| Route B's calibration network is right on average, wrong in the tail | Coverage lives at τ = .95, where data is thinnest | Weight the pinball loss toward high τ; report coverage, never mean squared error |
| Adversarial training does not converge | Known failure mode of minimax training | Route A and B are independent of it; the λ = 0 end of the sweep *is* Route A |
| Breadth dilutes power | More nuisance variation means a harder learning problem, and the leakage lesson says these features partly *use* design | P4 exists precisely to detect this; a tie still ships |
| Corpus generation dominates the schedule | It always has | `fast_battery.R` (verified identical), plus small-n weighting |

## §9. Order of work

1. Freeze this document. Write `07_UNIVERSAL_PREDECLARATION.md` fixing the exact architectures,
   seeds, and the held-out grid **before** generating anything.
2. Build and verify the wide corpus generator (extend `families_shared.R` with the new design axes;
   run `verify_fast_battery()` first).
3. Generate the corpus. Deposit everything.
4. Route A (free — it is one training run on the new corpus).
5. Route B — the primary bet.
6. Install CUDA torch; Route C λ-sweep; report the full frontier.
7. Evaluate all routes on held-out classes against P1–P5. Report as it falls.

**Do not write any manuscript text before step 7.** Every previous framing in this project changed
after a measurement.
