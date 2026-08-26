# 23 — GATE SEL: THE LEARNED SELECTOR (pre-declaration, frozen 2026-08-22 before any data)

**Frozen.** Authorised by the author 2026-08-22 ("option 1 and 2"). Written while S3 was still
running, from S3's corpus-1 numbers only; nothing in later corpora may amend it.

## 1. The idea, and why it survives every measured failure

Six constraints were bought with this project's failures (`22` — the thinking is summarised in
MANIFEST phase 23 when this gate closes). The two that killed the deep route:

- **F1** the sensitivity–validity tradeoff (~.07 price) is structural: any score powerful near the
  null manifold can also see the bootstrap-vs-truth gap;
- **F3** training-corpus variance (range .157 on one cell) exceeds every design effect measured.

The selector sidesteps both by a factorisation:

> **Classical statistics measure the MAGNITUDE of misfit** — they are near-LR-optimal in their own
> directions and carry ZERO corpus variance. **The learner infers only the DIRECTION** — *which
> question to ask of this dataset* — a ~12-way choice whose output space is tiny, whose features
> are pattern SHAPES (far nearer to nuisance-free than magnitudes, which is what F1 punishes).

**The measured prize (why this is not hope):** the Cauchy combination loses **.0565 on average**
(per-cell up to ~.08–.14) to its own best member — measured independently in gate PB, gate EM, and
S3 corpus 1. Recovering HALF clears the frozen +.03 bar on the classical tests' own turf.

**The floor (why this cannot repeat S3's failure):** if the selector learns nothing, its weights
are uniform and the statistic IS the Cauchy combination — the pre-declared comparator. The
construction shrinks toward the comparator, never below it, up to selection noise which SEL3 caps.

## 2. The construction (fixed)

**Arms (12):** the eleven classical tests, each **bootstrap-recalibrated** (its exact rank p-value
within its own B=199 bootstrap distribution — this alone frees the conservative components, whose
realized analytic sizes ran .014–.036 on real FDIC cells), plus the K=6 deployment-matched map
score as arm 12 (the location-departure specialist; gates S4/PB).

**Selector:** a small softmax model w(f) over the 12 arms reading DIRECTION features only —
the normalized shape of the residual pattern, never its magnitude: (i) the 12 arms' bootstrap
RANKS transformed to z-scores and then L2-normalized (shape of the evidence profile, not its
size); (ii) the map's top-3 normalized cell coordinates. Trained on the random-departure library
(kind-1 datasets) to predict which arm ranks the dataset most extreme; L2-shrunk toward uniform.

**Statistic:** the weighted Cauchy `C = Σ w_j(f) · tan((0.5 − r_j)π)` where r_j is arm j's
within-block bootstrap rank. **The identical pipeline — features, weights, combination — is applied
to the observed row and to EVERY bootstrap replicate, and C is ranked within its own block:**
`p = (1+#{C_b ≥ C_obs})/(B+1)`. This is gate EM's arm-C device, whose exact level was verified in
all three corpora (sizes .039–.057); the adaptivity is absorbed by the bootstrap by construction.

## 3. Design

n ∈ {50, 100, 200} × designs D1–D4 (S3's, unchanged) × 3 corpora = **36 blocks**, each 2000 nulls
(1000 train / 1000 size) + 1000 training departures + 2000 held-out departures, amplitude
`2·√(100/n)`, B=199, Firth policy — everything identical to S3 except that **every bootstrap
replicate now also carries the fast-battery p-values** (the cost driver: ~6.7 s per dataset).
Estimated wall-clock ≈ 1.5–2 nights.

## 4. Criteria (frozen now)

| # | Criterion | Threshold |
|---|---|---|
| **SEL1** | mean matched-level gain over the NA-robust CCT-11 across all cells | ≥ **+.03** |
| **SEL2** | deployed size ∈ [.035,.065] | every cell |
| **SEL3** | no-harm: gain vs CCT-11 | ≥ **−.01 in every cell** (selection noise cap) |
| **SEL4** | across-seed range of SEL power per (n, design) | ≤ .10; else mean-and-range only |
| **SEL5** (reported, no threshold) | fraction of the CCT→best-of-11 gap recovered | the mechanism number |

## 5. Fallbacks (frozen now)

- SEL1+SEL2+SEL3 pass → this becomes the paper's power section, honestly framed: *the learner's
  productive role is selection among classical questions, not replacement of classical statistics.*
- SEL1 fails but SEL3 holds → the selector cannot read direction at these n; report, keep the
  calibration-led framing, the CCT→best gap stands as an open inefficiency.
- SEL2 or SEL3 fails anywhere → that cell reported; if SEL3 fails the shrinkage was insufficient
  and the construction is reported as not-yet-safe. **No re-tuning of the shrinkage after results.**

## 6. Forbidden

Feeding the selector any magnitude information. Comparing against anything but the pre-declared
CCT-11 (best-of-11 is reported context, not the criterion). Re-running with adjusted shrinkage,
features, or arms after seeing results. Quoting single-seed numbers if SEL4 fails.
