# 00 — VISION: the Deep-GOF paper (target: Q1)

**Status:** vision + gates. Written 2026-08-15, the day the recipe passed its three gates.
**Working title candidates** (decide at build time, SEO rules of `academic-writing` §2):
1. *Pretraining a goodness-of-fit test: simulation-based inference for logistic regression*
2. *A pretrained neural goodness-of-fit test with exact Monte Carlo calibration*
3. *Beyond the battery: reading the residual curve for logistic-regression goodness of fit*

---

## 1. The one-sentence thesis

A goodness-of-fit test can be **pretrained once** on a large randomized library of simulated
misspecifications — reading the model's residual curve alongside the classical battery — and
shipped frozen with Monte Carlo calibration, so that a user hands it a fitted model and
receives a decision and a p-value that are **more powerful than any classical rule at every
sample size studied**, including detections (probit) that no classical test achieves.

## 2. Why this can carry a Q1 journal

1. **A positive, quantified headline** — the thing the LEGofT paper honestly lacked:
   beats the equal-weight Cauchy rule at ALL sizes (+.048/+.024/+.050), with controlled size,
   verified stable (20/20 inits), and generalizing to shape kinds absent from training
   (step keeps 88% of its gain). At n=1000 it exceeds the hindsight selection oracle.
2. **A qualitative first:** probit misspecification — max .078 for every individual test,
   "undetectable" in our own published words — detected at .147 (double the ceiling).
3. **A methodological identity Q1 editors currently want:** this is *simulation-based
   inference / amortized inference* (the SBI literature: Cranmer–Brehmer–Louppe PNAS 2020)
   brought to classical GOF testing — pretrain on simulations, deploy amortized, calibrate
   exactly. That framing places the paper in a hot literature while the application
   (logistic GOF) is classical and clinically motivated.
4. **The negative results are part of the method:** four documented failure modes (design
   leakage ×2, family memorization, location instability) each with its diagnostic — this is
   the "how to do this without fooling yourself" section reviewers reward, and nobody in the
   GOF literature has written it.
5. **A shipped artifact:** `deepgof(fit)` → decision + p-value, in `ebrahim.gof` (the
   deployment shell — battery, Monte Carlo calibration, CRAN pipeline — already exists
   from the LEGofT work).

## 3. Venue ladder (Q1 target per the author)

| Rank | Venue | Why | Risk |
|---|---|---|---|
| 1 | **JRSS-B** (reach) | pretraining + exact finite-sample calibration is methodology; SBI-for-classical-testing is novel there | wants more theory than we have; desk-reject risk |
| 2 | **JASA T&M** (reach) | the positive headline the earlier JASA attempt lacked now exists; BAGofT (the JASA 2023 neighbor) is the natural foil | needs the theory section built out (consistency of the pretrained rule) |
| 3 | **Biometrika** | short-format, method-first; the calibration exactness + holdout generalization is a clean story | length discipline; theory expectations |
| 4 | **JCGS** or **Statistics and Computing** | computational statistics home; simulation-heavy papers welcome | Stat&Comp desk-rejected EDGE on topic fit — JCGS safer of the two |
| 5 | **JMLR / TMLR** (alternative track) | if framed as amortized inference with a statistics application | leaves the biostat audience |

**Recommendation:** build to **JASA T&M** standards (the strictest apparatus we already own —
the `jasa-submission` skill and template exist from the LEGofT campaign), submit there first;
Biometrika second; JCGS the solid landing. Do NOT send to SiM — this is a methods paper, and
LEGofT already owns the SiM lane for this program.

## 4. Contribution list (draft)

1. A **pretraining framework** for GOF tests: a randomized misspecification library
   (the random-family generator), a residual-curve representation, and a deep ensemble —
   with Monte Carlo calibration making validity **independent of the learning** (Theorem
   imported from the LEGofT paper's machinery; cite, don't re-prove).
2. **Power:** dominates the strongest classical combination at every n studied; per-family
   anatomy with honest trades.
3. **The probit detection** — first calibrated detection of near-logistic link
   misspecification at these sizes (verify this claim in the lit sweep; scope carefully).
4. **The failure-mode catalogue** with diagnostics: design-leakage check (split by design),
   ceiling check (a result above its ceiling is a bug), location-instability fix
   (dual coding), memorization fix (diversity).
5. **Software:** `deepgof()` in ebrahim.gof.

## 5. What is NOT claimed

- No claim at n > 1000 (gains shrink as the battery saturates; scalability documented to 10⁴).
- No claim that learning beats equal weights **on the battery alone** (it does not — LEGofT
  settled that; the gain here comes from the residual curve, i.e., MORE INFORMATION, not
  better weighting. This sentence defuses the obvious referee attack).
- Per-family trades (cloglog, sine) are corpus-seed variable — claim means and robust wins.
- The two-component honesty: generic curve-reading (survives kind-holdout) + prior coverage
  (adds on top). Both stated.

## 6. GATES before writing (the LEGofT lessons, paid once, applied here)

| Gate | Requirement | Status |
|---|---|---|
| **C1** | **Confirmatory rerun under a pre-declared protocol**: fresh corpus seeds, fresh eval replicates, n ∈ {200, 500, 1000} declared in advance, headline = the deep ensemble, comparator = CCT-full, success criterion written BEFORE the run (mean gain > 0 at every n, sign consistent). The current numbers are EXPLORATORY | ☐ |
| C2 | **PyTorch stage** (50-bin curves, CNN/Deep Sets): if it beats the nnet ensemble, IT becomes the shipped rule and C1 runs on it instead — decide the architecture BEFORE C1 | ☐ |
| C3 | Literature sweep: SBI/amortized inference (Cranmer 2020), neural ratio tests, C2ST, NPLM, BAGofT/GRASP — position against each; verify the probit-first claim | ☐ |
| C4 | Theory section: validity theorem (imported), consistency sketch, and the honest statement of what is NOT proven | ☐ |
| C5 | Real-data demos: Pima2 + LBW through `deepgof()` (the pipelines exist) | ☐ |
| C6 | `deepgof()` shipped in ebrahim.gof with the frozen ensemble weights | ☐ |
| C7 | Every number in a manifest; every figure distinctive (the EK bar — no default ggplot) | ☐ |

## 7. The figures (planned; distinctive per the standing quality bar)

1. **The money figure:** power at all three n — ensemble vs CCT-full vs best individual,
   with the probit inset ("the undetectable, detected")
2. The residual-curve gallery: what the net sees (curves for quad/cubic/step/probit vs null)
   — this is the paper's visual identity, nobody in the GOF literature draws these
3. The failure-mode panel: 4 small multiples, one per documented failure, each with its
   diagnostic — the "don't fool yourself" figure
4. Kind-holdout: retained gain per family (the 88% bar)
5. The pipeline schematic: library → curves+battery → ensemble → Monte Carlo calibration →
   `deepgof(fit)` → decision

## 8. Relationship to the rest of the family

- **LEGofT (SiM, accepted internally):** owns the battery, the calibration theorem, the
  two-domain localization. This paper CITES it for machinery and never re-claims it.
  LEGofT's honest negative ("learning does not beat equal weights on the battery") is this
  paper's SETUP: the gain required new information, not new weights.
- **EDGES/EDGE:** untouched. The single-self-citation discipline relaxes here (LEGofT will
  be published/preprinted by then), but cite only what is used.
- The random-family **library** itself may be spun out as the standalone benchmark resource
  paper if reviewers ask for it (50–100 documented families was already identified as a
  field-wide need).
