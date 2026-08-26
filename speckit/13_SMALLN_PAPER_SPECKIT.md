# 13 — THE SMALL-SAMPLE PAPER: full spec-kit (2026-08-21)

**The author's reframing:** build the paper on n = 50–200, where classical GOF tests cannot lean on
their asymptotics. This document is the executor's plan for that paper — written for handoff to a
second agent (Opus) who will run the remaining gates and write the manuscript. It supersedes the
venue/thesis sections of `11_MASTER_PLAN_V2.md`; the mechanism/diagnostic/construction sections
there remain the technical foundation.

---

## §0. Why this scope is right — the evidence, assembled

Every measurement in sixteen manifest phases points to the same boundary:

| evidence | value | source |
|---|---|---|
| deployed gain at n=200 | **+.022** | ladder (phase 16.1) |
| matched-level gain at n=200 | **+.034** | pocket analysis (phase 15) |
| deployed gain at n=400 / 800 / 1600 | −.005 / +.006 / +.001 | ladder — the advantage dies as n grows |
| selection-oracle headroom at n=200 vs n=1000 | **+.067 vs +.017** | phase 0 — headroom grows as n shrinks |
| classical size on a realistic small-n design (n=250, compound) | **.1008 — double nominal** | hard cell (phase 16.3) |
| ⭐ matched-level gain, same compound cell at MID-POWER | **+.0710 (7 SE)** — learned .3833 vs classical .3123 | phase 17.1 — the pockets stack once off the power floor |
| learned size, same cell | .0695 | hard cell |

And the theory agrees from both sides: classical GOF tests are asymptotically near-optimal in
their directions (so beating them at large n is capped — measured, not assumed), while every one
of their references is a large-sample approximation (χ², normal) that degrades exactly on this
range. **A test trained at finite n and calibrated by simulation has no asymptotic step to fail.**
The scope is not a retreat; it is the regime the method is *for* — pilot studies, rare-disease
cohorts, subgroup analyses, phase-I/II trials.

## §1. The thesis, in the form it will be defended

> At the sample sizes where logistic regressions are most often over-interpreted (n = 50–200),
> the standard goodness-of-fit tests neither hold their stated level nor deliver their nominal
> power. A test **pretrained on simulation, with its null class matched to its deployment
> calibration**, holds level at any n by construction and matches or exceeds the classical battery
> where the battery's own approximations fail.

Claims ledger (each falsifiable, each with its gate):
- **T1** Classical analytic references are materially miscalibrated at n ≤ 200 → gate S1
- **T2** The deployment-matched learned test holds level across n = 50–200 and across designs → S3
- **T3** Its matched-level power meets or exceeds classical, with the gap growing as n falls → S3
- **T4** The mechanism (bootstrap-vs-truth separability) and its one-number diagnostic → already measured (phases 13, 15)
- **NOT claimed:** superiority at n ≥ 400 (measured tie); "exceeds a lot" anywhere until S3 says so.

## §2. The construction that ships (fixed; every piece has a measurement)

1. **Residual map** over covariate **ranks** (scale/distribution-free), cell = Σ(y−π̂)/√(Σπ̂(1−π̂))
2. **CNN score**, small (the C2 lesson: depth does not pay; resolution and geometry do)
3. **⭐ Deployment-matched training** — the null class IS bootstrap replicates (phase 15:
   AUC .5011, size .052, zero penalty machinery). *Train against the null you calibrate against.*
4. **Distillation** of the 11-test battery as teacher only (+.013), never an inference input
5. **Per-dataset parametric bootstrap**, exact rank p-value, B = 199 (exact .05 grid)
6. Optional reported variant: **Cauchy-combined with classical** (valid under arbitrary dependence;
   confirmed +.021 at 2.7 SE but below our +.025 bar — report as secondary, claim nothing)

## §3. Gates for the executor (criteria frozen HERE, before results)

### S1 — the classical breakdown table *(RUNNING — `smalln_breakdown.R`)*
Analytic size of all 11 tests + CCT at n ∈ {50,75,100,150,200} × event rate {50%,15%}, 4000 nulls
per cell, plus NA/failure and separation rates. **This is the paper's motivating Table 1.**
No criterion — it is a measurement of the premise. If classical calibration is FINE at small n,
the premise dies and the paper reverts to `11_MASTER_PLAN_V2.md` outcome 2. Report either way.

### S2 — n=100 deployed head-to-head *(RUNNING — the recon chain)*
Deployment-matched map vs CCT, both from the same B=199 refits, mid-power amplitudes (×2).
**GO criterion: gain ≥ +.03 at n=100** (the n=200 gain was +.022; the thesis predicts growth as
n falls). GO → build S3. NO-GO (gain ≤ +.01) → the small-n thesis fails its first extrapolation;
report, and the paper reverts to the price-of-validity framing.

### S3 — the confirmatory grid (the paper's main table; ~2 nights of compute)
n ∈ {50, 75, 100, 150, 200} × designs {textbook p=2; p=3 with one 3-level factor; p=3 with
ρ=.7; 15% events} × 2000 nulls + 2000 alternatives per cell, amplitudes scaled ∝ 1/√n to hold
mid-power. One deployment-matched rule per n (or one n-pooled rule — decide by a 2-cell pilot,
pre-declare before the full grid). **Criteria:** size ∈ [.035,.065] in every cell; matched-level
gain > 0 in ≥ 80% of cells; mean gain ≥ +.03; gain at n=50 > gain at n=200 (the trend claim).
**Fallbacks written now:** all pass → the paper as thesised; size holds but mean gain ∈ (0,.03) →
the paper leads with CALIBRATION (T1+T2) and reports power parity honestly; size fails anywhere →
that cell is reported and the scope shrinks to where it holds.

### S4 — tiny-n engineering — **COMPLETE 2026-08-21. SELECTED K = 6.**

> **S3 RUNS AT K = 6 (36 cells), NOT K = 12.** Full result: `20_S4_PREDECLARATION.md`.
> The shipped 144-cell grid was LAST in all four cells (deficit ~.07, ~7 SE) and is 47-69% empty at
> these n; {4,6,8} are statistically indistinguishable and K=6 follows from the frozen tie rule.
> Separation policy: Firth refit when `any(|beta|>8)` or non-convergence, applied to the observed
> fit AND every bootstrap replicate — rate is 0.03%, not the 1.14% originally quoted.
> ⚠️ Gate S2's +.090 was measured at K=12 AND at size .062 vs .056 (not matched level). Those two
> corrections pull in opposite directions — **quote no revised headline until S3 measures both.**

#### (original S4 text, superseded)
At n=50: map resolution K ∈ {4,6,8,12} (144 cells over 50 points is sparse — measure, don't
assume); grouped-feature G ∈ {5,10}; **separation policy**: measure frequency (S1 gives it);
policy = Firth refit for BOTH arms when |β̂|>8, pre-declared, applied identically.

### S4b — GATE S-THEORY (added 2026-08-21; do this in parallel with S3's compute)

Three results, in priority order — proofs to the SUPPLEMENT, statements in the body:

1. **Proposition 1 (exchangeability = validity; ~days).** The exact-rank bootstrap p-value is
   exactly Uniform iff (S_obs, S_1..S_B) are exchangeable under H0; the bootstrap breaks
   exchangeability only through beta-hat; hence if the score's distribution is invariant to the
   bootstrap-vs-true null flavour (distributional AUC = .50), the deployed test is EXACTLY sized.
   **This turns the fig-4 AUC diagnostic into the theorem's one measurable condition, and
   deployment-matched training into the procedure that enforces it.** Cite Barnard (1963) for the
   exchangeable-rank device; state the condition honestly as sufficient, empirically approximated
   (measured AUC .5011, not .5000).
2. **Proposition 2 (Pi-optimality; ~days).** The training library is a prior Pi over departures;
   the score estimates the Pi-mixture likelihood ratio on the map; Neyman-Pearson gives approximate
   Pi-average-power optimality among map-measurable tests (Wald 1943; King 1988;
   Andrews & Ploberger 1994). Converts per-family losses into corollaries.
3. **Consistency (supplement-grade; 1-2 weeks).** Power -> 1 against fixed alternatives with
   non-degenerate map signal under standard GLM regularity. Skip any regret/finite-sample
   optimality bound — vacuous at these n (standing decision from the 2026-08-17 review).

### S5 — GATE PB: LINK MISSPECIFICATION — **COMPLETE, FAILED** (2026-08-21)

> **PB1 FAIL · PB2 FAIL · PB3 FAIL (cloglog) · PB4 PASS.** The showcase does not exist: on link
> departures the classical directional tests win, by more as n grows (−.1375 at n=200 cloglog).
> Full result table and the fired fallback: `17_PROBIT_PREDECLARATION.md`, RESULT section.
> **Executor: do not re-open this gate.** The manuscript prose is ALREADY DRAFTED in
> `18_LINK_SECTION_DRAFT.md` — drop it in, verify the numbers against `results/probit/VERDICT.csv`,
> and move on. It carries TWO findings of different kinds and both must appear: probit is a
> **shared blind spot** (no test of twelve beats chance at n<=200 — a statement about the problem),
> cloglog is a **real boundary** (directional tests reach .244, ours .106 — a statement about our
> method). Publishing the first without the second is selective reporting.

**Full pre-declaration: `17_PROBIT_PREDECLARATION.md`. Read it before touching any PB result.**

The author's proposal: a figure and simulation showing the classical battery failing to detect a
probit misfit while the pretrained test succeeds — eleven famous tests against one. Scientifically
it is the strongest single demonstration available, because probit-vs-logit is the acknowledged
*hardest* link alternative (both links symmetric, differing only in tail weight).

**This is the claim that was RETRACTED on 2026-08-17** (`03_C1_DEFECTS.md` D4). Gate PB re-asks it
honestly: three corpus seeds with the spread reported, deployed bootstrap calibration, probit and
cloglog held out of training entirely, criteria frozen before the data, and a randomized
matched-level device so no arm is credited for over-rejecting.

**For the executor:**
- **Do not quote any pre-2026-08-21 probit number.** They are retracted. Only `results/probit/`.
- **PB1 may fail.** `DEF.stukel` is purpose-built for link departures while the map is a location
  detector; a loss here is a legitimate outcome and goes in the boundary section, not the bin. The
  expectation was recorded before the result specifically so it could not be revised afterwards.
- **`probit_extreme` is a CONSERVATISM cell, not a power cell** — the timing pilot showed the
  battery rejecting .000 there because near-deterministic fitted probabilities make these tests
  conservative. Never present it as power.
- **mix is pinned at 1.0.** mix 0.5 is deposited for information; no claim may cite it.
- If PB3 (across-seed range ≤ .10) fails, **no single-seed number may be quoted anywhere** — mean
  and range only, and the figure moves to the supplement.
- Figure 5 already exists (`code/figs/fig5_probit.py`) and draws the across-seed range as its
  falsification object. Do not remove that band to make the picture cleaner.

**Where it sits in the manuscript.** If PB1–PB4 pass, this is the paper's *demonstration* figure
and earns roughly one page inside the S3 simulation budget — it does not get its own section, and
the page budget in S6b does not move. If it fails, one honest paragraph in the boundary section.

### S6b — PAGE BUDGET AND THE TWO SUBMISSION FORMATS (added 2026-08-21)

The MANUSCRIPT is <= 35 pp double-spaced (JASA hard rule; healthy discipline for EJOR too):
Intro+Table 1 (4pp) · method (5pp) · Props 1-2 statements (2pp) · simulation S3, ONE figure + ONE
table (5pp) · FDIC application (5pp) · mechanism/diagnostic (3pp) · discussion+boundary (3pp) +
references. Proofs, per-cell tables, the ablation history, and the full breakdown grid -> the
SUPPLEMENT. **The 19-phase MANIFEST is the reproducibility archive, never manuscript text.**
Build BOTH shells at S6: EJOR (primary) and JASA Applications & Case Studies (the application
LEADS there — use the `jasa-submission` skill's working template; double-anonymous: cite the
ebrahim.gof family neutrally in third person, NEVER withhold a citation — the ADAC lesson).
Final venue decision at the S6 checkpoint with both PDFs in hand.

### S5 — real data
Two small clinical datasets run through the full pipeline (candidates: the ICU study n=200 —
already vendored in the package work, 198 covariate patterns; a pilot-sized subset analysis of
Pima). Report whatever they show.

### S6 — the manuscript
Sections: Intro (the small-sample gap; Table 1 = S1) → The test (§2, each piece with its
measurement) → The mechanism + diagnostic (phases 12–15 — this is the methods contribution) →
Simulation (S3; **size before power**, per-design, never pooled) → Real data → Honest boundary
(n ≥ 400: measured tie; the ladder table goes IN the paper). **Figures** (all data exists or
comes from S1–S3): the breakdown table as a dot-plot with nominal band; the map itself
(null vs misspecified, one glance); the AUC diagnostic curve; the ladder/price-of-validity;
the main grid forest. House style: EK palette, no default ggplot (standing rule).

## §3b. THE FINANCE APPLICATION (added 2026-08-21) — see `14_FINANCE_APPLICATION.md` in full

**The problem:** PD-model calibration backtesting happens at rating-grade level on small samples
(EBA data: 5th percentile of obligors per grade = 34 — **agent-verified, Opus must confirm**), the
Basel-documented model-level test is Hosmer–Lemeshow (BCBS WP14 p.52; MATLAB/PDtoolkit implement
it), its small-sample unreliability is acknowledged qualitatively but **never quantified** — we
quantify it and supply the replacement. ⚠️ Referee-proofing rule: NEVER say "regulators mandate
HL" — the ECB grade-level instruction is the Jeffreys test; we COMPLEMENT it at model level.
**Primary dataset (endpoints re-verified live by hand):** FDIC BankFind API — Georgia banks at
2008-12-31 (n=335) vs failures 2009–2011; pre-declared 5-covariate CAMELS-proxy logit; replication
= every state with ≥15 failures, ALL reported. Backup: German credit by loan purpose;
multiplicity engine: Taiwan cells 50≤n≤300. **Both branches pre-declared** — if the real-design
bootstrap does not reproduce the simulated size distortion, that transfer failure is reported
per cell.
**S5 is upgraded accordingly:** the finance cells REPLACE the Pima/ICU real-data plan as primary
(ICU stays as a secondary biostat example if space allows).

## §4. Venue — REVISED by the finance research

1. **EJOR — submit first** (Q1, IF ~6.0, "Stochastics and Statistics" track, the strongest
   credit-scoring tradition; the LDP framing defuses its pure-statistics desk rule). No Egypt APC
   deal → subscription route + arXiv green OA (zero cost).
2. **IJF** — if the final application reads as PD forecasting; mandatory code+data suits us.
3. **ESWA** — ML-friendly fallback. (The biostat ladder — Biometrical Journal/SMMR/JCGS — remains
   the alternative if the finance framing is dropped; decision at S6.)

## §4-old. (superseded)

SiM is off-limits (LEGofT's lane — standing decision). Options, with my recommendation first:
1. **Biometrical Journal** — small-sample methodology is house style; mandatory audited
   reproducibility is a STRENGTH for us (every p-value already deposited); Q1; skill exists.
2. **SMMR** — clinical-methods audience that lives at these sample sizes.
3. **JCGS** — the computational home; safe fallback, less clinical reach.

## §5. Standing rules (unchanged, violations already paid for)
Oracle numbers are never results · both arms get the same calibration on the same replicates ·
power ≥ size asserted in code · size per design, never pooled · every per-replicate p-value
deposited · every trial into the MANIFEST with numbers and verdict, failures included ·
**verify every sed/path edit by grep before launch** (three silent failures this project).

## §6. Handoff note for Opus
Read first: `MANIFEST.md` (17 phases — the full history including every retraction),
`11_MASTER_PLAN_V2.md` (mechanism + construction), this document. The recon results land in
`logs/hard_mid_verdict.log`, `logs/smalln_breakdown.log`, `logs/n100_depmatch.log`. Execute S4 →
S3 → S5 → S6 in that order, honouring every frozen criterion; record everything in the MANIFEST
as new phases. The one unbreakable rule: **a criterion that fails is reported as failed.** This
project has retracted three headline numbers by its own audits; the record of that discipline is
part of the paper's credibility.
