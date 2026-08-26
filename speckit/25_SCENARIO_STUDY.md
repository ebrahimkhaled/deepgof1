# 25 — GATE SCN: THE SCENARIO x SEVERITY STUDY (pre-declaration, frozen 2026-08-22)

**Frozen before any data.** Author-requested (2026-08-22): a scenario-by-scenario power comparison
at several severity levels and exact sample sizes, against the nearest classical competitors.

## 1. Why this gate exists, and what S3 could not answer

S3 compared against a single aggregate comparator (CCT-11) over a RANDOM MIXTURE of departures and
returned one number per cell (mean -.1003). It could not say WHICH departures the learned test wins
or loses on, or at what severity, because the mixture averages them together. Gate PB then showed
the answer is family-dependent (link departures: decisive loss), and the S3 by-design breakdown
showed it is representation-dependent. **This gate resolves the mixture into its components.**

## 2. ⛔ COMPARATOR SET — `A.poly` and `A.spline` are EXCLUDED

Those two are `cdef.gof`, i.e. the author's own EDGE/EDGES family, still unpublished preprints.
Benchmarking this paper against them compares a preprint to a preprint and is not a meaningful
external comparison. **They are excluded from every criterion and every figure in this gate**, and
their exclusion is stated in the caption. *(For the record, in S3 they were the strongest arms:
A.spline .3878, A.poly .3718 — that fact is reported in the text, never used as the comparator.)*

**The NINE comparators**, in two mechanistically distinct families declared here:
- **PARTITION** (bin the fitted risk): `HL`, `HL-equalwidth`, `Pigeon-Heyse`, `Tsiatis`, `Xie`
- **DIRECTED / SMOOTH IN THE INDEX**: `EF`, `DEF.poly2`, `DEF.poly3`, `DEF.stukel`
- plus their **Cauchy combination of the nine** (`CCT-9`), which is the pre-declared HEADLINE
  comparator for the criteria below.

## 3. Scenarios — declared in two blocks BY MECHANISM, not by result

The split is the paper's own pre-existing claim (phase 11, Figure 3): a covariate-space partition
finds departures that occupy a REGION; risk-binning averages them away.

- **BLOCK A — localised in covariate space** (the map's predicted home): `bump`, `step`,
  `interaction`, `abs`
- **BLOCK B — functions of the linear predictor / smooth** (the directed tests' home): `indexq`,
  `indexs`, `quad`, `cubic`, `sine`

Generators are the existing `make_family`/`gen_random` term kinds, one kind per scenario, unchanged.

## 4. Severity — an interpretable, scenario-independent axis

`sd(delta)` on the LOGIT scale, which is exactly the `amp` parameter the generators already use.
**Severity s in {0, 0.3, 0.6, 1.0, 1.5, 2.2}**; s = 0 IS the null and yields each arm's SIZE, so
every power curve begins at its own measured level. **Severity is NOT rescaled by n**, so power
visibly rises with n down each column — that comparison is the point.

## 5. Design

n in {50, 100, 200} · 9 scenarios · 6 severities · **2 corpora** (101000, 102000) ·
1000 evaluation datasets per (n, scenario, severity, corpus) · 2000 nulls per (n, corpus)
(1000 train / 1000 size) · 1000 random-departure training datasets per (n, corpus) ·
B = 199 · K = 6 (gate S4) · deployment-matched training, mix pinned at 1.0 ·
randomized matched level at alpha = .05 · Firth separation policy · design D1 (p = 2, textbook).

**Two corpora, not three:** S3 established that grid MEANS agree across corpora to .008 while
single cells range .0437. This gate reports per-cell curves, so the across-corpus range IS DRAWN as
a ribbon on the learned arm and no single-corpus number is quoted.

## 6. Criteria (frozen now)

| # | Criterion | Threshold |
|---|---|---|
| **SCN1** | mean matched-level gain over **CCT-9** across BLOCK A scenarios | > **+.03** |
| **SCN2** | mean gain over the **PARTITION family's best** across BLOCK A | > **0** |
| **SCN3** | deployed size (severity 0) | in [.035,.065] in every (n, corpus) |
| **SCN4** | BLOCK B is REPORTED, not required to win | no threshold — the mechanism predicts losses here |
| **SCN5** | across-corpus range of the learned arm | drawn on every panel; if > .10 anywhere, that panel is annotated |

## 7. Fallbacks

- SCN1 + SCN2 pass -> a bounded, mechanistically predicted positive claim: *the pretrained
  covariate partition beats the risk-partition family on departures localised in covariate space.*
  This becomes the paper's power section, scoped to BLOCK A and stated with BLOCK B's losses.
- SCN1 passes, SCN2 fails -> the gain is over the COMBINATION only, not over the best partition
  test; report as such, no family-level claim.
- Both fail -> the mechanism does not deliver even on its predicted home ground. Report, and the
  paper is calibration-and-lessons only.

## 8. Forbidden

Re-drawing the scenario blocks after seeing results (they are fixed in section 3 by mechanism).
Adding or removing comparators after results. Quoting `A.poly`/`A.spline` as comparators.
Quoting a single-corpus number. Selecting severities post hoc.

## 9. Execution

`code/scn_grid.R` -> `code/scn_score.py` -> `code/figs/fig7_scenarios.py` + `code/scn_verdict.py`.
Deposit `results/scn/`.
