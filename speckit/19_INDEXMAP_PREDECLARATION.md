# 19 — GATE EM: THE INDEX-ROTATED MAP (pre-declaration, frozen 2026-08-21 before any data)

**Frozen.** Reported against verbatim. Authorised by the author on 2026-08-21: *"no problem to have
two versions of the test if the other could have specific power in the link specification."*

## 1. The question

Gate PB established that the pretrained test does not detect link misspecification, and named the
mechanism: the residual map grids over covariate RANKS and finds departures that occupy a REGION of
covariate space, whereas a link change is monotone in the linear predictor and occupies no region.

That leaves one question a referee will certainly ask, and which gate PB cannot answer:
**is the boundary a property of LEARNED GOF TESTS, or only of THIS GRID GEOMETRY?**

Both answers improve the paper, which is why this gate is worth running and why it is not a rescue:

| outcome | what the paper gets |
|---|---|
| the index-rotated map DETECTS link departures | the boundary is **geometric, not methodological** → a design principle: *the map's axes decide which departures it can see*. Strictly more general than a limitation. |
| it does NOT | the limitation becomes **bulletproof** — the obvious alternative was pre-declared, run, and reported. The referee's question dies before it is asked. |

## 2. The construction (fixed)

**Arm A — covariate map (unchanged).** K x K grid over (rank x1, rank x2). The shipped rule.

**Arm B — index-rotated map (NEW).** The SAME cell statistic and the SAME standardisation, on a grid
rotated onto the fitted index: axes are (rank eta-hat, rank v), where eta-hat = b1*x1 + b2*x2 is the
fitted linear predictor and v = b2*x1 - b1*x2 is the orthogonal direction in the design. This is a
rotation of covariate space aligned to the fitted model, so the index becomes the primary axis.
Ranks keep it scale-free; nothing else about the construction changes.

**Arm C — the combined test (the one that would SHIP if both are kept).**
⚠️ **This is the honesty-critical piece.** Two geometries must NOT mean "report whichever won."
Both maps are scored on the observed dataset AND on every bootstrap replicate; their two scores are
combined into a SINGLE statistic by the Cauchy combination of their within-bootstrap rank
p-values; that combined statistic is then ranked against its own bootstrap distribution:

    p_C = (1 + #{ C_b >= C_obs }) / (B + 1)

**The level is exact by construction and needs no multiplicity correction**, because the same
combination rule is applied to the observed data and to every replicate alike. This is the only
form in which a two-geometry test may be reported.

## 3. Design

- **Cells:** n = 200 ONLY — the decisive cell, where gate PB's gap is widest (-.138 on cloglog) and
  the signal strongest. Probit is excluded from the criteria: gate PB established it is a blind spot
  for all twelve arms, so it can discriminate nothing. It is still generated and reported.
- **Seeds:** the same three corpora, 101000 / 102000 / 103000.
- Everything else identical to gate PB: 2000 nulls (half train / half deployed size), 1000 random
  departures as the training alternative class, 2000 each of probit / probit-extreme / cloglog as
  held-out KINDS, B = 199, deployment-matched training with **mix pinned at 1.0**, randomized
  matched-level comparison at alpha = .05, NA treated as a failed test that cannot reject.
- **Arm B and Arm C are trained on the same random-departure library as Arm A.** No link family
  enters training, in either geometry. This remains a held-out-kind evaluation.

## 4. Criteria (declared now, on the cloglog cell at n=200)

| # | Criterion | Threshold |
|---|---|---|
| **EM1** | Arm B matched-level power exceeds Arm A | > +.05, mean over 3 seeds AND in >= 2 of 3 seeds |
| **EM2** | Arm B deployed size | in [.035, .065] at every seed |
| **EM3** | Arm C size | in [.035, .065] — the combination must not break the level |
| **EM4** | Arm C power >= max(Arm A, Arm B) - .02 | i.e. combining costs at most 2 points against the better single geometry |
| **EM5** | Across-seed range of Arm B power | <= .10, else mean-and-range only, no single-seed number |

## 5. Fallbacks (declared now)

- **EM1 + EM2 pass** → the boundary is GEOMETRIC. Report both geometries; if EM3 + EM4 also pass,
  **Arm C is the shipped test** and the paper gains a design principle. Section 18's boundary prose
  is rewritten accordingly (the blind-spot half about probit is unaffected and stands).
- **EM1 fails** → the boundary is METHODOLOGICAL and deep. Arm B is reported as a tested-and-rejected
  alternative in one sentence; section 18 stands as drafted and becomes much stronger for it.
- **EM2 or EM3 fails** → that geometry is not validly calibrated; report and claim nothing from it.
  A geometry that cannot hold its level does not ship regardless of its power.

## 6. Forbidden

Selecting between Arm A and Arm B per cell after seeing results — that is the multiplicity Arm C
exists to handle. Quoting Arm B power on probit as evidence of anything (blind spot, established).
Re-running with a third geometry and reporting it under this gate. Adding seeds to move a criterion.
Reporting Arm C without reporting that its level was verified (EM3).

## 7. Execution

`code/index_map.R` (both maps per dataset, one block per seed) -> `code/index_map_score.py`
(three arms, shared `probit_common.py` matched-level device) -> `code/index_verdict.py`.
Deposit `results/indexmap/` — every per-replicate p-value for all three arms, all seeds.

---

# RESULT (2026-08-21, all 3 corpora — criteria above were frozen before any data)

**EM1 FAIL · EM2 PASS · EM3 PASS · EM4 FAIL · EM5 PASS.**
**The boundary is NOT geometric.** Rotating the map onto the fitted index does not reliably help.

## cloglog at n=200, matched level

| seed | A covariate | B index | C combined | best of 11 | B - A | C - max(A,B) |
|---|---|---|---|---|---|---|
| 101000 | .0219 | .0674 | .0302 | .2485 | **+.0455** | -.0372 |
| 102000 | .1259 | .1603 | .1427 | .2535 | **+.0343** | -.0175 |
| 103000 | .1786 | .0698 | .1287 | .2295 | **-.1088** | -.0500 |
| **mean** | **.1088** | **.0992** | **.1005** | **.2438** | **-.0097** | — |

⭐ **THE DECISIVE OBSERVATION, and the reason this gate answers a different question than it asked.**
Arm A alone swings **.0219 - .1786 across corpora, a range of .1567** — LARGER than any B-A
difference observed (+.046, +.034, -.109). **The corpus-to-corpus noise exceeds the geometry effect,
so a geometry effect is UNIDENTIFIABLE at three seeds.** The first two corpora looked like a clean
+.04 gain; the third reversed it by -.109. Anyone reading two corpora would have reported a
tripling that does not exist.

**Consequence for the manuscript, and it generalises beyond this gate:** *hyper-parameter and design
comparisons for a pretrained test cannot be made at the seed counts normally used, because the
training-corpus variance dominates.* Distinguishing a .04 effect against a .157 range needs on the
order of dozens of corpora, not three. This is a methodological finding in its own right and belongs
beside the instability paragraph in the limitations.

## The fallbacks that fire

- **EM1 fails** -> section 5, second bullet, verbatim: *"the boundary is METHODOLOGICAL and deep.
  Arm B is reported as a tested-and-rejected alternative in one sentence; section 18 stands as
  drafted and becomes much stronger for it."* Even the best single corpus put arm B at .0674
  against DEF.poly2's .2485 — the rotation never came close to closing the gap.
- **EM4 fails in all three seeds** (-.037, -.018, -.050). **The two geometries are NOT freely
  combinable**: the Cauchy combination costs .02-.05 against the better single geometry, because it
  dilutes when only one component carries signal. This is the SECOND independent measurement of that
  dilution — gate PB measured CCT-11 at .1037 against its own best component's .2438.
  ⛔ **A two-geometry test does not ship as an automatic combination.** If both geometries are ever
  offered, it is as a PRE-SPECIFIED choice by the analyst's stated concern, never as a blend and
  never as a post-hoc pick.
- **EM2, EM3 pass** — both geometries and the combination hold their level (arm B .047/.037/.053;
  arm C .057/.039/.051). The construction is VALID; it is simply not more powerful. Arm C's exact
  level is worth keeping as a demonstration that the exchangeable-combination device works as the
  theory says, even though the combination is not adopted.
- **EM5 passes** at .0929, but only barely, and against a .10 bar set before we knew arm A itself
  ranged .157. Read it as uninformative rather than reassuring.

## Reproducibility check (unplanned, and the strongest one in the project)

Gate EM regenerated its data from scratch with a different script. Arm A on cloglog read
**.0219 / .1259** for seeds 101000 / 102000; gate PB read **.0226 / .1249** on the same seeds and n.
Independent regeneration, same numbers to within .001 — the whole pipeline reproduces.

## What the manuscript takes from this gate

One sentence in the boundary section: *the obvious alternative geometry — a grid rotated onto the
fitted linear predictor — was pre-declared, run on three corpora, and did not reliably improve
detection of link departures; the training-corpus variance exceeded the effect.* That sentence kills
the referee question for the cost of a line, which is exactly what the gate was bought for.
