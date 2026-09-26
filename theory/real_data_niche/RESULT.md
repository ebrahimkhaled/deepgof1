# Real-data check of the DeepGOF-1 niche: RESULT

Date: 2026-09-25. Design and reading rule: `PREDECLARATION.md`, hashed at **2026-09-25 04:49:40 EDT**
(`PREDECLARATION.sha256`, sha256 `f016bc63...`) before any goodness-of-fit test on any candidate. The
calibration output and the Part B code were hashed at 05:00:12 EDT, before any Part B study replicate
(`RUN.sha256`). Both hash files were re-checked after the run: all files match.
Full tables: `result_tables.md`, `summary_size.csv`, `summary_power.csv`, `summary_diff.csv`,
`partA_pvalues.csv`. Per-replicate p-values: `rows/` (800 files, 16,000 rows).

## Verdict (pre-declared, primary): NOT CONFIRMED -- NEGATIVE

On the two real designs (myopia and Pima), in the niche outcome model, DeepGOF-1 was above every
calibrated comparison arm by more than 2 paired SE in **0 of 12** cells, and in **12 of 12** cells at
least one comparison arm was above DeepGOF-1 by more than 2 paired SE. DeepGOF-1 itself held the level
(size .048 and .057, SE .007). The same result holds per data set (0 of 6 each), at raw level, and in the
dense control. **Recommendation: do not add a real-data niche section to the paper.** The simulated niche
(`theory/niche/`) stands as a simulation result about ten independent Gaussian covariates; this study
found no real covariate design among the installed data sets on which it carries over.

## Row order (the known tie-breaking defect)

Every real file was permuted once before any fit (`set.seed(20260925 + k)`; `perm_<data>.csv`). Mroz is
stored sorted by the outcome. Part A also reports DeepGOF-1 on 20 further random orders (declared
secondary). In Part B the outcome is simulated independently of the fixed row order, so the defect cannot
leak y there.

## Candidate list (a)

201 data objects scanned; 40 passed the coarse metadata screen; three met all criteria:
**myopia** (aplore3; n = 618, 81 events; 12 covariates, 8 continuous), **Pima** (MASS Pima.tr + Pima.te;
n = 532, 177 events; 7 covariates, all continuous), **Mroz** (carData; n = 753, 428 events; 7 covariates,
3 continuous). Reasons for each exclusion are in `PREDECLARATION.md` (a). No enriched (published
interaction) model could be verified for any of the three, so Part A has main-effects models only.
Part B used myopia and Pima (most covariates; Pima over Mroz by the declared tie-break).

## Part A: real outcomes (descriptive; the truth is unknown)

| test | myopia | Pima | Mroz |
|---|---|---|---|
| **DeepGOF-1** (B = 199, K = 6, permuted rows) | .770 | .080 | .205 |
| DeepGOF-1, 20 other row orders: median (range), # p <= .05 | .405 (.010-.935), 4 | .048 (.020-.095), 13 | .463 (.050-.965), 1 |
| HL | .737 | .702 | .001 |
| Stukel SstBoth (marginal) | .400 | .002 | .001 |
| Stukel joint (secondary) | .463 | .015 | .006 |
| Stute-Zhu (B = 200) | .870 | .215 | .000 |
| GiViTI | .423 | .012 | .010 |
| mHL (F-test) | .960 | .102 | .000 |
| Osius-Rojek | .601 | .209 | .558 |
| projection, raw X (B = 999) | .204 | .183 | .031 |
| projection, standardised X (secondary) | .746 | .028 | .001 |
| BAGofT, p.value (fast engine) | .190 | .040 | .000 |
| BAGofT, p.value3 (secondary) | .070 | .150 | .000 |

DeepGOF-1 axes: myopia (DADMY, MOMMY), two binary covariates; Pima (glu, bmi); Mroz (k5, age), k5 a count
with 4 values. At .05: myopia, no test rejects. Pima, Stukel (both forms), GiViTI and BAGofT reject (and the secondary standardised projection, .028);
DeepGOF-1 does not on the declared order (.080) but does on 13 of 20 other orders. Mroz, every test except
DeepGOF-1 and Osius-Rojek rejects. On myopia and Mroz the DeepGOF-1 p-value moves over almost the whole
unit interval with the row order, because both selected axes are heavily tied: on such designs the test is
partly a randomised test until the tie-breaking fix lands. Runtimes: DeepGOF-1 1-2 s, projection (B = 999)
5-15 s, BAGofT 290-325 s, the rest under 1 s.

## Part B: real-design simulation (confirmatory)

**Designs.** myopia: a = SPORTHR, b = ACD (standardised; fitted b_a = -.362, b_b = .346); Pima: a = glu,
b = bmi (b_a = 1.095, b_b = .569). Calibrated C (HL pilot, `calibration.csv`): myopia niche inter
1.4/2.2/3.0, quad .5/.6/.7; Pima niche inter .4/.6/.7, quad .2/.3/.5; dense in `scenarios.csv`.
Replicates: 1,000 per null, 500 per alternative; 16,000 datasets; no fit needed the Firth fallback.

**Size at .05 (1,000 nulls, SE about .007; * = not calibrated, above .0707)**

| data, design | DG | HL | Stukel marg. | Stukel joint | Stute-Zhu | GiViTI | mHL | Osius-Rojek | proj |
|---|---|---|---|---|---|---|---|---|---|
| myopia niche | .048 | .036 | .012 | .033 | .062 | .040 | .028 | .117* | .033 |
| Pima niche | .057 | .033 | .068 | .052 | .070 | .048 | .462* | .064 | .045 |
| myopia dense | .040 | .047 | .050 | .039 | .063 | .038 | .050 | .050 | .029 |
| Pima dense | .045 | .053 | .066 | .056 | .064 | .051 | .692* | .055 | .046 |

Stukel marginal is undefined (NA, so p = 1) on 68 % of myopia niche nulls (no fitted risk above .5), hence
its low size. mHL is badly inflated on the Pima design (.46, .69): reported, not investigated further here.

**Matched-level power, niche design (SE .006-.022)**

| data | family | C | DG | HL | Stukel marg. | Stukel joint | Stute-Zhu | GiViTI | mHL | Osius-Rojek | proj |
|---|---|---|---|---|---|---|---|---|---|---|---|
| myopia | inter | 1.4 | .160 | .146 | .070 | .222 | .155 | .180 | .146 | .668 | 1.000 |
| myopia | inter | 2.2 | .442 | .270 | .154 | .360 | .322 | .318 | .280 | .838 | 1.000 |
| myopia | inter | 3.0 | .749 | .274 | .332 | .368 | .338 | .336 | .296 | .868 | 1.000 |
| myopia | quad | .5 | .540 | .164 | .246 | .368 | .246 | .348 | .180 | .428 | .938 |
| myopia | quad | .6 | .770 | .308 | .524 | .616 | .513 | .610 | .384 | .592 | .983 |
| myopia | quad | .7 | .909 | .576 | .804 | .836 | .787 | .846 | .670 | .774 | .990 |
| Pima | inter | .4 | .127 | .184 | .334 | .240 | .215 | .368 | .176 | .298 | .131 |
| Pima | inter | .6 | .276 | .336 | .580 | .508 | .411 | .598 | .342 | .536 | .281 |
| Pima | inter | .7 | .392 | .426 | .688 | .600 | .533 | .710 | .426 | .660 | .355 |
| Pima | quad | .2 | .070 | .170 | .270 | .192 | .201 | .282 | .156 | .204 | .080 |
| Pima | quad | .3 | .101 | .228 | .462 | .330 | .298 | .466 | .272 | .322 | .097 |
| Pima | quad | .5 | .336 | .520 | .808 | .720 | .668 | .818 | .608 | .566 | .200 |

**DeepGOF-1 minus the strongest arm, niche design (matched level, paired SE)**

| cell | strongest other arm | difference |
|---|---|---|
| myopia inter 1.4 / 2.2 / 3.0 | projection | -.840 (.016) / -.558 (.022) / -.251 (.019) |
| myopia quad .5 / .6 / .7 | projection | -.398 (.024) / -.213 (.019) / -.081 (.013) |
| Pima inter .4 / .6 / .7 | GiViTI | -.241 (.026) / -.322 (.028) / -.318 (.028) |
| Pima quad .2 / .3 / .5 | GiViTI | -.212 (.021) / -.365 (.024) / -.482 (.025) |

Against the other arms: on myopia, DeepGOF-1 is above HL, both Stukel forms, Stute-Zhu, GiViTI and mHL by
more than 2 SE in the quadratic cells and at inter C = 2.2 and 3.0 (for example quad C = .5: +.17 (.03) over
Stukel joint, +.19 (.03) over GiViTI), but the projection test is far above it everywhere, and Osius-Rojek
(not calibrated on this design) is above it in the interaction cells. On Pima, DeepGOF-1 is below both
Stukel forms, Stute-Zhu, GiViTI and Osius-Rojek by more than 2 SE in every cell, below HL and mHL by more
than 2 SE in 5 of 6 cells, and level with the projection test (which it exceeds only at quad C = .5,
+.136 (.025)).

**Dense control.** No advantage was expected. Myopia dense: DeepGOF-1 is above all arms except projection
and Osius-Rojek (projection 1.000 on the interaction). Pima dense: DeepGOF-1 is below Stukel SstBoth and
GiViTI by more than 2 SE in every cell.

**Why the niche did not carry over (diagnostics, `result_tables.md`).**
- Myopia: the two strongest continuous effects are weak on the standardised scale (|b| about .35, event
  rate 13 %), and binary covariates (DADMY .50, MOMMY .39) rank above them in the real fit. DeepGOF-1
  selected the pair carrying the departure in only 60 % of null datasets, 36-38 % under the quadratic and
  2 % under the interaction. SPORTHR is right-skewed (max 4.15 sd), so the departure sits on a few
  children with extreme values; the projection test, which weighs individual directions, picks this up
  almost always.
- Pima: DeepGOF-1 selected the right pair in 97-99 % of datasets, yet had the lowest power. The departure
  in glu (strongest effect, b = 1.1) acts mainly as a curvature in the risk scale that the calibration
  tests (Stukel, GiViTI) read directly; this repeats the GLOW S2 finding. With seven covariates, all
  continuous and none strongly collinear, the eight-irrelevant-covariate dilution that drives the
  simulated niche is absent.
- So the simulated niche depends on conditions (ten independent Gaussian covariates; strong, symmetric
  effects on the two active ones; many irrelevant covariates) that neither real design reproduces.

## Deviations from the pre-declaration (all reported)

1. Calibration rule gap: for myopia dense inter, HL never reached .10 on the grid, and the rule did not
   say what C1 is then. The code set C1 = 3.0 (grid maximum) and the declared +0.1 rule gave C2 = 3.1,
   C3 = 3.2. This affects the dense control only.
2. `calibrate.R` was edited twice after the pilot rows were written, only to fix the reporting lookup of
   HL power at C2/C3 (a floating-point match that failed, then a C value off the grid). The pilot rows and
   the rule were not changed; the final file was hashed in `RUN.sha256` before Part B.
3. Part B ran with 1, then 3, then 12 workers (stopped and resumed twice; partial chunks deleted; each
   replicate's seed is fixed, so the rows do not depend on this).
4. Not a deviation but worth noting: on the myopia niche design Osius-Rojek is not calibrated (.117), and
   mHL is not calibrated on both Pima designs; per the rule they count in "below none" but not in "above
   every calibrated arm". Removing them changes nothing: the projection test (calibrated) and Stukel /
   GiViTI (calibrated) are above DeepGOF-1 in all 12 cells.

## What this means

- The pre-declared real-data confirmation failed, clearly: on both real designs another published test
  had higher power than DeepGOF-1 in every niche cell. The paper should not claim that the niche has been
  seen on real data. No `INSERT_MAIN.tex` was written.
- What can be said honestly, if anything: on the myopia design DeepGOF-1 had higher power than HL, Stukel,
  Stute-Zhu, GiViTI and mHL for the quadratic departure, but the projection test had higher power still.
  That is not a niche result and is not recommended for the paper.
- Part A adds a practical point for the limitations section: on real data with binary or count covariates
  on the selected axes (myopia, Mroz), the DeepGOF-1 p-value depends strongly on row order until the
  tie-breaking fix is released; Mroz, where every other test except Osius-Rojek rejects, shows that
  DeepGOF-1 can miss a misfit the other tests agree on.

## Runtime

Part A: about 5 min (main arms and 60 extra DeepGOF-1 runs) plus 290-325 s per BAGofT run.
Calibration pilot (48,000 HL fits, 1 worker): about 3 min. Part B: 16,000 datasets, 13.5 CPU hours
(2.4-3.6 s per dataset), wall 05:06-06:24 EDT (3 workers until the BAGofT run finished at 05:17, then 12).
