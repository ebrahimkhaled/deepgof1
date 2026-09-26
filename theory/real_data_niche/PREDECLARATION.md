# Real-data check of the DeepGOF-1 niche: PRE-DECLARATION

Declared 2026-09-25, before any goodness-of-fit test was computed on any candidate data set. The file
is hashed (sha256, with a timestamp) into `PREDECLARATION.sha256` together with the metadata scan
(`scan_candidates.R`, `candidate_scan.csv`) and the separation gate (`check_separation.R`,
`check_separation.log`). Nothing below may be changed after hashing; any later deviation is reported
as a deviation in RESULT.md.

What was computed before hashing (and nothing else):
1. `scan_candidates.R`: dimensions and variable types of every data object in the six named packages
   (201 objects, `candidate_scan.csv`). No model, no outcome.
2. `check_separation.R`: for the data sets that pass the metadata screen and have a natural binary
   outcome, ONE glm fit of the main-effects model, printing only convergence, the number of extreme
   fitted risks, the largest |b_j| sd_j and the event count (`check_separation.log`); plus two
   design-only facts (myopia AL against ACD + LT + VCD; whether the file is sorted by the outcome).
   No coefficient table was printed and no goodness-of-fit statistic was computed.

Context. The simulated niche (`theory/niche/`, confirmed) has ten continuous covariates, two active
with an omitted interaction or quadratic and eight irrelevant; DeepGOF-1 selects the two covariates with
the largest |b_j| sd_j and reads a 6 x 6 rank-cell residual map. The projection test's efficiency falls
like 1/p when irrelevant covariates are added (`theory/local_power/irrelevant/`). The GLOW trial
(`theory/real_data_glow/`) failed because GLOW is not niche-type (6 covariates, all relevant, 4 binary).
This study asks whether the simulated advantage survives on REAL covariate designs.

## Known defect, neutralised

`deepgof1()` bins each axis by `rank(v, ties.method = "first")`. With tied covariate values the bins
depend on row order, and a file sorted by the outcome lets the map see y (GLOW, Section 4 of its
RESULT.md). A fix is being written elsewhere. Until then, in THIS study:
- Every real data set is randomly permuted ONCE before any fit: `set.seed(20260925 + k); o <- sample(n)`,
  with k = 1 (myopia), 2 (Pima), 3 (Mroz). The permutation is recorded in `perm_<data>.csv`. The
  same permuted design is used in Part A and Part B. (Mroz is stored sorted by the outcome:
  `check_separation.log`.)
- Part A reports, as a declared SECONDARY, the DeepGOF-1 p-value over 20 further random row orders
  (seeds 20261000 + 100 k + j, j = 1..20; median and range), so the order sensitivity is visible.
- In Part B the outcome is generated afresh in every replicate, independently of the (fixed, permuted)
  row order, so tie-breaking by row order cannot carry information about y.

## (a) Inclusion criteria and the complete candidate list

Criteria (all must hold):
1. A natural binary outcome of the study: one row per independent subject; not a censoring indicator of
   a survival endpoint, not a covariate, not grouped/binomial counts; no matching, clustering,
   repeated measures, case-cohort or complex-survey design.
2. n >= 200 after complete cases on the model's variables.
3. The natural main-effects model (documentation or a standard reference) has >= 6 covariates, of which
   >= 3 are continuous. "Continuous" = numeric with more than 10 distinct values.
4. No separation: glm(maxit = 100) converges, no fitted risk below 1e-8 or above 1 - 1e-8, and the
   largest |b_j| sd_j is at most 15.
5. Among several versions of the same study, the version with the largest complete-case n that meets
   1-4; a version with documented invalid values is excluded.

Screen. 201 data objects were scanned (`candidate_scan.csv`). A coarse metadata screen (data frame,
complete-case n >= 200, at least one two-level variable, at least 3 continuous variables, at least 7
columns) leaves 40. The other 161 fail it (too small, no two-level variable, fewer than 3 continuous
variables, or not a data frame) and are listed with their metadata in `candidate_scan.csv`. The 40:

| package | data | decision | reason |
|---|---|---|---|
| aplore3 | aps | exclude | outcome (placement) has 4 levels; the binary columns are covariates (multinomial example, HL Ch. 8) |
| aplore3 | burn1000 | exclude | outcome death; the model has 2 continuous covariates (AGE, TBSA) |
| aplore3 | burn_eval_1, burn_eval_2 | exclude | subsets of burn1000 (same reason) |
| aplore3 | burn13m | exclude | 1-3 matched case-control |
| aplore3 | glow11m | exclude | 1-1 matched |
| aplore3 | glow500 | exclude | the book's main-effects model (HL Table 4.9) has 2 continuous covariates; already studied (`real_data_glow/`) |
| aplore3 | glow_bonemed, glow_rand | exclude | the same GLOW cohort as glow500 |
| aplore3 | glow_mis_comp | exclude | missing-data teaching version of GLOW |
| aplore3 | icu | exclude | n = 200 and 3 continuous columns, but the only documented model (HL Exercise 2.2: AGE, CAN, CPR, INF, RACE) has 1 continuous covariate; no documented model with >= 3 continuous covariates (all 19 covariates give 40 events for 20 columns) |
| aplore3 | **myopia** | **INCLUDE** | n = 618, 81 events; 12 covariates (HL Exercise 4.6 set minus AL, see (b)), 8 continuous; no separation |
| aplore3 | nhanes | exclude | complex survey design (PSU, strata, weights; HL Sec. 6.4) |
| MASS | Aids2 | exclude | status is the censoring indicator of a survival endpoint |
| MASS | Boston | exclude | the binary column (chas) is a covariate |
| MASS | Melanoma | exclude | survival data |
| MASS | OME | exclude | grouped binomial counts, repeated measures on children |
| MASS | **Pima.tr + Pima.te** | **INCLUDE** (combined, n = 532) | the two are a random split of the same 532 complete records (MASS documentation); 177 events; 7 covariates, all continuous; no separation |
| MASS | Pima.tr, Pima.te, Pima.tr2 separately | exclude | versions of the same study (criterion 5) |
| MASS | crabs | exclude | separation (species: 200/200 fitted risks extreme; sex: 109/200) |
| MASS | epil | exclude | count outcome, repeated measures |
| mlbench | BostonHousing, BostonHousing2 | exclude | the binary column (chas) is a covariate |
| mlbench | Ionosphere | exclude | separation (43 extreme fitted risks) |
| mlbench | PimaIndiansDiabetes | exclude | same study as Pima; zeros code missing values (documented invalid values) |
| mlbench | PimaIndiansDiabetes2 | exclude | same study as Pima; complete-case n = 392 < 532 (criterion 5) |
| mlbench | Sonar | exclude | separation (208/208 extreme) |
| carData | GSSvocab | exclude | no binary outcome (vocabulary score) |
| carData | KosteckiDillon | exclude | repeated measures |
| carData | **Mroz** | **INCLUDE** | n = 753, 428 events; 7 covariates, 3 continuous (age, lwg, inc); no separation |
| carData | Rossi | exclude | survival (arrest by week) |
| carData | Wong | exclude | no binary outcome; repeated measures |
| survival | cgd, rhDNase | exclude | recurrent events |
| survival | diabetic, retinopathy | exclude | paired eyes, survival |
| survival | flchain, pbc | exclude | survival endpoints |
| survival | nwtco | exclude | case-cohort design; 1 continuous covariate |

rms ships no data sets. Included: **myopia, Pima (n = 532), Mroz**.

## (b) Models

All covariates enter linearly; binary factors as 0/1 indicators. No enriched model is used: for none of
the three data sets could a published interaction or nonlinear term be verified from a source available
here (HL's myopia interaction, Figure 10.14, involves SPHEQ, which HL exclude from the model because it
defines the outcome; Exercise 4.6 asks the reader to find interactions but gives no answer).

- **myopia** (Hosmer, Lemeshow and Sturdivant 2013, Applied Logistic Regression, 3rd ed., Wiley;
  Sec. 1.6.6, Table 1.10, Exercise 4.6; data aplore3 0.9). Outcome MYOPIC. Covariates = the Exercise 4.6
  set ("all others in Table 1.10 except SPHEQ ..., DIOPTERHR and STUDYYEAR") minus AL: AGE, GENDER,
  ACD, LT, VCD, SPORTHR, READHR, COMPHR, STUDYHR, TVHR, MOMMY, DADMY (12). AL is dropped for the same
  reason the book drops DIOPTERHR: it is a composite. AL = ACD + LT + VCD to rounding
  (correlation 1.0000, max absolute difference 0.008 mm; design-only check in `check_separation.log`).
- **Pima** (Smith et al. 1988, Proc. Symp. Comput. Appl. Med. Care, 261-265; Ripley 1996, Pattern
  Recognition and Neural Networks, CUP; Venables and Ripley 2002, MASS, 4th ed.; data MASS Pima.tr +
  Pima.te). Outcome type (diabetic). Covariates npreg, glu, bp, skin, bmi, ped, age (7): all the
  recorded predictors of the documented data set.
- **Mroz** (Mroz 1987, Econometrica 55, 765-799; Long 1997, Regression Models for Categorical and
  Limited Dependent Variables, Sage; Fox and Weisberg 2019, An R Companion to Applied Regression, 3rd
  ed., Sage; data carData). Outcome lfp. Covariates k5, k618, age, wc, hc, lwg, inc (7): the standard
  labour-force-participation logit of Long (1997) and Fox and Weisberg (2019).

Fit: glm (logit), with the Firth-type bias-reduced fallback (brglm2 AS_mean) only if glm fails to
converge (`fit_pol` of `theory/niche/run_niche.R`).

## (c) Part A: real outcomes (descriptive only)

On each included main-effects model (permuted rows), every paper test:

| arm | definition |
|---|---|
| DeepGOF-1 | `ebrahim.gof::deepgof1(fit, B = 199, K = 6)` (installed 2.8.0), set.seed(20260926 + 10 k + 1) |
| HL | `ebrahim.gof:::gof_hl` on `.gof_context(fit, G = 10)` |
| Stukel SstBoth (primary form) | `gof_stukel(ctx, list(form = "marginal"))`; the joint form is reported as secondary |
| Stute-Zhu | `gof_stutezhu(ctx)` (B = 200), set.seed(20260926 + 10 k + 2) |
| GiViTI | `givitiR::givitiCalibrationTest(o = y, e = clamped fitted risk, devel = "internal")` (the call in `gof_giviti`) |
| mHL | `gof_ftest(ctx)` (deviance residual on decile group, ANOVA F) |
| Osius-Rojek | `gof_osius(ctx)` |
| projection | `theory/addon_arms/proj_liu.R`, Liu et al. (2024) weight on RAW covariates (the model matrix without intercept), `proj_p_multi` with B = 999, set.seed(20260926 + 10 k + 3); the standardised-X column is secondary |
| BAGofT (optional) | `theory/addon_arms/bagoft_fast_pre.R` at the paper arm's settings, set.seed(20260926 + 10 k + 4); headline p.value |

All datasets are reported. The truth of a real model is unknown: Part A claims only which tests reject
at .05, and nothing about which test is right.

## (d) Part B: real-design simulation (CONFIRMATORY)

**Designs.** The two included data sets with the most covariates: myopia (12). Pima and Mroz tie at 7;
the tie is broken by the number of continuous covariates (Pima 7, Mroz 3), because the simulated niche
has continuous covariates. So Part B uses **myopia and Pima**; Mroz is Part A only.

For each of the two, X = the real covariate matrix of the main-effects model (permuted rows as above;
complete cases), with every CONTINUOUS covariate standardised to mean 0 and sd 1 (sample sd); binary
and other non-continuous covariates as recorded. X is fixed across replicates.

Fit the main-effects model to the REAL outcome on this X. Let a and b be the two CONTINUOUS covariates
with the largest |b_j| sd_j (a the larger). Two outcome designs:
- **niche**: eta0 = c0 + b_a x_a + b_b x_b, with b_a, b_b the fitted coefficients and every other
  coefficient 0; c0 solved so that mean(plogis(eta0)) equals the observed event rate.
- **dense (control)**: eta0 = the full fitted linear predictor (all fitted coefficients).

Alternatives: eta = eta0 + C g, with g = x_a x_b - mean(x_a x_b) ("inter") or g = x_a^2 - mean(x_a^2)
("quad"). y ~ Bernoulli(plogis(eta)). The fitted model in EVERY replicate is the main-effects model on
all covariates (12 or 7), with the `fit_pol` fallback.

**Calibration of C (fixed rule, pilot arm HL, not DeepGOF-1).** For each (data set, design, family):
grid C = 0.1, 0.2, ..., 3.0; 200 pilot replicates per grid value, seeds 91000000 + 100000 pid +
1000 gi + r (pid = 1..8 the (data, design, family) index in the order myopia/Pima x niche/dense x
inter/quad; gi the grid index); HL raw rejection at .05 on the fitted main-effects model.
C1 = the smallest grid C whose pilot power, and that of the next grid value, are both >= .10;
C3 = the same with >= .40 (C3 = 3.0 if never reached); C2 = the grid value nearest (C1 + C3)/2 (lower
on a tie). If C2 <= C1 then C2 = C1 + 0.1; if C3 <= C2 then C3 = C2 + 0.1. (In the simulated niche, HL
.10 to .40 corresponded to best-test power of about .3 to .8.) The pilot rows are written to
`pilot/` and are not used in any other estimate.

**Replicates.** 1,000 nulls per (data set, design); 500 per alternative cell (2 data sets x 2 designs
x 2 families x 3 C = 24 cells). Total 16,000 datasets. Seed of replicate r in scenario id: 31000000 +
100000 id + r; the scenario table (ids) is written to `scenarios.csv` before the run.

**Arms (all on the same dataset).** HL, Stukel marginal (SstBoth) and Stukel joint, Stute-Zhu (B = 200),
GiViTI, mHL, Osius-Rojek: the Part A calls. Projection: Liu weight on raw X (as given; the continuous
columns are already standardised), `proj_p_multi` with **B = 199** (B = 999 is Part A only; cost), the
weight matrix computed once per design since X is fixed. DeepGOF-1: the selected-axes reading
(two largest |b_j| sd_j, 6 x 6 rank-cell map, the shipped network weights), parametric bootstrap
B = 199, p = (1 + #{T* >= T})/(B + 1), failed refit scored +Inf, refits with `fit_pol` (the niche
study's protocol; the fast batched forward pass of `run_niche.R`, checked against
`ebrahim.gof::deepgof1` on smoke datasets before the run). BAGofT is not run in Part B (cost).
GiViTI is called in-process (no callr session) to keep memory low; a crash stops only its chunk,
which is rerun (resumable).

**NA policy.** A test that returns NA cannot reject: p = 1, identically for nulls and alternatives. NA
rates are reported.

**Estimands.** Raw size at .05 on the 1,000 nulls (binomial SE .0069); an arm is "calibrated" on a
design if its raw size is <= .05 + 3 SE = .0707. Matched-level power: each arm randomized to exactly
.05 on its own 1,000 nulls of the same (data set, design), by `matched_phi()` of
`theory/niche/analyze_niche.R` (copied verbatim). Paired difference DeepGOF-1 minus arm per cell:
mean(phi_DG - phi_arm), paired SE = sd(phi_DG - phi_arm)/sqrt(500), conditional on the null-estimated
thresholds. Raw power reported as sensitivity. Diagnostics: rate at which DeepGOF-1 selects the pair
(a, b); Firth rate; events.

**Pre-declared reading (primary: niche designs, matched level).** The 12 niche alternative cells
(myopia and Pima x inter and quad x C1..C3). Comparison arms: HL, Stukel marginal, Stukel joint,
Stute-Zhu, GiViTI, mHL, Osius-Rojek, projection.
- **NICHE CONFIRMED on real designs** if (i) DeepGOF-1 is calibrated on both niche nulls, (ii) in at
  least 6 of the 12 cells DeepGOF-1 exceeds EVERY calibrated comparison arm by more than 2 paired SE,
  and (iii) in none of the 12 cells is ANY comparison arm (calibrated or not) above DeepGOF-1 by more
  than 2 paired SE.
- **NOT CONFIRMED** otherwise; labelled **NEGATIVE** if in at least 6 of the 12 cells some comparison
  arm is above DeepGOF-1 by more than 2 paired SE.
- The same computation is reported per data set (6 cells each; rule "at least 3 of 6"), at raw level
  (sensitivity), and for the dense control. None of these can replace the primary verdict, and no
  advantage is claimed or expected in the dense control.
Everything is reported, including failures.

## Output

`partA_pvalues.csv`, `partA_axes.csv`, `partA_perm_sensitivity.csv`, `perm_<data>.csv`, `scenarios.csv`,
`calibration.csv`, `pilot/`, `rows/` (one CSV per scenario x 20 replicates, written .part then renamed;
resumable), `summary_size.csv`, `summary_power.csv`, `summary_diff.csv`, `result_tables.md`, RESULT.md.
CPU: at most 3 R worker processes until `theory/addon_arms/bagoft_waiter.log` shows the BAGofT run and
its analysis finished; then up to 12.
