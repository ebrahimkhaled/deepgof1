# DeepGOF-1: reproducibility archive, version 2.4

Archive for the paper

> Ebrahim, E. K., Hussein, O. A., and El-Kotory, A. *Auditing a Frozen Neural Goodness-of-Fit Test for
> Logistic Regression: Validity, Certified Consistency and Local Power of DeepGOF-1.*

and its Supporting Information; the paper is now titled *Where Does a Logistic Regression Model Fail? A Pretrained Neural Goodness-of-Fit Test with Bootstrap Validity*, and its supplement is Online Resources 1 (proofs) and 2 (supporting studies). DeepGOF-1 is a goodness-of-fit test for logistic regression
whose statistic is a small convolutional network (18,273 frozen parameters) that reads the fitted
model's *residual map*, a 6 x 6 grid of standardized residual sums over the ranks of the two strongest
covariates (or, in the all-pairs reading, of every pair of covariates, taking the largest score). The p-value is the rank of the observed score inside the analyst's own parametric
bootstrap (B = 199), so the level is a property of the calibration, not of what the network learned.

| DOI | what it identifies |
|---|---|
| 10.5281/zenodo.22113219 | concept DOI: always the newest version |
| 10.5281/zenodo.22113220 | version 1.0.0, the first archive |
| (shown on this record's page) | this version, 2.4.0: cite it to reproduce a number of the paper above |
| 10.32614/CRAN.package.ebrahim.gof | the R package that distributes the deployable test as `deepgof1()` |

## What changed from version 2.3.0

* **The pair score test** (`theory/tier3/`): a score test for the squares and the product of the pair DeepGOF-1
  selects, calibrated by the same bootstrap with the pair selected again in every sample; it is `pairscore.gof()`
  in `ebrahim.gof` 3.0.0. On the Liu et al. (2024) benchmark, Settings 2-4 (43,000 data sets with the design with
  ten covariates), together with the chi-square statistic of the map, the GRP test (GRPtests) and Stukel's test with
  the one-sided fallback. `20_run_AB.R` reproduces the stored DeepGOF-1 p-values exactly on every data set.
* **Stukel and Osius-Rojek on all 144,000 benchmark data sets** (`theory/tier3/30_run_C.R`).
* **Networks trained on reduced corpora** (`theory/tier3/10_corpus.R`, `11_train_export.py`, `nets/`): without
  interaction terms, without quadratic terms, and a stricter variant, three seeds each, with three full-corpus
  controls; every exported network is verified against the training framework (`NET_VERIFY.csv`).
* **Departures of other shapes** (`theory/tier3b/`): threshold, local bump, oscillation and corner, in the design
  of Table 2 and the design with ten covariates (32,000 data sets), and the limiting power of the pair score test
  for the departures of Table 2 (`50_limit_score.py`). `REPORT.md` summarises the study.
* `tables/` in each folder holds every number of the paper's Tables 2, 3 and 5 and of Online Resource 2, Section 29;
  `40_analyze.R` recomputes them from the per-data-set rows.

## What changed from version 2.2.0

* **The network against simple statistics of its own map, and the combined reading** (`theory/ablation/`).
* **Locating the misfit**: the map against Lin-Wei-Ying cumulative residuals and a spline screen (`theory/localization/`).
* **External validation of a published risk model**: exact Monte Carlo tests of frozen predictions (`theory/external/`);
  the test is `deepgof1.external()` in `ebrahim.gof` 2.9.0.
* **SUPPORT**: the ladder read by the final 2.9.0 readings, a plasmode on the study's own covariates up to all
  8,873 patients, the computing times of DeepGOF-1 and the projection test, and a split-sample external validation
  over 20 random splits (`application/support/14`-`18`).
* **Corrupted records at n = 500 and 1,000**: DeepGOF-1, the projection test and BAGofT on the stored datasets of an
  earlier corrupted-record study, beside that study's stored p-values (`theory/corruption/`).
* Sections 4 and 5 below follow the paper's and the Supporting Information's current numbering.

## What changed from version 2.1.0

* **Very small samples and wrong links** (`theory/smalln_link/`): the size of five tests at n = 20 to 100, their matched power against an omitted quadratic or interaction at n = 20 to 50, and against four wrong links (probit, complementary log-log and two Stukel links) at n = 20 to 100; the Supporting Information tables and the main-text figure are generated from these files.
* Figure scripts brought up to date: `figures/fig_smalln_size.py` (new), `figures/fig_support_maps.py` (the three-panel SUPPORT figure), `figures/fig_are_irrelevant_sc.py`.

## What changed from version 2.0.0

* **Bug fix in the deployed test.** Up to archive 2.0.0 (and `ebrahim.gof` 2.8.0), `deepgof1()` refitted
  each bootstrap sample by evaluating the model formula on the model frame. For a term that transforms a
  covariate (`log(x)`, `ns(x, 3)`, `poly(x, 2)`) the model frame holds the transformed column and not
  `x`, so every refit failed, every replicate scored `+Inf` and the p-value was 1. `R/deepgof1.R` now
  refits on the fitted model's design matrix, as `ebrahim.gof` 2.9.0 does. For models without such terms
  the p-value for a given seed is unchanged. No simulation of the paper fits a transformed term; in the
  SUPPORT application it affected the spline repairs, whose p-values are recomputed here.
* **The all-pairs reading.** `deepgof1(fit, reading = "allpairs")` takes the largest score over the maps
  of every pair of covariates, the same maximum in every bootstrap replicate. The default reading is
  unchanged. With one covariate both readings use the 36-cell map along it. `R/check_against_package.R`
  checks that `R/deepgof1.R` and the package give the same p-value for the same seed.
* **New studies**: the shipped network on the level-study nulls and on the benchmark, with both readings
  (`theory/shipped_network/`, `theory/allpairs_study/`); a U-shaped covariate among linear ones
  (`theory/ushape/`); the projection test on the corrupted-record datasets
  (`theory/small_n/proj_addon/`); and the SUPPORT application (`application/support/`).

## What changed from version 1.0.0

* The archive is organised around the paper above. Every table, figure and supplementary result maps
  to files below (Sections 4 and 5 of this README).
* **Ties.** `R/deepgof1.R` breaks ties among covariate values at random, once per call: a random
  permutation of the rows is drawn only when some covariate column has a tied value, and the same
  permutation serves the observed map and every bootstrap map. This is the rule of `ebrahim.gof` 2.8.0
  and of the paper's Algorithm 1. With continuous covariates nothing is drawn and the p-value is
  bit-identical to version 1.0.0 (checked against the package on untied, binary and rounded designs:
  identical statistic, bootstrap scores and p-value for the same seed). `ties = "roworder"` restores
  the earlier rule and exists only to reproduce numbers computed before the fix (see Section 6).
* The benchmark harnesses are self-contained (`library(ebrahim.gof)`; no other project needed) and
  compute the paper's arms. A smoke test regenerates datasets from their seeds and reproduces the
  deposited p-values exactly.
* The per-replicate files carry the paper's arms: Hosmer-Lemeshow (HL), HL equal-width,
  Pigeon-Heyse, Stukel (the SstBoth form of Liu et al. 2024), Osius-Rojek, Stute-Zhu, mHL, GiViTI, the
  projection test of Liu et al. (2024), BAGofT and DeepGOF-1.
* New: `theory/` (the paper's theoretical results and their numerical checks), `tables/` (one script
  regenerates every generated number of the tables), `bagoft/`, and the FDIC application in
  `application/fdic/`.

## 1. Quick start

Pure base R for the test itself; no Python and no GPU.

```r
source("R/deepgof1.R")                      # from the archive root
fit <- glm(y ~ x1 + x2, data = d, family = binomial())
set.seed(1); deepgof1(fit)                  # B = 199; report the seed when an axis has ties
```

or `install.packages("ebrahim.gof"); ebrahim.gof::deepgof1(fit)`. Use version 2.9.0 or later: it refits
the bootstrap on the design matrix and offers `reading = "allpairs"`; 2.8.0 and later break ties at
random, earlier versions by row order.

| check | command (from the archive root) | expected |
|---|---|---|
| deployed test end to end | `Rscript smoke_test.R` | misspecified model p = .005; correct model p = .915 |
| base R against the training framework | `Rscript training/VERIFY.R` | max difference 4.6e-11, PASS |
| benchmark pipeline on tiny blocks | `sh benchmark/smoke_test.sh` | every regenerated p-value equals the deposit (difference 0) |
| every generated number of the tables | `python tables/make_tables.py` | `tables/tables_generated.tex`, identical to the manuscript's |
| archive copy of the test against the package | `Rscript R/check_against_package.R` | the same p-value for every model and reading |
| FDIC adequacy calls, both tie rules | `Rscript application/fdic/fdic_tie_check.R` | fixed (random) rule identical to the deposit and the paper; row-order rule gives the pre-fix values |

Software used: R 4.4.1 with `ebrahim.gof` 2.8.0 and 2.9.0, `brglm2`, `givitiR` 1.3 and `BAGofT` 1.0.0; Python 3
with numpy, pandas, scipy, scikit-learn and PyTorch (training and the per-block networks only).

## 2. Layout

| folder | contents |
|---|---|
| `R/` | `deepgof1.R` (the test) and `deepgof1_forward.R` (the forward pass in base R) |
| `inst/extdata/` | frozen weights, standardizer, and the cross-language verification maps |
| `training/` | `01_corpus.R` (corpus generator), `ship_corpus.csv` (the 8,400 labelled training maps), `02_train.py` (network and training loop), `03_export_to_R.py` (export; documents the row-major/column-major port trap), `VERIFY.R` |
| `benchmark/` | `liu/` (Simulation II harness), `grid/` (Simulation I harness), `common/` (matched-level device, network), `theory/` (certificate and blind-cone checks), `smoke_test.sh` |
| `results/` | `liu/` (24 blocks), `s3/` (60 blocks), `smalln/`, `fdic/`, `theory_blind_direction.csv`. The pooled summaries of the first run (`results/liu/AUGSCORE_INTERIM.csv`, `SUMMARY_AUG.csv`, `VERDICT.csv`, `AUGMENT_SUMMARY.csv`, and each block's `summary.csv`) are kept with the rows of the paper's arms; `tables/make_tables.py` recomputes Table 3 from the per-replicate files and checks it against `AUGSCORE_INTERIM.csv`. The column `proj_superseded` of `augment_pvalues.csv` is the first run's projection statistic (a different weight); it serves only as the pairing key of the corrected projection test |
| `bagoft/` | `confirm/` (the paired comparison of Section 6.2), `explore/` (the earlier exploratory trial, the fast implementation and the repair) |
| `application/fdic/` | FDIC data (from two keyless API calls) and the three scripts of the application |
| `tables/` | `make_tables.py` (writes `tables_generated.tex` and `results/liu/TABLE2_SCORES.csv`), `stukel_na.py` |
| `figures/` | the scripts of Figures 1, 2, 4 and 5 (Figure 3 is in `theory/local_power/irrelevant/`) |
| `theory/` | one folder per theoretical question (Section 5); `_support/rdeploy/` holds the deployed forward pass, weights and scaler the theory scripts load |

### How the simulation p-values were produced

* A **block** is one (design or setting, n, corpus seed). Datasets are generated from a deterministic seed
  per replicate (written at the top of each harness), so any single replicate can be re-derived.
* The classical arms are `ebrahim.gof` implementations; Stukel is the SstBoth form of Liu et al. (2024)
  (`gof_stukel(form = "marginal")`), which is what the benchmark's Stukel column holds throughout.
* The DeepGOF-1 column applies the deployed test (map, network, rank among B = 199 bootstrap maps) with
  the network trained on the block's training corpus: half of the block's nulls, as bootstrap
  replicates, against its random training departures (`kind = 1` rows), by `benchmark/liu/liu_score.py`
  and `benchmark/grid/s3_score.py`. The corpus seeds (101000, 102000 and, in the grid, 103000) are the
  paper's independent training corpora. In the grid, one network per (n, seed) is trained on design D1
  and applied unchanged to D1-D4. The paired BAGofT comparison (`bagoft/confirm/`) reads these
  per-block DeepGOF-1 p-values of corpus 101000 (`results/liu/*_s101000/deep_pvalues.csv`), and the
  classifier (AUC) diagnostic of the paper's Section 3.1 is the `score_auc` column of the grid's
  `summary.csv` files.
* **Which network.** The benchmark simulations (Simulations I and II, with the BAGofT comparison)
  therefore use per-block networks with the architecture and training loop of the shipped network, not
  the shipped weights. Every other DeepGOF-1 result of the paper (Sections 3-4 and 7-9: theory checks,
  local power, certificates, the simulation with two active covariates among ten and its real-design
  check, the corrupted-record studies, FDIC) and of the Supporting Information uses the shipped weights, loaded
  from `inst/extdata/` by `R/deepgof1.R`, from the package, or from `theory/_support/rdeploy/weights.csv`
  (the same file).
* NA policy: in Table 3 (the Liu benchmark) each test is evaluated on the datasets where it returns a
  p-value, for nulls and alternatives alike (`tables/make_tables.py`; the original convention, no
  p-value = no rejection, is kept in the `*_na1` columns of `results/liu/TABLE2_SCORES.csv` and gives
  Stukel's test 0.643 overall instead of 0.610). Elsewhere, and in each block's `summary.csv`, a test
  that returns no p-value cannot reject (p = 1).
  Matched power: each arm is randomized to exactly .05 on the block's held-out nulls
  (`benchmark/common/matched_power.py`).

## 3. Reproducing

| step | command | time |
|---|---|---|
| one Liu block | `LU_SET=2 LU_N=100 LU_SEED=101000 Rscript benchmark/liu/liu_grid.R`, then `python benchmark/liu/liu_score.py --dir results/liu/S2_n100_s101000`; Stute-Zhu: `benchmark/liu/liu_augment.R` | CPU-hours per block; n = 500 is the slowest |
| mHL, GiViTI, projection, BAGofT | `theory/addon_arms/` (`liu_addon.R`, `proj_addon.R`, `niche_addon.R`, `niche_bagoft_fast.R`) | about 0.4 CPU-s per dataset for mHL and GiViTI; 0.5-3.5 s for the projection test |
| one grid block | `S3_DESIGN=D1 S3_N=100 S3_SEED=101000 Rscript benchmark/grid/s3_grid.R` (all four designs), `python benchmark/grid/s3_score.py --n 100 --seed 101000`, `python benchmark/grid/s3_collect.py` | CPU-hours per block |
| Table 3, Table 4, Table 8 and the text macros | `python tables/make_tables.py` | minutes |

All harnesses read the environment variable `DEEPGOF_ARCHIVE` (default: the current directory) and
run with two worker processes unless told otherwise. The deposited per-replicate files exist so that
nobody has to rerun the multi-CPU-day simulations to check a number.

**Scripts added in version 2.1.0** (`theory/shipped_network/`, `theory/allpairs_study/`, `theory/ushape/`,
`theory/small_n/proj_addon/`, `application/support/`) keep the absolute paths of the author's project
tree, which map as below. Some of them load `ebrahim.gof` 2.9.0 from a private library (`lib290`,
`lib290c`) because they ran before that version reached CRAN; with 2.9.0 or later installed, drop the
`lib.loc` argument. `application/support/` needs no path changes: run its scripts in order from that folder.

**Path constants of the theory scripts.** Scripts under `theory/` were run in place in the author's
project tree; their path constants now read `/path/to/project/paper_deepgof/...`. Map them as follows:
`.../theory/<folder>` is `theory/<folder>` of this archive; `.../code/rdeploy` is `theory/_support/rdeploy`;
`.../release/deepgof1/results` and `.../results` are `results/`; `.../confirm_bagoft` and
`.../explore_bagoft` are `bagoft/confirm` and `bagoft/explore`; `/path/to/ebrahim.gof` is a clone of the
package (version 2.8.0). The hash-covered scripts are deposited byte-identical so that they still match
their recorded hashes; their path constants therefore keep the author's original absolute project path
(ending in `.../paper_deepgof/...`), which maps in the same way. Scripts and comments that name files of the project's `code/` folder refer to
their archive copies: `liu_grid.R`, `liu_augment.R` -> `benchmark/liu/`; `s3_grid.R`, `s3_score.py` ->
`benchmark/grid/`; `universal_design.R` -> `benchmark/grid/s3_designs.R` (the same generator
functions); `map_train.py` -> `benchmark/common/`; `probit_common.py` (`randomized_matched_power`) ->
`benchmark/common/matched_power.py`; `fdic_deepgof.R` -> `application/fdic/`. The scripts in `benchmark/theory/` run with the working directory
`theory/_support/rdeploy`. `theory/_support/rdeploy/deepgof.R` is the wrapper the theory simulations
used; it breaks ties by row order, which matters only for tied covariates, and every theory simulation
uses continuous covariates. For data analysis use `R/deepgof1.R` or the package.

## 4. The paper: tables, figures and in-text results

Numbers are those of the paper's current version (twelve sections, eight tables, seven figures).

| paper item | what | files |
|---|---|---|
| Figure 1 (residual maps) | simulated maps of five misfits | `figures/fig3_maps_sc.py` (with sign marks: `fig3_maps_sc_signs.py`), style `ekml_style_sc.py` |
| Figure 2 (architecture) | the network and the deployed test | `figures/fig1_architecture_sc.py` |
| Section 2.1 (the readings: axis rule, all pairs, combined) | the deployed code | `R/`; the package `ebrahim.gof` 2.9.0 (`deepgof1()`, `deepgof1.external()`) |
| Section 3.1 (validity; gap estimate) | the tail-localized gap, 90th percentiles .030 (n = 50) and .007 (n = 200) | `theory/calibration_gap/` (`per_dataset.csv`, `cell_summary.csv`, `PILOT_RESULTS.md`) |
| Section 3.2, Table 1 (certificate margins, blind cone) | recession slopes along five directions; the most negative direction found | `benchmark/theory/theory_verify_numbers.R`, `theory_sphere_check.R`, `theory_blindcone_probe.R`, `theory_recession_check.R`; `results/theory_blind_direction.csv` |
| Section 4.1, Table 2 (local sample-size ratios) | HL and Stukel columns; projection column | `theory/local_power/FORMULA.json` (`ARE_vs_DG`), `01_local_power.py`; projection: `theory/local_power/irrelevant/PROJ_LIU_RAW.json` (k = 0); SstBoth: `theory/stukel_sstboth/TAB_ARE_SSTBOTH.json` |
| Section 4.2, Figure 3 (covariates that carry no signal) | ratio against q irrelevant covariates | `theory/local_power/irrelevant/fig_are_irrelevant.py` (reads `ARE_IRRELEVANT.json`, `PROJ_LIU_RAW.json`); SstBoth: `theory/stukel_sstboth/` |
| Section 5, Figure 4 (level at very small n) | size of five tests at n = 20 to 100 | `theory/smalln_link/` (see Supporting Information Section 18) |
| Section 5, Figure 5 (the sixty-cell level study) and the shipped network's twenty cells | per-replicate p-values of every arm and cell | `results/s3/<design>_n<n>_s<seed>/`, `results/s3/VERDICT.csv`, figure `figures/fig6_s3grid_sc_v2.py`; the shipped network: `theory/shipped_network/` (`run_shipped_level.R`, `rows/`, `SHIPPED_LEVEL.csv`) |
| Section 6, Table 3 and Figure 6 (the extended Liu benchmark) | per-replicate p-values; per-block scores | `results/liu/<block>/` (`blk_pvalues.csv`, `augment_pvalues.csv`, `deep_pvalues.csv`), `theory/addon_arms/` (mHL, GiViTI, projection), `bagoft/confirm/` (BAGofT); table `tables/make_tables.py`; figure `figures/fig8_liu_sc_v2.py` |
| Section 6.1 (the shipped network and its all-pairs reading) | both readings on the same datasets, ranking, paired standard errors | `theory/shipped_network/`, `theory/allpairs_study/` |
| Section 6.1 (the network or the map?) | the network against the sum of squared cells and the largest cell of the same map | `theory/ablation/` (see Supporting Information Section 21) |
| Section 7, Table 4 (two active covariates among ten) | one row per dataset, the decision rule, BAGofT at n = 200; the network against simple map statistics | `theory/niche_paper/` (`NICHE_ROWS.csv`, `analyze_niche_paper.py` writes `niche_generated.tex`); `theory/ablation/`; the U-shape study `theory/ushape/`; the real designs `theory/real_data_niche/` |
| Section 8 (locating the misfit) | the map, Lin-Wei-Ying and the spline screen on the same datasets | `theory/localization/` (see Supporting Information Section 22) |
| Section 9 (a published risk model on new patients) | exact Monte Carlo tests of frozen predictions | `theory/external/` (see Supporting Information Section 23) |
| Section 10, Tables 5 and 6 (one corrupted record at n <= 50; the repair rules) | per-replicate rows of Studies 1 and 2 | `theory/small_n/` (`rows/`, `rows_trim/`, `RESULTS.json`, `RESULTS_TRIM.json`, `sstboth_addon/`, `touch_rate/`, `proj_addon/`) |
| Section 10.3, Table 7 (corrupted records at n = 1,000) | DeepGOF-1, projection and BAGofT on the stored datasets of the earlier corrupted-record study, beside its stored p-values | `theory/corruption/` (see Supporting Information Section 25) |
| Section 11, Table 8, Figure 7 (the SUPPORT application) | scripts in the order they were run, the public data, every result; the 2.9.0 ladder, the plasmode, the computing times and the split-sample validation | `application/support/` (`01`-`13` as before; `14_ladder_290d.R`, `15_plasmode.R`, `16_projection_timing.R`, `17_deepgof_timing.R`, `18_external_split.R`; results in `application/support/results/`) |
| the generated numbers of Sections 6.1, 7-11 and Supporting Information Sections 21-25 | macros and tables written from the result files | `tables/make_round2_si.py` |
| the deployed test | deployed code, weights, corpus, generator, training and export code, harnesses | `R/`, `inst/extdata/`, `training/`, `benchmark/`, `bagoft/` |
| ties | the tie-breaking study: FDIC in 30 random row orders, Liu Setting 3, GLOW | `theory/tie_fix/`, `theory/fdic_roworder/` |

## 5. The Supporting Information: results and files

S-numbers and the section numbers below are those of the Supporting Information. Results that are proved
analytically have no data files; the files listed are the numerical checks and inputs.

| Supporting Information | results | files |
|---|---|---|
| Table S1 | correspondence of main-text and supplement results | (none) |
| Sections 1-3 | setting, assumptions A1-A7, Lemmas A.1-A.2 | (none) |
| Section 4 | Lemma S.1 (recession function) and Remark R1 (slope check) | `benchmark/theory/theory_recession_check.R`, `theory_verify_numbers.R` (run in `theory/_support/rdeploy`) |
| Section 5 | Lemma S.2, Proposition S.3 (validity) | size evidence: `results/s3/` |
| Section 6 | Lemmas A.3-A.5, Propositions S.7, S.8, the exact maximized Monte Carlo variant, the AUC diagnostic, the gap pilot | `theory/calibration_gap/`, `theory/exact_mmc/` |
| Section 7 | Theorem S.11 (data-computable gap bound), coverage | `theory/level_gap/` |
| Section 8 | the level error against n | `theory/level_gap/rate/` |
| Section 9 | Theorem S.4 (conditional consistency) | firth rates: `results/s3/*/summary.csv` (`firth_rate`) |
| Section 10 | Corollaries S.5, S.6 (blind cone, representation limit) | `benchmark/theory/theory_sphere_check.R`, `theory_blindcone_probe.R`, `results/theory_blind_direction.csv` |
| Section 11 | Theorem S.9, Corollary S.10 (local asymptotic power), numerical corroboration | `theory/local_power/` |
| Section 12 | covariates that carry no signal | `theory/local_power/irrelevant/` |
| Section 13 | power certificates | `theory/finite_power/` (`smoothing/`, `bab/`) |
| Section 14 | one corrupted record: theory | `theory/small_n/SPIKE_RESPONSE.json`, `03_trim_miss.py` |
| Section 15 | one corrupted record: Studies 1 and 2 | `theory/small_n/rows/`, `rows_trim/`, `RESULTS.json`, `RESULTS_TRIM.json`; SstBoth: `sstboth_addon/` |
| Section 16 | the deletion rule with more, or non-normal, covariates | `theory/small_n/touch_rate/` |
| Section 17 | Theorem S.13 (Neyman-Pearson gap) | `theory/optimality/` |
| Section 18 | very small samples and wrong links | `theory/smalln_link/` |
| Section 19 | the paired comparison with BAGofT | `bagoft/confirm/`, `bagoft/explore/` |
| Section 20 | level transfer to bank failures after 2008 | `application/fdic/`, `results/fdic/` |
| Section 21 | the network against simple statistics of the same map; the combined reading | `theory/ablation/` (`run_ablation.R`, `rows/`, `analyze_ablation.R`, `ABLATION_POWER.csv`, `ABLATION_SIZE.csv`) |
| Section 22 | locating the misfit | `theory/localization/` (`run_localization.R`, `rows/`, `analyze_localization.py`, `LOCALIZATION.csv`) |
| Section 23 | external validation of a published risk model | `theory/external/` (`run_external.R`, `rows/`, `analyze_external.py`, `EXTERNAL_LEVEL.csv`, `EXTERNAL_POWER.csv`); `deepgof1.external()` in `ebrahim.gof` 2.9.0 reproduces the DeepGOF-1 columns exactly |
| Section 24 | the SUPPORT design: level, power, computing times, split-sample validation | `application/support/` (`15_plasmode.R`, `15_analyze_plasmode.py`, `16_projection_timing.R`, `17_deepgof_timing.R`, `18_external_split.R`), `application/support/results/` (`plasmode2_rows/`, `PLASMODE2_LEVEL.csv`, `PLASMODE2_POWER.csv`, `projection_M0.csv`, `deepgof_timing_M0.csv`, `external_split.csv`, `dg290d_ladder.csv`) |
| Section 25 | corrupted records at n = 500 and 1,000 | `theory/corruption/` (`run_dg_on_edge_data.R`: DeepGOF-1 on the rebuilt datasets, each checked against its stored events and slope; `rows_edge_mt/`; `run_rivals_n1000.R`, `rows_rivals/`: projection and BAGofT; `input_battery/`: the stored p-values of the other tests; `analyze_dg_on_edge_data.py`, `CORRUPTION_FALSE_ALARMS_MT.csv`; `check_clean_level.R`: the level check on fresh seeds) |
| Section 26 | benchmark bookkeeping: missing p-values and failed refits | `tables/stukel_na.py`; failed projection refits: `proj_nfail` in `theory/addon_arms/proj_rows_*` and `NICHE_ROWS.csv` |
| Section 27 | deposit index, training protocol | `training/` |

## 6. Ties and the FDIC numbers

The paper's FDIC adequacy calls (the "DG-1 observed p" column) use the fixed rule of `R/deepgof1.R` and
of `ebrahim.gof` 2.8.0. `application/fdic/fdic_deepgof.R` defaults to it (`FDG_TIES=random`) and
reproduces the deposited `results/fdic/deepgof_observed.csv` exactly: p = .350, .175, .785, .655, .330,
.905 for GA, FL, IL, CA, WA, MN (`fdic_tie_check.log`; `FDG_REPS=0` runs these calls alone, in seconds).
`FDG_TIES=roworder` remains as an option and gives the values computed before the fix, .320, .185, .805,
.710, .350, .905: Monte Carlo movement of the B = 199 bootstrap only; no cell rejects under either rule,
and `theory/tie_fix/` reports 30 random row orders per cell under both rules. `eq_ta` and `roa` have no
ties; `npa_ta`, `bro_dep` and `cre_ta` are tied only at zero. The deposited size column (2,000 nulls per
cell on the real design; `deepgof_size.csv`, `deepgof_pvalues_<S>.csv`, the paper's "DG-1 size") was
computed with ties broken by row order; the script's size part keeps that rule by default
(`FDG_SIZE_TIES=roworder`) so that it reproduces the deposit, and `FDG_SIZE_TIES=random` reruns it under
the fixed rule. In the simulations the covariates are continuous
or their rows are in random order with respect to y, where the two rules agree in law
(`theory/tie_fix/RESULT.md`). The real-design check in `theory/real_data_niche/` permuted each data
file once at random before any fit (`perm_<data>.csv`), which is the random rule with one fixed
permutation.

## 7. Not in this archive

The project's working notes and draft analyses are not included. Every number of the paper and of its
Supporting Information can be regenerated from the files above.

## Licence

Code: MIT. Data, results, corpora and weights: CC BY 4.0. See `LICENSE`.

## Citing

> Ebrahim, E. K. (2026). *Auditing a Frozen Neural Goodness-of-Fit Test for Logistic Regression:
> Validity, Certified Consistency and Local Power of DeepGOF-1* (Version 2.2.0) [Software and data].
> Zenodo. https://doi.org/10.5281/zenodo.22113219
