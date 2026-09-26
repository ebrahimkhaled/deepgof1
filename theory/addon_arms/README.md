# Add-on arms: mHL, GiViTI, Stukel (SstBoth), the projection test, and BAGofT on the niche

These scripts add arms to the two stored simulations -- the extended Liu et al. (2024) benchmark
(`results/liu/`, 24 blocks) and the simulation with two active covariates among ten
(`theory/niche_paper/`) -- on the SAME datasets: every dataset is regenerated from its own seed and
the pairing is checked on every row.

## The arms

| arm | code | definition |
|---|---|---|
| mHL | `ebrahim.gof:::gof_ftest` (run.all.gof row "F-test") | `lm(deviance residual ~ G)`, G the ten equal-frequency groups of the fitted risk (the HL groups), then the ANOVA F test (Liu et al. 2024, Sec. 2.2) |
| GiViTI | `givitiR::givitiCalibrationTest(o = y, e = p, devel = "internal")$p.value` (givitiR 1.3), the call inside `ebrahim.gof:::gof_giviti` | one `callr` session per worker, 120 s timeout; `giviti_status` 0 ok, 1 givitiR error, 2 timeout, 3 session died |
| Stukel | `ebrahim.gof:::gof_stukel`, `form = "marginal"` (SstBoth, the form Liu et al. report; the paper's Stukel) and `form = "joint"` | both forms are carried on every row |
| projection | `proj_liu.R` (`proj_A_liu()`, compiled copy `proj_liu_c.c`; build with `R CMD SHLIB proj_liu_c.c`) | Liu et al.'s statistic after Escanciano (2006): A_ijl = (pi - angle(x_i - x_l, x_j - x_l)) / (2 pi), 1/2 if one difference is zero, 1 if both, on the RAW covariates (`proj_liu_raw`, the paper's arm); `proj_liu_std` is the same weight on standardized covariates; B = 199 parametric bootstrap with the benchmark's refit rule, a failed refit scoring -Inf |
| BAGofT (niche, n = 200, 100 replicates per scenario) | `bagoft_fast_pre.R` | BAGofT 1.0.0 settings (`testGlmBi(logit)`, `parRF(ntree = 60)`, `nsplits = 100`, `nsim = 100`, `ne = floor(5 sqrt(n))`); the fast copy includes the package's variable pre-selection step (`dcPre()`, five covariates by distance correlation) and reproduces the package `identical()`-ly on 20 datasets (`verify_fast_pre.R`, `verify_fast_pre.csv`, `verify_fast_pre/`) |

## Files

| file | what it does |
|---|---|
| `addon_common.R` | mHL, GiViTI (session), both Stukel forms, the HL pairing check, an atomic CSV writer |
| `liu_addon.R` | the Liu grid: regenerates each dataset from `SEED0 + SET*500000 + kind*50000 + ci*5000 + i` and writes `liu_rows/<block>.csv` |
| `niche_addon.R` | the niche: regenerates each dataset from `123000 + id*100000 + r`; writes `niche_rows/` (and, in `AD_MODE=bagoft`, `niche_rows_bagoft/`) |
| `proj_addon.R` | the projection arm on both simulations, with the same bootstrap outcomes as the stored run (the random-number state left by the arm computed just before it is replayed); writes `proj_rows_liu/<block>/chunk_*.csv`, `proj_rows_niche/` |
| `validate_proj_liu.R`, `check_proj_liu_c.R` | the closed-form weight against 200,000 random directions, and the compiled weight against the R definition (logs beside them) |
| `verify_pairing.R` | regenerates datasets and compares HL and Stukel with the stored rows; checks the GiViTI session against `run.all.gof(tests = "GiViTI")` |
| `analyze_addon.R`, `analyze_proj.R`, `analyze_bagoft_niche.R` | raw rejection rates with binomial SE (`*_rates.csv`), the projection pairing tables (`proj_pairing_*.csv`) |
| `niche_bagoft_fast.R` | the BAGofT niche run (refuses to start unless the stage-B verification shows 20 identical datasets) |

The niche scripts read the stored niche rows; in this archive those rows are the joined table
`theory/niche_paper/NICHE_ROWS.csv` (paper arms only), so point their input at that file.

## Pairing and the Stukel form (verify_pairing.log)

* Liu grid: HL identical to the stored rows on 210 of 210 regenerated datasets (six blocks), and on all
  6,000 rows of the one full block rerun; the stored Stukel column equals the marginal (SstBoth) form on
  210 of 210 datasets. The benchmark's Stukel column is therefore SstBoth throughout.
* Niche: HL identical on 100 of 100 datasets and on all 1,000 rows of one full scenario; the stored Stukel
  column of the niche run is the joint form, and `NICHE_ROWS.csv` carries both forms (`Stukel_marginal`
  is the one the paper reports).
* Projection: the recomputed first weight (`proj_old`) equals the stored first-run column on 90 of 90 Liu
  and 28 of 28 niche datasets in the checks, and wherever it was recomputed in the full runs.

`liu_rows_smoke/`, timing runs and superseded two-column runs are not deposited.
