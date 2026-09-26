# Tie-breaking fix for deepgof1(): RESULT

Date: 2026-09-25. Trigger: `theory/real_data_glow/RESULT.md` section 4 (GLOW M2 rejected at p = .005 only
because `glow500` is sorted by the outcome and `.dg_map()` broke ties by row order).

Files: `common.R` (old and new statistic), `fdic_roworder.R`, `glow_roworder.R`, `sim_liu3.R`,
`test_discriminates.R`, `summarise.R` -> `summarise.log`. Per-unit CSVs: `fdic_roworder.csv`,
`glow_roworder.csv`, `liu3_{null,power,power2,sorted}.csv`. Tie/order diagnostic of the FDIC file:
`../fdic_roworder/diag_ties_order.R` -> `diag_ties_order.log`.
(`liu3_null.log` and `liu3_power.log` end with a parse error printed AFTER "COMPLETE": the script file
was edited while those processes were finishing. All 1,000 rows of each were written before it.)

## 1. FDIC (paper Section "Application: bank failures after 2008"): nothing changes

- Row order of `data/fdic_2008_banks.csv` is the API order (not sorted by failure, cert or date).
  Row index vs failure within cells: Wilcoxon p .03 (GA), .24 (FL), .34 (IL), .06 (CA), .03 (WA), .84 (MN),
  i.e. weak association in two cells.
- Ties: `eq_ta` and `roa` have none. `npa_ta`, `bro_dep`, `cre_ta` are tied only at zero
  (bro_dep 25-59 % zeros; npa_ta 3-20 %; cre_ta 0-4 %).
- Axes actually selected: eq_ta+cre_ta (GA, FL, IL, WA), eq_ta+npa_ta (CA), eq_ta+roa (MN). Tied rows on
  selected axes: GA 6, FL 4, IL 7, WA 0, MN 0; CA 63 (npa_ta = 0, all 63 non-failures, so order inside the
  tie cannot carry y).
- Old statistic, stored order, paper seed: reproduces all six published p-values exactly
  (.320, .185, .805, .710, .350, .905). The observed statistic S is identical across 30 random orders
  in five cells; in CA it moves between -1.771 and -1.734.
- p over 30 random row orders (old / new): GA .305 / .307 median (range .225-.380), FL .21 / .22
  (.165-.29), IL .80 / .80 (.74-.87), CA .725 / .71 (.64-.765), WA .315 / .333 (.245-.40),
  MN .877 / .865 (.815-.925). No p <= .05 anywhere. The spread is bootstrap Monte Carlo noise (B = 199).
- New statistic, stored order, paper seed: GA .350, FL .175, IL .785, CA .655, WA .330, MN .905.
- The size column (2,000 nulls on the real design) is unaffected in law: the simulated outcomes are
  independent of row order, so row-order ties are then a y-independent tie-break.
- Conclusion: every FDIC number and the sentence "DeepGOF-1 rejects no cell (p-values .185 to .905)"
  stand. If the table is regenerated with 2.8.0 the DG-1 p-values move by bootstrap noise only.

## 2. Fix (ebrahim.gof 2.8.0, worktree feature-proj-bagoft)

`.dg_map(fit, K, tb)` ranks each axis by `order(v, tb)`. `deepgof1()` draws `tb <- sample.int(n)` once,
only if some design column has a tie, and passes the same `tb` to the observed map and to all B bootstrap
maps. The bootstrap already used the same design matrix in the same row order, so with row-order ties
the bootstrap binning was identical to the observed one; the defect was only that the tie-break was
correlated with the observed y and not with y*. With no ties nothing is drawn and results are
bit-identical to 2.7.0 (unit test checks the RNG stream).

## 3. Verification

Liu Setting 3 (x2 ~ Bin(.5)), n = 100, B = 199, paired (same data, same seed), rejection at .05:

| condition | R | old | new | new - old (paired SE) | discordant old-only / new-only | McNemar p |
|---|---|---|---|---|---|---|
| null, random row order | 1000 | .050 (.007) | .041 (.006) | -.009 (.008) | 33 / 24 | .29 |
| C = 0.5, random order | 1000 | .096 (.009) | .099 (.009) | +.003 (.011) | 56 / 59 | .85 |
| C = 1.0, random order | 600 | .300 (.019) | .302 (.019) | +.002 (.020) | 74 / 75 | 1.00 |
| null, rows SORTED by y | 500 | .282 (.020) | .046 (.009) | -.236 (.021) | 129 / 11 | < .001 |

KS on the p-value distributions (old vs new): .86, .72, .89 in random order. The paired SE of about
.008-.02 bounds a size or power change to roughly +-.02 (2 SE); the equal-in-law argument is exact.

GLOW (`aplore3::glow500`), new statistic, B = 199:
- M2 file order: p = .515 (was .005). 30 random orders: median .438, range .080-.965, none <= .05.
- M1 file order: p = .860; 30 random orders median .827, none <= .05.
- In law, M2: file order with 300 seeds vs 300 random orders with one seed: rejection .040 (SE .011)
  vs .047 (SE .012), medians .455 vs .425, KS p .25. M1 (100 each): KS p .81, no rejections.

## 4. Package checks

- testthat, full suite with NOT_CRAN = true: 82 test blocks, 0 failures, 0 errors.
- The new unit test fails on the old statistic (all 8 seeds p = 1/20; Wilcoxon sorted vs shuffled
  p = .0001) and passes on the new one (mean p .44; Wilcoxon p .51) (`test_discriminates.R`).
- `R CMD check --as-cran --no-manual` on the built 2.8.0 tarball: Status 1 NOTE ("unable to verify
  current time", offline clock check only).
