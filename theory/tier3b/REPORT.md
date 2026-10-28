# Tier 3b: score test vs DeepGOF-1 on departures that are not polynomial (2026-10-24)

All numbers are in `tables/`. The scripts are `b_common.R`, `00_hmeans.R`, `20_run.R` (smoke / pilot / main), `11_choose_C.R`, `30_verify.R`, `40_analyze.R` and `50_limit_score.py`.

## Design
- **Designs.** D2: two independent N(0,1) covariates, eta = -0.5 + x1 + 0.6 x2 (the design of Table 2). D10: ten covariates, eta = -0.5 + x1 + 0.8 x2 + 0.15 sum x3..x10. In both, the departure C*h(x1, x2) is added and centred to mean zero, and the fitted model is linear in every covariate.
- **Shapes.** threshold 1{x2 > 1}; bump exp(-((x1-0.5)^2+(x2-1)^2)/0.5); oscillation sin(2 x1); corner 1{x1 > 0.8, x2 > 0.8}. The corner uses C < 0, because with C > 0 it barely moves the data.
- **Sizes.** n = 200 and 500. Three values of C per (design, shape) were chosen from a pilot at n = 316, aiming at DeepGOF-1 power of 0.2, 0.5 and 0.8 (`C_CHOSEN.csv`). In D10 the largest corner value reached only about 0.68.
- **Data sets.** 500 misspecified data sets per setting (24,000 in all) and 2,000 correct-model data sets per (design, n) (8,000).
- **Tests.** One bootstrap with B = 199 per data set gives DeepGOF-1 (distributed network, default version), the chi-square of the 36 cells and the score test (bootstrap and asymptotic). GRP, Hosmer–Lemeshow and Stukel (with the one-sided fallback) were run on every data set. No p-value was missing.
- **Checks.**
  - The D2 loop matched `deepgof1()` 2.8.0 exactly on the smoke data sets.
  - 400 regenerated data sets matched their fingerprints.
  - The RNG state was asserted unchanged after every third-party call.

## Type I error rate at 0.05 (`B_TYPE1.csv`)
| | D2 n=200 | D2 n=500 | D10 n=200 | D10 n=500 |
|---|---|---|---|---|
| DeepGOF-1 | 0.052 | 0.051 | 0.055 | 0.057 |
| chi-square of map | 0.054 | 0.046 | 0.054 | 0.055 |
| score test (bootstrap) | 0.058 | 0.046 | 0.054 | 0.046 |
| score test (asymptotic) | 0.064 | 0.046 | 0.071 | 0.052 |
| GRP | 0.006 | 0.008 | 0.001 | 0.003 |
| HL | 0.056 | 0.049 | 0.056 | 0.054 |
| Stukel | 0.067 | 0.061 | 0.080 | 0.057 |

## Size-adjusted power
**Mean over all 48 settings** (`B_MEAN_POWER.csv`): DeepGOF-1 0.480, GRP 0.479, score test (asymptotic) 0.452, score test (bootstrap) 0.446, chi-square of map 0.400, Stukel 0.355, HL 0.229.

**Settings where one test leads by more than 2 paired SE** (`B_COUNTS.csv`), as DeepGOF-1 ahead / other test ahead / within 2 SE:

| compared with | count |
|---|---|
| score test (bootstrap) | 18 / 23 / 7 |
| chi-square of map | 31 / 10 / 7 |
| GRP | 16 / 17 / 15 |
| Stukel | 33 / 12 / 3 |
| HL | 43 / 0 / 5 |

**By shape**, 12 settings each, in the same order:

| shape | vs score test | vs chi-square of map |
|---|---|---|
| oscillation | 12 / 0 / 0 (mean +0.29) | 12 / 0 / 0 |
| bump | 6 / 0 / 6 (mean +0.08) | 9 / 0 / 3 |
| threshold | 0 / 11 / 1 (mean -0.09) | 9 / 0 / 3 |
| corner | 0 / 12 / 0 (mean -0.14) | 1 / 10 / 1 |

## Limiting power (`LIMIT_SCORE_TABLE2.csv`)
The ratio is the number of observations the score test needs for limiting power 0.5, divided by the number DeepGOF-1 needs, in the design of Table 2:

| departure | ratio |
|---|---|
| quadratic | 0.38 |
| interaction | 0.33 |
| cubic | 2.11 |
| local bump | 3.25 |
| threshold | 0.92 |

The same code reproduces the Stukel column of `local_power/FORMULA.json` to 6 decimals.

## Summary
1. Neither test dominates on non-polynomial departures.
2. The network wins on oscillating departures in every setting, and on the local bump it is ahead or tied.
3. The score test wins on thresholds and corner effects, and on quadratic and interaction terms (Tier 3).
4. The network beats the chi-square of its own map in 31 of 48 settings; all 10 settings where the chi-square wins are corner settings.
5. GRP matches DeepGOF-1 on mean power, leading in D2 and trailing in D10.
