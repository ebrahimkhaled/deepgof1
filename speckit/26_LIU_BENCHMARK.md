# 26 — GATE LIU: THE PUBLISHED-BENCHMARK REPLICATION (frozen 2026-08-22, before any data)

**Frozen.** Author-directed. The point is that our numbers land on a table a referee can compare to
a published one, in the journal we are targeting.

## 1. The benchmark being matched

Liu, Li, Chen, Härdle & Liang (2024), *Statistics and Computing* **34**:175,
"A comprehensive comparison of goodness-of-fit tests for logistic regression models."
Their DGPs, reproduced EXACTLY (their §4, C = departure parameter):

| setting | linear predictor | covariates |
|---|---|---|
| **1** | `-1.963 + 0.982*x1 + C*x1^2` | `x1 ~ U[-3,3]` |
| **2** | `-2 + 0.3*x1 + 0.3*x2 + 0.3*x3 + C*x1^2` | `x1 ~ U[-3,3]`, `x2 ~ N(0,1.5^2)`, `x3 ~ t4` |
| **3** | `-1.972 + 0.135*x1 + 1.117*x2 + C*x1*x2` | `x1 ~ U[-3,3]`, `x2 ~ Bin(.5)` |
| **4** | `0.3*x1 + 0.3*x2 + C*x1*x2` (NO intercept) | `x1, x2 ~ U[-3,3]` |

**C in {0, 0.25, 0.5, 0.75, 1.0}** (C = 0 is the null and gives each arm its size);
**n in {50, 100, 500}** — their sample sizes; 1000 replicates per cell; 2 corpora; B = 199; K = 6.
The fitted model always omits the C-term, as in their study.

**⭐ Setting 2 is the decisive cell.** Liu et al. report that all nine classical tests FAIL there —
"their power functions are not even increasing with C" — and only their projection test `proj`
gains power (>90% at C>0.9). The departure `C*x1^2` is nearly orthogonal to the fitted index
because x1 enters it with coefficient 0.3. **This is precisely the covariate-space corner the
residual map was built for (phase 11 / Figure 3), so it is where this project's central mechanism
claim is most directly testable.** It is also where a failure would be most informative.

## 2. Comparators — EXTERNAL, PUBLISHED tests only

`HL`, `HL-equalwidth`, `Pigeon-Heyse`, `Tsiatis`, `Xie`, `Stukel`, `Osius-Rojek`, `le-Cessie`
(+ their NA-robust Cauchy combination, `CCT-8`).

- **The author's own family is reported SEPARATELY and never as the comparator**: `EF`,
  `DEF.poly2/3/stukel` (EDGE), `A.poly`, `A.spline` (EDGE/EDGES, unpublished).
### AMENDMENT (2026-08-23): both missing comparators are ADDED

**Stute-Zhu — my error, corrected.** It was recorded as "absent, prohibitively slow". It is in the
package as `gof_stutezhu`, reachable via `run.all.gof(fit, tests="Stute-Zhu", include_slow=TRUE)`,
and costs **0.29 s per dataset**. The earlier judgement came from an `include_slow=TRUE` FULL-battery
call timing out, wrongly attributed to Stute-Zhu without testing it in isolation. **It is now a
first-class comparator.**

**`proj` (Escanciano RMEP) — implemented, not quoted.** Liu et al. report ~2.5 h per simulation
setting, which appeared to put it out of reach. It is not, because of one structural fact:

> **The A matrix depends only on X, never on y.** Compute it ONCE per dataset — O(n^3 p), vectorised
> over the third index — and every bootstrap replicate then costs only O(n^2) for `eps' A eps`.
> The 199 refits already exist for the residual map, so proj rides on the same loop.

Statistic: `T = eps-hat' A eps-hat / n^2` with
`A_ij = sum_r ||X_i - X_r|| ||X_j - X_r|| (pi - arccos(a_ijr))`,
`a_ijr = (X_i-X_r)'(X_j-X_r) / (||X_i-X_r|| ||X_j-X_r||)`, X standardised, intercept dropped.
Calibrated by the SAME model-based bootstrap as everything else, exact rank p-value at B = 199.

**VALIDATION REQUIRED BEFORE USE (frozen here):** our `proj` must reproduce Liu et al.'s reported
qualitative behaviour — high and MONOTONE power on Setting 2, where their nine classical tests
flatten. If it does not, the implementation is wrong and no proj number may be reported.
Liu's published curves are the external check; they are NEVER copied into our tables.

**Cost control:** if the n = 500 cells prove too slow to complete, proj is run at n = 50 and 100
only and its absence at n = 500 is stated per cell. Reps are never silently reduced.

## 3. ⚠️ The p = 1 map, declared before use

The residual map is a K x K grid over TWO covariate axes. **Setting 1 has one covariate**, so the
two-axis grid degenerates. **Declared construction for p = 1: a ONE-dimensional partition of
rank(x1) into K^2 = 36 quantile cells, reshaped row-major into 6 x 6** so the network input shape is
unchanged. This is a finer 1-D partition, not a 2-D map; it is reported as such and Setting 1's
result is never used to support a claim about covariate-SPACE detection.
For p >= 2 the standard rule applies (top two predictors by |beta_j| * sd(x_j), data-column order).
Setting 3's `x2` is binary, so its axis has two populated bands — expected, not a defect.

## 4. Criteria (frozen now)

| # | Criterion | Threshold |
|---|---|---|
| **LIU1** | deployed size at C = 0 | in [.035,.065] in every (setting, n, corpus) |
| **LIU2** | ⭐ Setting 2: power increases with C | rejection at C=1.0 minus at C=0 > **+.15**, mean over corpora — the cell where all nine classical tests fail |
| **LIU3** | mean matched-level gain over **CCT-8** across all four settings | > **+.03** |
| **LIU4** | gain over the best EXTERNAL test, per setting | reported for each; no threshold |
| **LIU6** | `proj` validation: Setting 2 rise, our implementation | > **+.60** and monotone in C, else the implementation is rejected and no proj number is reported |
| **LIU5** | across-corpus range | drawn/reported; no single-corpus number quoted |

## 5. Fallbacks

- **LIU2 passes** -> the paper has a genuine positive headline on a published benchmark: *the
  pretrained map detects the off-index departure that defeats the classical battery.* `proj` also
  detects it, so the honest claim is "a second, pretrained route to the same corner", never
  "the only test that works".
- **LIU2 fails** -> the central mechanism claim (Figure 3, phase 11) fails its most direct external
  test. That is reported prominently, and the covariate-space story is retired.
- **LIU3 fails but LIU2 passes** -> scope the claim to the off-index corner alone.

## 6. Forbidden

Copying `proj` numbers from Liu et al. into our tables — ours are run, or the cell is blank.
Reporting any proj number if LIU6 fails. Using our own EDGE family as the comparator. Changing C,
n, or the DGPs after seeing results. Single-corpus numbers. Silently reducing reps for cost.

## 6b. THE CLAIM THIS GATE IS BUILT TO SUPPORT (agreed with the author 2026-08-23)

> *The first pretrained goodness-of-fit test for logistic regression. It holds correct size on an
> external published benchmark; it detects the off-index departure that the partition family and
> Stukel miss; it is NOT the most powerful — le-Cessie and the author's own EDGE family exceed it.*

Every clause is a measurement, not a position: size from LIU1, the off-index detection from LIU2,
and the ranking from LIU4. **The concession is part of the claim, not an apology attached to it** —
it belongs in the abstract and the results table, never only in a limitations section.

## 7. Execution

`code/liu_grid.R` -> `code/liu_score.py` -> `code/figs/fig8_liu.py` + verdict appended here.
Deposit `results/liu/`.
