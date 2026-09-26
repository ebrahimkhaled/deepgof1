# Level-rate study — result (2026-09-24)

**One line.** At every $n\ge100$ the deployed test's size, pooled over the three null designs, is within
two Monte Carlo standard errors of the nominal level ($|{\rm size}-\alpha|\le.0025$ at $\alpha=.05$, SE
$.0014$); single designs at $n\ge100$ exceed $\alpha$ by at most $.0038$, and the three single-design cells
more than two SE from $\alpha$ (D4 $n=100$, $200$; D3 $n=400$; at $\alpha=.05$, largest $|{\rm err}|=.0067$)
all lie *below* it. So the *rate* of the level error cannot be
estimated from 7,600 replicates per cell: the fitted log–log slope is $-0.34$ with SE $0.40$ (pooled,
$\alpha=.05$) and its replicate bootstrap interval $[-1.34,+0.54]$ contains both $-1/2$ and $-1$. The one
cell with a detectable error is the rare-event design at $n=50$ (15% events, about 7–8 events per
dataset): size $.058$ (SE $.003$) at $\alpha=.05$ and $.116$ (SE $.004$) at $\alpha=.10$.

## Design
* Released `ebrahim.gof::deepgof1(fit, B = 199L, K = 6L)` (installed 2.8.0; its `deepgof1()` selects the
  two strongest covariates, the same rule as CRAN 2.7.0 and the paper; verified from the installed
  function body). The analyst's own plain `glm` fit is passed; no separation policy (as deployed).
* Null designs of the pre-declared grid (Simulation I, generators copied verbatim from
  `code/universal_design.R`): D1 ($p=2$, 50% events), D4 ($p=2$, 15% events), D3 ($p=3$, AR(1) $\rho=.7$).
  Fresh covariates and a fresh random $\beta$ per replicate.
* $n\in\{50,100,200,400,800\}$, 7,600 null replicates per cell (SE $\le.0025$ at $\alpha=.05$), 15 cells,
  114,000 replicates; per-replicate p-values in `RATE_RAW/` (1,140 chunk files of 100, all verified
  complete, no NA, no duplicate seeds). Seeds deterministic per chunk (`run_rate.R`).
* $\alpha=.05$ and $.10$ are both on the p-value grid $k/200$ ($k=10,20$), so the comparison is with
  the attainable level, not a floor of it.
* The 20,000-replicate version was not run: at $\approx3.5$ s per replicate it needed $\approx26$ h of wall
  time against the 24 h budget; 7,600 was the pre-set SE target.

## Size, |size − α| and SE (7,600 replicates per cell; pooled rows 22,800)

| design | $n$ | size at .05 | SE | \|err\| | $z$ | size at .10 | SE | \|err\| | $z$ | sep. flag |
|---|---|---|---|---|---|---|---|---|---|---|
| D1 | 50 | .0514 | .0025 | .0014 | 0.6 | .1018 | .0035 | .0018 | 0.5 | .0000 |
| D1 | 100 | .0514 | .0025 | .0014 | 0.6 | .1017 | .0035 | .0017 | 0.5 | .0000 |
| D1 | 200 | .0504 | .0025 | .0004 | 0.2 | .1037 | .0035 | .0037 | 1.1 | .0000 |
| D1 | 400 | .0522 | .0026 | .0022 | 0.9 | .1024 | .0035 | .0024 | 0.7 | .0000 |
| D1 | 800 | .0501 | .0025 | .0001 | 0.1 | .0946 | .0034 | .0054 | −1.6 | .0000 |
| D4 | 50 | **.0578** | .0027 | .0078 | **2.9** | **.1162** | .0037 | .0162 | **4.4** | .0100 |
| D4 | 100 | .0433 | .0023 | .0067 | −2.9 | .0984 | .0034 | .0016 | −0.5 | .0001 |
| D4 | 200 | .0449 | .0024 | .0051 | −2.2 | .0972 | .0034 | .0028 | −0.8 | .0000 |
| D4 | 400 | .0488 | .0025 | .0012 | −0.5 | .0979 | .0034 | .0021 | −0.6 | .0000 |
| D4 | 800 | .0484 | .0025 | .0016 | −0.6 | .0986 | .0034 | .0014 | −0.4 | .0000 |
| D3 | 50 | .0491 | .0025 | .0009 | −0.4 | .1014 | .0035 | .0014 | 0.4 | .0000 |
| D3 | 100 | .0476 | .0024 | .0024 | −1.0 | .1008 | .0035 | .0008 | 0.2 | .0000 |
| D3 | 200 | .0538 | .0026 | .0038 | 1.5 | .1020 | .0035 | .0020 | 0.6 | .0000 |
| D3 | 400 | .0449 | .0024 | .0051 | −2.2 | .0970 | .0034 | .0030 | −0.9 | .0000 |
| D3 | 800 | .0480 | .0025 | .0020 | −0.8 | .0976 | .0034 | .0024 | −0.7 | .0000 |
| pooled | 50 | .0528 | .0015 | .0028 | 1.9 | .1065 | .0020 | .0065 | 3.2 | |
| pooled | 100 | .0475 | .0014 | .0025 | −1.8 | .1003 | .0020 | .0003 | 0.2 | |
| pooled | 200 | .0497 | .0014 | .0003 | −0.2 | .1010 | .0020 | .0010 | 0.5 | |
| pooled | 400 | .0486 | .0014 | .0014 | −1.0 | .0991 | .0020 | .0009 | −0.5 | |
| pooled | 800 | .0489 | .0014 | .0011 | −0.8 | .0969 | .0020 | .0031 | −1.6 | |

"sep. flag" = fraction of observed fits with a non-finite or $|{\rm coef}|>8$ coefficient. Bootstrap
refit failures (scored $+\infty$ by the released code): zero in all 114,000 replicates.

Joint test that the pooled size equals $\alpha$ at all five $n$ ($\sum z^2$ over five independent
cells, **5** degrees of freedom): $8.3$ at $\alpha=.05$ ($p\approx.14$) and $13.0$ at $\alpha=.10$
($p\approx.02$), the latter carried by $n=50$. (Corrected 2026-09-24: the first version called this a
heterogeneity test on 4 df, $p\approx.08$ and $.01$; the statistic is a sum of five squared $z$'s about
$\alpha$, not about their mean. `analyze_rate.py` still labels the field `chi2_4df`; the value is right,
the label is not.) The size is averaged over random covariates and coefficients, so it is the
unconditional size, not the conditional size of Theorem R.

## The fitted rate (log|error| on log n)

| series | $\alpha$ | WLS slope (SE) [delta-method weights] | OLS slope (SE) | replicate-bootstrap 95% of OLS slope | cells with \|err\| < 2 SE |
|---|---|---|---|---|---|
| pooled | .05 | −0.34 (0.40) | −0.35 (0.42) | [−1.34, +0.54] | 5 of 5 |
| pooled | .10 | −0.32 (0.25) | −0.06 (0.62) | [−1.16, +0.41] | 4 of 5 |
| D1 | .05 | +0.22 (0.93) | −0.63 (0.50) | [−3.70, +2.87] | 5 of 5 |
| D1 | .10 | +0.39 (0.54) | +0.36 (0.15) | [−0.76, +1.32] | 5 of 5 |
| D4 | .05 | −0.42 (0.34) | −0.71 (0.20) | [−3.45, −0.03] | 2 of 5 |
| D4 | .10 | −1.08 (0.49) | −0.66 (0.36) | [−1.61, +0.09] | 4 of 5 |
| D3 | .05 | +0.28 (0.58) | +0.33 (0.29) | [−0.95, +3.24] | 4 of 5 |
| D3 | .10 | +0.25 (0.92) | +0.34 (0.19) | [−1.07, +3.06] | 5 of 5 |

**Verdict on the rate: not identifiable here, and I do not conclude one.** The slope is a regression on
$\log|{\rm err}|$ where, in 34 of the 40 cell–level pairs (15 design cells and 5 pooled rows, each at two
levels; the pooled rows reuse the design replicates, so the 40 are not independent; 25 of 30 for the
design cells alone), $|{\rm err}|$ is inside two standard errors of zero; its value is then mostly the log of Monte Carlo noise, which is why the OLS standard errors (which
assume the points are informative) are smaller than the honest replicate-bootstrap spread, and why the
D1 and D3 slopes are *positive*. The pooled bootstrap interval at $\alpha=.05$, $[-1.34,+0.54]$,
contains $-1/2$ (Beran's first-order rate), $-1$ (second-order), and $0$.

**What would be needed.** Anchoring both candidate rates at the pooled $n=50$ error ($.0028$), the
$n^{-1/2}$ and $n^{-1}$ curves differ by $.0006$, $.0007$, $.0006$ and $.0005$ at $n=100$, $200$, $400$
and $800$. Separating them at two standard errors needs about $530{,}000$, $390{,}000$, $460{,}000$ and
$700{,}000$ null replicates per $n$, i.e. $17$–$31$ times this study's $22{,}800$ pooled replicates per $n$
(about $2.1$ million replicates in all, roughly $2{,}000$ CPU-hours at $3.5$ s per replicate; corrected
2026-09-24 from "50–90 times, 5,000 CPU-hours", which compared with a single design's 7,600).
And that presumes the sign of the error is stable in $n$, which the table does not support (D4 changes
sign between $n=50$ and $n=100$). A rate study of the level error of this test is therefore not a
feasible route to a second-order claim; the theory (THEOREM_RATE.md) gives $n^{-1/2}\log n$ and the
data give an *upper bound on the error*, not a rate.

## What the data do establish
1. **An upper bound on the level error from $n=100$.** Pooled over designs, every $|{\rm size}-\alpha|$ at
   $n\ge100$ is $\le.0025$ at $\alpha=.05$ (95% upper limit $\approx.005$) and $\le.0031$ at $\alpha=.10$
   (upper limit $\approx.007$). The largest single-cell excess at $n\ge100$ is $+.0038$ (D3, $n=200$, $z=1.5$).
2. **The one real excursion is small-$n$, rare-event.** D4 at $n=50$: $+.008$ at $\alpha=.05$ ($z=2.9$) and
   $+.016$ at $\alpha=.10$ ($z=4.4$). It is *not* the observed fits flagged for separation (the 76
   flagged fits, any coefficient non-finite or $|{\rm coef}|>8$, reject at $1.3\%$; the 7,524 others at
   $.058$ / $.117$, recomputed from `RATE_RAW/`). Two limits of that statement: the flag misses
   quasi-separation with coefficients below 8, and separation inside the bootstrap replicates was not
   recorded. With about 7.5 events in 50 records the cells of the $6\times6$ map hold about $0.2$ events
   each, where one expects the Gaussian surrogate that Theorem R and Proposition R.6 rest on to be far
   from the lattice truth; this is an interpretation, not something the study tests. The excess is
   confined to $n=50$: at $n=100$ (about 15 events) the same design is *conservative* by $.007$
   ($z=-2.9$), still conservative at $n=200$ ($z=-2.2$ at $\alpha=.05$), and within noise from $n=400$.
3. **The level is on the attainable level, not below it.** With $\alpha$ on the grid $k/200$ the
   asymptotic size is exactly $\alpha$ (Proposition S.3(i)); the pooled means at $n\ge100$, $.0475$–$.0497$
   at $\alpha=.05$, are consistent with that. Pooling all $n\ge100$ (91,200 replicates) gives $.0487$
   (SE $.0007$, $z=-1.9$): at most a slight conservatism, carried by the D4 and D3 cells above.

## Consistency with the paper's earlier numbers
* Sixty-cell grid (Section 5 of the main text; $n\le200$, Firth in the harness): mean deployed size at
  $\alpha=.05$ in band $[.035,.065]$ in 58 of 60 cells. The present cells at $n\le200$ ($.043$–$.058$) fall
  in the same band, without Firth.
* Theorem G coverage study (300 datasets per $n$): measured sizes $.030/.067/.063$ at $n=50/100/200$, SE
  $.013$ — those were within noise of $.05$ and are now superseded by SE $.0025$.
* Proposition S.7's measured tail gap ($q_{90}$ $.030$ at $n=50$, $.007$ at $n=200$) bounds the size
  *inflation* by $\varepsilon+\eta$; the sizes here are far inside that bound, as the bound's own text says.

## Runtime
About $100$ CPU-hours: 390 chunks on 4 cores (until the shared simulation finished), then 147 and 603
chunks on 12 cores (last stint 3.4 h). Logs: `run_rate_4cores.log`, `run_rate_12cores.log`,
`run_rate_12cores_b.log`. The study is **complete**; `python analyze_rate.py` regenerates `RATE.json`
and the tables above from `RATE_RAW/`.

## Files
`rate_common.R` (designs, one replicate), `run_rate.R` (resumable driver), `time_one.R`,
`analyze_rate.py`, `RATE.json`, `RATE_RAW/` (114,000 p-values), `THEOREM_RATE.md`, `INSERT_MAIN.tex`,
`INSERT_ESM.tex`.
