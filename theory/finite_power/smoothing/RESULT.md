# W4 — smoothing certificate: results (2026-09-24)

Theorem and proofs: `THEOREM_SMOOTH.md`. Code: `smooth_core.py`, `01_null_quantiles.py` (null draws and
critical-value bounds), `02_mc_net.py` (Monte Carlo at the net), `03_tail.py` (CROWN interval chain to
$t_g$, then Lemma 2), `04_truth.py` (the truth by bisection, with SEs), `05_certify.py` (assembly →
`SMOOTH_CERT.json`), `06_fill_inserts.py` (fills the TeX inserts and the tables below). Driver
`run_all.py`; everything resumes from the jsonl files.

Revised 2026-09-24 after `ADVERSARIAL_REVIEW.md` (M1–M5, m1–m9): certified $n$ now rounded up, per-rule
vs joint probability stated, global constants stated in the ESM insert, fixed tags S.12a/S.12b, matched
truth everywhere, three slips corrected. Manuscript edits for the merge: `MERGE_NOTES.md`.

## What is certified (and what is not)

Every number below is a statement about the **Gaussian model of the map** used by the paper's
Proposition S.12: $m=tv+e$, $e\sim N(0,I_{36})$, $S$ the frozen network, $n=4t^{2}$. It is **not** a
guarantee for the finite-$n$ Bernoulli law of the deployed test (referee point B2): the deployed map is
Gaussian only in the limit, its null covariance is $\Sigma$ rather than $I$, and its critical value is a
bootstrap rank. The theorem is proved for a general singular $\Sigma$ along rays in $\mathrm{range}(\Sigma)$
(the local-model case of Theorem S.9), but the computation is for $\Sigma=I$ so that it compares cell by
cell with the paper's Table of Proposition S.12.

Two rejection rules are certified in the model: the **fixed threshold** $S>\bar t=2.4865$ (the paper's
table), and the **bootstrap rule** ($B=199$, $\alpha=.05$) through Dufour's binomial form of its power and
a layer cake over 14 critical-value levels (Corollary in `THEOREM_SMOOTH.md`).

The certificate is a Monte Carlo certificate: it holds with probability $\ge1-\eta$ over its own draws,
$\eta=.001$ for the net (Bonferroni over the net points; over net points × levels for the bootstrap rule),
plus $.001$ for the critical values ($\bar t$ carries the paper's own $.001$; the 14 levels of the
bootstrap rule share $.001$ from 2,000,000 fresh null draws). **The level $.998$ is per direction and per
rule.** Jointly (union bound): all five fixed-rule rows $\ge1-5(.001)-.001=.994$; all twenty cells (two
rules × two targets × five directions; the two targets of one rule and direction share one CP event)
$\ge1-10(.001)-.001-.001=.988$. The tail beyond $t=250$ is deterministic
(CROWN over balls of radius $8.74$, $P\{\|e\|>8.74\}=10^{-4}$, chained to $t_g$, then the certified
constants $c\le100.70$ (exact 100.6965; offset range $[-100.70, 82.83]$) and $L\le1735.64$ (exact
1735.636) with Lemma 2). Both constants come from interval propagation in `../03_certify.py` (offsets
through every layer for $c$; the gradient through every activation pattern for $L$, so the same $L$ holds
for $S$ and $S_\infty$) and were checked on 240,000 sampled maps at six scales plus a gradient attack
($c_{emp}=7.51$, $L_{emp}=3.82$). The ESM did not state them before; `INSERT_ESM.tex` now does. The
chains and $t_g$ use $c^\circ=\max(\bar t,\hat c_1,\dots,\hat c_{14})=5.8321$ in place of $\bar t$ (one
chain serves both rules). Nothing in the certified inequality is a point estimate.

**Rounding rule.** Every certified number is rounded in the safe direction: $n_{cert}=\lceil4t_{cert}^2\rceil$,
$t_{cert}$ rounded up to the printed digit, $c$, $L$, $t_g$ rounded up (`06_fill_inserts.py` checks all 40
printed certified numbers). Truth estimates and ratios use ordinary rounding; ratios are from unrounded values.

## Main table: fixed threshold, power .95 (the paper's setting)

`n_CROWN` is the paper's usage-aware relaxation certificate (Online Resource 1, Table
"Certified against true sample sizes"); the CROWN ratio is recomputed against the truth below.

| direction | true t* (SE) | true n* | n* CROWN (paper) | ratio CROWN | t_cert smoothing | n* smoothing | ratio smoothing [2 SE] |
|---|---|---|---|---|---|---|---|
| step | 4.274 (0.006) | 73 | 3,789 | 51.9 | 4.327 | 75 | **1.02** [1.02, 1.03] |
| diffuse | 5.105 (0.006) | 104 | 4,311 | 41.4 | 5.154 | 107 | **1.02** [1.01, 1.02] |
| stripe | 5.547 (0.006) | 123 | 5,881 | 47.8 | 5.606 | 126 | **1.02** [1.02, 1.03] |
| saddle | 7.156 (0.007) | 205 | 8,044 | 39.3 | 7.213 | 209 | **1.02** [1.01, 1.02] |
| random | 13.529 (0.013) | 732 | 13,918 | 19.0 | 13.701 | 751 | **1.03** [1.02, 1.03] |

($t_{cert}$ rounded up, $n^*$ smoothing $=\lceil4t_{cert}^2\rceil$; ratios from unrounded values.)

**The smoothing certificate is within 1.6–2.6 % of the truth in $n$ (0.8–1.3 % in $t$), against 19–52× for
the relaxation certificate.** In every cell the certified $t$ exceeds the true $t$ (as it must with probability
$\ge.998$), by 8–13 SE of the truth (9.5, 8.3, 9.7, 8.8, 13.1), so the excess is real slack (net spacing +
Clopper–Pearson), not noise.

The truth here is the Monte Carlo $t$ at which power reaches .95 **with the same threshold $\bar t$ the
certificate uses**, from 200,000 fresh draws per bisection step, SE by the delta method. The paper's
"true $n^*$" (72, 103, 122, 204, 729) used the Monte Carlo point estimate $2.4616$ of the null quantile
rather than $\bar t$, and 20,000 draws; with that threshold the fresh estimate is $t=4.257, 5.078, 5.531,
7.154, 13.487$ (SE $.005$–$.013$), consistent with the paper's $4.241, 5.066, 5.516, 7.134, 13.498$ at its
SE (about $.016$). The paper's CROWN ratios (ESM tab:certS prints 52.7, 42.0, 48.3, 39.5, 19.1 from
unrounded $t$; 52.6, 41.9, 48.2, 39.4, 19.1 from the rounded $n$) become 51.9, 41.4, 47.8, 39.3, 19.0
against the $\bar t$-truth. Note that the 2,000,000-draw point estimate of the
$.95$ null quantile is $2.4753$, above the paper's $2.4616$ (400,000 draws); the paper's $\bar t=2.4865$
remains a valid $99.9\%$ upper bound (my own from the fresh draws is $2.4855$).

**Consistency with the ESM (review M3).** Every certificate, Proposition S.12's included, is for the
threshold $\bar t$, so the matched truth is the $\bar t$-truth ($n=73, 104, 123, 205, 732$). ESM Table
tab:certS and the "19–53" / "21–57" / "47.4" / "5.52" statements must move to it at merge time; the exact
edits are in `MERGE_NOTES.md`.

## Secondary tables

**Fixed threshold, power 1/2**

| direction | true t* (SE) | true n* | t_cert | n* smoothing | ratio [2 SE] |
|---|---|---|---|---|---|
| step | 2.312 (0.003) | 21.4 | 2.342 | 22 | **1.03** [1.02, 1.03] |
| diffuse | 3.048 (0.004) | 37.2 | 3.088 | 39 | **1.03** [1.02, 1.03] |
| stripe | 3.341 (0.004) | 44.6 | 3.389 | 46 | **1.03** [1.02, 1.03] |
| saddle | 4.689 (0.005) | 88.0 | 4.774 | 92 | **1.04** [1.03, 1.04] |
| random | 8.839 (0.008) | 312.5 | 9.001 | 325 | **1.04** [1.03, 1.04] |

**Bootstrap rule B = 199, power .95** (truth: Dufour form with $W$ against the 2,000,000 null draws)

| direction | true t* (SE) | true n* | t_cert | n* smoothing | ratio [2 SE] |
|---|---|---|---|---|---|
| step | 4.302 (0.005) | 74.0 | 4.429 | 79 | **1.06** [1.05, 1.06] |
| diffuse | 5.145 (0.006) | 105.9 | 5.268 | 112 | **1.05** [1.04, 1.05] |
| stripe | 5.594 (0.006) | 125.1 | 5.730 | 132 | **1.05** [1.04, 1.05] |
| saddle | 7.232 (0.007) | 209.2 | 7.432 | 221 | **1.06** [1.05, 1.06] |
| random | 13.657 (0.013) | 746.0 | 14.083 | 794 | **1.06** [1.06, 1.07] |

**Bootstrap rule B = 199, power 1/2**

| direction | true t* (SE) | true n* | t_cert | n* smoothing | ratio [2 SE] |
|---|---|---|---|---|---|
| step | 2.326 (0.004) | 21.6 | 2.438 | 24 | **1.10** [1.09, 1.11] |
| diffuse | 3.066 (0.004) | 37.6 | 3.206 | 42 | **1.09** [1.09, 1.10] |
| stripe | 3.356 (0.004) | 45.1 | 3.507 | 50 | **1.09** [1.09, 1.10] |
| saddle | 4.719 (0.005) | 89.1 | 4.914 | 97 | **1.08** [1.08, 1.09] |
| random | 8.886 (0.009) | 315.9 | 9.305 | 347 | **1.10** [1.09, 1.10] |

(In all tables: $t_{cert}$ rounded up, $n^*$ smoothing $=\lceil4t_{cert}^2\rceil$, true $n^*$ an estimate.
The B = 199 $t_{cert}$ are grid points of spacing $\Delta t=.001$, at which the lower bound
$\underline\beta\ge\beta+0.2\Delta t$; between grid points it is $0.4$-Lipschitz.)

The bootstrap rule costs more (5–10 %) because the layer cake discards the mass of $W$ strictly between
consecutive levels and because 14 levels share the Bonferroni budget; at power 1/2 the layer $w\in(.045,.05]$
carries weight .46 and is where the loss sits. Finer levels would shrink it. (Efficiency note, review m7:
the levels $w\le.007$ carry total weight about $10^{-3}$ but take 6/14 of the Bonferroni budget; fewer
of them, or more levels in $(.03,.05]$, would tighten the 1.08–1.10 at power 1/2. Not re-run: the
certificate is valid as it stands, and changing the levels now would be a post-hoc design change.)

Sensitivity to the Monte Carlo budget (fixed threshold, .95): certified $t$ (rounded up) at $\eta=.001$ vs
$.01$ is 4.327/4.323, 5.154/5.150, 5.606/5.602, 7.213/7.210, 13.701/13.697 — under 0.1 %. The certificate is limited by
the net spacing ($0.1$, i.e. $0.05$ on the $\Phi^{-1}$ scale at the worst point) and the draw count
($M=10^{5}$), not by the confidence level.

## Monte Carlo budget and runtime

| direction | net points J | draws (M) | net range | tail intervals (CROWN calls) | t_global (rounded up) | min CROWN lb in chain | max Λ |
|---|---|---|---|---|---|---|---|
| step | 260 | 6.3 | [1.27, 250] | 19 (57) | 6,944 | 247.9 | 3.28 |
| diffuse | 271 | 7.4 | [1.52, 249] | 21 (63) | 10,430 | 97.4 | 3.44 |
| stripe | 278 | 8.0 | [1.66, 250] | 22 (62) | 8,835 | 39.5 | 3.42 |
| saddle | 300 | 10.0 | [2.14, 250] | 31 (75) | 12,862 | 5.95 | 3.66 |
| random | 387 | 18.2 | [4.05, 249] | 50 (100) | 20,908 | 10.7 | 3.65 |

$t_g=(c^\circ+c+LR)/\kappa$ with $c^\circ=5.8321$ (the largest critical value in use, $\ge\bar t$), not $\bar t$.

Net design (fixed before any draw): spacing 0.1 with $M=10^{5}$ on $[0.3,1.5]\times$ the paper's $t^*$;
spacing 0.5 with $M=2\times10^{4}$ to $4\times t^*$; spacing 1.25 with $M=4{,}000$ to $t=250$. Null draws:
2,000,000 (131 s). Truth: about 4 million draws per direction (5 targets × ~17 bisection steps × 200,000).
Wall time on a shared box (12–14 R processes running throughout, 4–8 Python workers with 2 threads
each): about 4 min for the null draws, 15–40 min per net, 5–9 min per tail, 45–75 min per truth run;
the whole pipeline about 3 h wall, of which the truth runs (not part of the certificate) were half.
Forward-pass throughput was 3,000–9,000 draws/s per worker under contention. The pipeline was killed
twice (a session exit and a power-off) and resumed from the jsonl files without loss; no line was
truncated (`clean_jsonl.py`).

## Caveats (honest list)

1. **Model, not test.** Gaussian map, identity covariance, fixed direction, $n=4t^{2}$, and a fixed or
   model-bootstrap critical value. Say so wherever the numbers appear, as the paper does for S.12.
2. **Monte Carlo certificate.** Holds with probability $\ge.998$ over the draws **per direction and per
   rule** ($\ge.994$ for the five fixed-rule rows jointly, $\ge.988$ for all twenty cells), treating the
   pseudo-random generator as ideal; the draws are on disk (`net_*.jsonl`, seeds recorded) and the
   certificate is recomputable from the counts alone.
3. **The truth is itself a Monte Carlo estimate** (SE $.003$–$.013$ in $t$); the ratio intervals in the tables
   move the truth by 2 SE. The certified $t$ is a bound, not an estimate, and is always above the truth.
4. **The tail uses CROWN soundness and the certified global constants**, both inherited from the paper's
   `03_certify.py`; the CROWN chain covers $[250,t_g]$ with every lower bound above $c^\circ=5.832$. The
   smallest lower bounds per direction are 247.9, 97.4, 39.5, **5.95** (saddle: margin 0.12 over $c^\circ$)
   and 10.7 (random). CROWN runs in double precision without outward rounding; the smallest margin,
   0.12, is many orders of magnitude above rounding error, so this does not matter, but the saddle
   margin is thin in absolute terms.
5. **The bootstrap-rule truth** uses the empirical null cdf from 2,000,000 draws inside the Dufour formula;
   its SE ignores the null-cdf error, which is of order $10^{-3}$ in $W$ and negligible at these $t$.
6. The relaxation certificate (S.12) remains in the paper as a diagnostic of verification tools; the
   smoothing certificate is the usable one. Neither replaces the local-power theorem or simulation.
