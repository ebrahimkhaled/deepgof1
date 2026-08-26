# 28 — THEORY PACKAGE: validity + conditional consistency with a computable certificate

**Status: PROOFS WRITTEN 2026-08-25 — `theory/theory_supplement.tex`, 8 pp, xelatex 0 errors.**
Numbering there: Lemma S.1 = Lemma 0 · Lemma S.2 = Lemma 1 · Proposition S.3 = asymptotic
validity (NEW — formalizes the Beran connection: our p-value IS the first-stage prepivoted
root computed by Monte Carlo, so ERP is O(n^-1/2) generally and O(n^-1) under asymptotic
pivotality) · Theorem S.4 = Theorem 2 · Corollaries S.5/S.6 = Corollaries 3/4.

⚠️ **REFINEMENT found while proving (supersedes Cor. 4 below):** “undetectable by ANY
map-based test” was an overclaim. The provable statement (Cor. S.6): with mu = 0 the map is
uniformly TIGHT, so no map-based test can be CONSISTENT — but nontrivial asymptotic power
below 1 can persist through second-order (covariance) differences, and the limit power has a
computable expression. Use the S.6 wording in the manuscript, never the old one.

**Original pre-declaration (engine verified before proof prose):** Companion scripts: `code/theory_recession_check.R`, `code/theory_sphere_check.R`,
`code/theory_blindcone_probe.R`. Deposit: `results/theory_blind_direction.csv`.

## 0. What this buys, in one paragraph

BAGofT (JASA 2023) has four theorems; Kim et al. (Stat&Comput 2025) has zero and says so; Simić
(2021) says "we do not have theoretical proofs" outright. This package gives DeepGOF-1 (i) an
honest exact-validity statement with the precise condition under which it is exact, (ii) a
consistency theorem whose key condition is **computable from the frozen weights in one forward
pass**, and (iii) a formal description of the blind cone — which turns out to be sub-resolution
oscillation, the same boundary every binned instrument has. No claim anywhere exceeds what is
proved or computed.

## 1. Setup and notation

Data $(x_i, y_i)$ iid, $y\mid x \sim \mathrm{Bern}(\pi^*(x))$. Null model: logistic
$\pi_\beta(x)$. QMLE $\hat\beta$. Map axes = the two covariates maximizing
$|\hat\beta_j|\,\hat\sigma_j$; cells $C_{k\ell}$ = the $K\times K$ grid over empirical ranks
($K=6$); cell statistic $m_{k\ell} = \sum_{C_{k\ell}}(y_i-\hat\pi_i)\big/\sqrt{\sum_{C_{k\ell}}
\hat\pi_i(1-\hat\pi_i)}$. Shipped scaler $(a, s)$; score $S(z)$, $z=(m-a)/s$, with $S$ the frozen
CNN (affine layers, ReLU, max-pool, global max/mean). Bootstrap: $B$ refits from $\hat\pi$;
$p = (1+\#\{S_b \ge S_{obs}\})/(B+1)$, $B = 199$.

## 2. Assumptions

- **(A1)** iid sampling; $0<\pi^*<1$ a.s.; the selected covariates have continuous marginals and
  every population quantile cell has positive mass.
- **(A2)** White (1982) QMLE regularity: unique KL projection $\beta^*$,
  $\sqrt n(\hat\beta-\beta^*)=O_p(1)$.
- **(A3)** Axis stability: the population criterion $|\beta^*_j|\sigma_j$ has a unique top-two
  set.
- **(A4)** Map-visible misfit: $\mu \ne 0$, where
  $\mu_{k\ell} = \mathbb E[(\pi^*-\pi_{\beta^*})\,\mathbf 1_{C_{k\ell}}]\big/
  \sqrt{\mathbb E[\pi_{\beta^*}(1-\pi_{\beta^*})\,\mathbf 1_{C_{k\ell}}]}$
  (cells = population quantile rectangles of the selected pair).
- **(A5)** Certificate: $S_\infty(\mu/s) > 0$ (Lemma 0 below; positively homogeneous, so only the
  direction matters).

## 3. The results

**Lemma 0 (recession function, closed form).** For $S$ built from affine maps, ReLU, max-pool and
global max/mean, $S_\infty(u) := \lim_{t\to\infty} S(z_0+tu)/t$ exists for all $u, z_0$, does not
depend on $z_0$, is continuous and positively homogeneous, and **equals the same network with all
biases set to zero**, evaluated at $u$. Moreover $\sup_z |S(z)-S_\infty(z)| \le C < \infty$ (bias
perturbations propagate through Lipschitz layers).
*Proof:* structural induction over layers ($\mathrm{relu}(at+b)/t \to \mathrm{relu}(a)$, max of
finitely many affine functions likewise). **Numerically verified: the closed form matches the
measured large-$t$ slope to 2.2e-16 on five directions.**

**Lemma 1 (pivotality ⇒ exactness).** If the law of the score under $P_\theta$ (generate → refit
→ map → score) is the same for every $\theta$ in the null family, then under $H_0$ the vector
$(S_{obs}, S_1,\dots,S_B)$ is iid, and $P(p \le k/(B+1)) = k/(B+1)$ exactly, ties handled
conservatively by the $\ge$ in the count (Besag–Clifford 1989).
*Why it matters:* this is the precise statement of what "deployment-matched training" aims at.
The bootstrap scores are independent of the data whenever the score law is $\theta$-free, so the
rank argument is exact — **the phase-24 AUC diagnostic (bootstrap-null vs true-null, AUC .501) is
the empirical check of this lemma's hypothesis**, rehabilitated with its correct meaning. Away
from exact pivotality the test retains the classical parametric-bootstrap asymptotic validity
(Beran 1988); the observed deployed sizes (58/60 cells in gate S3, 23/24 external, DEEP .0506 on
GATE LIU) are the finite-sample evidence.

**Theorem 2 (conditional consistency).** Under (A1)–(A5), for any fixed $\alpha \ge 1/(B+1)$,
$P(p \le \alpha) \to 1$ as $n \to \infty$.
*Proof sketch:* (i) axis selection is eventually the population pair (A2, A3, LLN). (ii)
empirical rank cells converge to quantile rectangles; counts $= n\,p_{k\ell}(1+o_p(1))$. (iii)
decompose $\sum_{C}(y-\hat\pi) = \sum_C (y-\pi^*) + \sum_C(\pi^*-\pi_{\beta^*}) +
\sum_C(\pi_{\beta^*}-\hat\pi)$: CLT term $O_p(\sqrt n)$, drift term
$n\,\mathbb E[\delta \mathbf 1_C](1+o_p(1))$, estimation term $O_p(\sqrt n)$ by (A2) and
smoothness; hence $m/\sqrt n \to_p \mu$. (iv) $S(z_n) = S_\infty(z_n) + O(1)$ (Lemma 0 bound)
$= \sqrt n\, S_\infty(z_n/\sqrt n) + O(1)$, and $z_n/\sqrt n \to \mu/s$ with $S_\infty$
continuous, so $S_{obs} = \sqrt n\,(S_\infty(\mu/s) + o_p(1)) \to +\infty$. (v) under the
bootstrap law the model is correctly specified with parameter $\hat\beta$, so each
$S_b = O_p(1)$; $B$ fixed ⇒ $\max_b S_b = O_p(1)$. (vi) hence $p = 1/(B+1) \le \alpha$
with probability → 1. ∎

**Corollary 3 (the blind cone is computable — and small).** $\mathcal B = \{u : S_\infty(u)\le
0\}$ is a closed cone computable from the frozen weights. For DeepGOF-1's shipped weights:

- 200,000 + 50,000 random unit directions: **zero hits** (spherical measure < 7.4e-5, 95%
  Clopper–Pearson);
- multi-start descent finds min $\approx -0.281$: the cone is nonempty but thin;
- **the worst direction is sub-resolution oscillation** — lag-1 spatial autocorrelation $-0.40$,
  no row/column structure — a checkerboard at exactly the cell scale. This is the resolution
  limit EVERY binned statistic shares (cf. EDGE's osc4/sawtooth concession), now formalized;
- natural departure directions carry margins an order of magnitude above zero:

| direction | $S_\infty$ | mirror |
|---|---|---|
| step / threshold | **2.200** | 1.075 |
| curvature stripe (omitted $x^2$) | **1.729** | 1.644 |
| local bump | 1.506 | 1.493 |
| diffuse drift | 1.464 | 1.463 |
| interaction saddle | 1.188 | 1.141 |

**Corollary 4 (representation limit).** If $\mu = 0$ the drift vanishes for every statistic
computed through this map: cell-balanced misfit is undetectable by ANY map-based test. This is
the formal statement of the p=3 projection cost (a discarded covariate makes $\mu$ small or
zero) and of the probit finding ($\mu$ tiny ⇒ consistency holds but the constant is poor — no
contradiction with the observed small-$n$ blindness).

## 4. Practitioner reading of the certificate

Given ANY hypothesized departure $\delta$: simulate once at large $n$ to estimate $\mu$, take one
zero-bias forward pass, read the sign. Positive ⇒ the test will catch it eventually; the
magnitude is the drift constant. This is a checkable, per-alternative consistency certificate —
something neither BAGofT's rate conditions nor Kim et al.'s cut-off construction provides.

## 5. Honesty ledger

- Consistency is **conditional** on (A5); the omnibus version is FALSE for the shipped weights
  (min $S_\infty = -0.281$). Stated, with the cone measured and characterized.
- Exactness is **conditional** on pivotality (Lemma 1); away from it we claim only Beran-type
  asymptotic validity plus the measured sizes. No "exact by construction" without the condition.
- No rates are claimed for Theorem 2.
- $\mu = 0$ alternatives are conceded as undetectable at the representation level (Cor. 4).
- The certificate is exact for the *limit*; finite-$n$ power is what GATE LIU / S3 measured —
  the theory does not replace those tables, it explains them.

## 6. Remaining work to publication grade

1. Write the supplement proofs in full (Lemma 0 induction; Theorem 2 steps (i)–(vi) with the
   rank-cell empirical-process argument spelled out; Lemma 1 exchangeability).
2. Cite-verify against White (1982), Beran (1988), Besag–Clifford (1989) — PDFs requested from
   the author 2026-08-25 for the MD library.
3. Adversarial verification pass (multi-agent) on the proofs once drafted.
4. Fold into the manuscript: Lemma 1 + Theorem 2 in the main text (~1.5 pp), proofs to the
   Online Resource, blind-cone figure optional.
