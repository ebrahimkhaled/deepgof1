# W3 — A data-computable bound on the calibration gap

**Status 2026-09-23: revised after two adversarial passes. Theorem G is an asymptotic result under an
explicitly stated Edgeworth-type assumption (G-A); its finite-sample coverage is measured.**

## 1. The objection
Proposition 2 of the manuscript bounds the size of the deployed test at every $n$ by
$\alpha+\varepsilon+\eta$, where $\eta=P(\Delta^+>\varepsilon)$ and
$\Delta^+=\sup_t(F_{\hat\theta}(t)-F_{\theta_0}(t))_+$ is the one-sided gap between the bootstrap law of the
score and its true null law, conditional on the covariates. The CSDA editor's objection: $\varepsilon$ is
unknown. Here it is computed from the data.

## 2. Two constructions, and why the obvious one is useless
**(a) Worst case over a confidence set.** For any set $C$ with $P(\theta_0\in C)\ge1-\gamma$,
$\hat\varepsilon_{wc}=\sup_{\theta\in C}\hat D^+(\hat\theta,\theta)+2\epsilon_{DKW}$ satisfies
$P(\Delta^+>\hat\varepsilon_{wc})\le\gamma+\delta_{MC}$. *Proof:* on $\{\theta_0\in C\}$,
$\Delta^+\le\sup_{\theta\in C}D^+(\hat\theta,\theta)$, and DKW controls every empirical cdf jointly (union
bound) with probability $\ge1-\delta_{MC}$. $\square$ This is the logic of Berger & Boos (1994) and of
the maximised Monte Carlo test the manuscript already built (half the power). **Measured, it is useless:**
the pilot gives $\hat\varepsilon_{wc}\approx.06$–$.21$ against true gaps of .003–.02. The supremum over a
95% set charges every dataset for the worst direction at the edge of the set.

**(b) The gradient-quantile bound (the one proposed).** The gap is driven by $\hat\theta-\theta_0$, which is
approximately $N(0,\mathcal I^{-1})$, and its sensitivity can be estimated from the draws already made at
$\hat\theta$. Conditional on $X$, the likelihood-ratio derivative gives
$$g(t)=\partial_\theta F_\theta(t)\big|_{\hat\theta}=E_{\hat\theta}\bigl[1\{T\le t\}\,U\bigr],\qquad U=X^\top(y-\pi_{\hat\theta}),$$
estimated by averaging over $N$ draws at $\hat\theta$ ($N/n\to\infty$ needed; see G-A) (centred, since $E[U]=0$). Define
$$\hat\varepsilon_q=\text{the }q\text{-quantile of }\sup_{t\in\mathcal T}\bigl(g(t)^\top Z\bigr)_+,\qquad Z\sim N(0,\mathcal I(\hat\theta)^{-1}),$$
over the tail set $\mathcal T=\{t\ge F_{\hat\theta}^{-1}(0.85)\}$. Only the tail enters the size at
$\alpha$; this is the tail-localised form, Proposition S.8.

## 3. Theorem G (revised after adversarial pass 2, 2026-09-23)

**Assumption G-A (a local expansion of the null law, stated as an assumption and not verified).** There
is a function $g_0$ on $\mathcal T$ such that, for every $M<\infty$,
$$\sup_{t\in\mathcal T,\ \|h\|\le M}\Bigl|F^{(n)}_{\theta_0+h/\sqrt n}(t)-F^{(n)}_{\theta_0}(t)-n^{-1/2}g_0(t)^\top h\Bigr|=o(n^{-1/2}),$$
and $n^{1/2}g_n(\hat\theta)\to g_0$ uniformly on $\mathcal T$ in probability, where
$g_n(\theta)=E_\theta[1\{T\le t\}U]$. For each fixed $n$, $\theta\mapsto F^{(n)}_\theta(t)$ is a finite sum
over $y\in\{0,1\}^n$ and is trivially smooth. The content of G-A is the **$n$-uniform** expansion with a
fixed limit, an Edgeworth-type condition. Lattice terms from integer cell counts could violate it at
order $n^{-1/2}$. It is the condition the coverage study tests empirically.

**Also assumed:** A1–A3 and A6; $\alpha\le1-0.85$, so that the rejection region lies in $\mathcal T$;
$q>P\{\sup_t(g_0^\top W)_+=0\}$ (the limit law has an atom at 0); and the $N$ draws used for $\hat g$
satisfy $N/n\to\infty$. The Monte Carlo error of $\hat g^\top Z$ is of order $N^{-1/2}$, while
$\hat\varepsilon$ itself is of order $n^{-1/2}$. **So the deployed $B=199$ draws do NOT suffice.** The
coverage study used $N=20{,}000$. (The earlier claim "at no extra simulation cost" is withdrawn.)

**Theorem G.** Under these conditions:
1. **Asymptotic coverage.** $\liminf_nP(\Delta^+_{\mathcal T}\le\hat\varepsilon_q)\ge q$.
2. **Asymptotic size of the gap-adjusted rule, for the ideal p-value** $p_\infty=1-F_{\hat\theta}(T-)$:
   $\limsup_nP(p_\infty+\hat\varepsilon_q\le\alpha)\le\alpha+(1-q)$. At finite $B$ an additional Monte Carlo
   term of order $B^{-1/2}$ appears, because $\hat\varepsilon$ depends on the data.
3. **Honest size statement for the deployed test.** For every constant $e$,
   $P(p\le\alpha)\le\alpha+e+P(\Delta^+>e)$. This is Proposition 2. Theorem G(1) supplies a data-based
   *estimate* of which $e$ makes $P(\Delta^+>e)$ small. It does not turn the random $\hat\varepsilon$ into a
   finite-sample bound.

*Proof.* (1) Under G-A and $\sqrt n(\hat\theta-\theta_0)=O_P(1)$,
$\sqrt n\,(F_{\hat\theta}-F_{\theta_0})(t)=g_0(t)^\top\sqrt n(\hat\theta-\theta_0)+o_P(1)$ uniformly on
$\mathcal T$. Hence $\sqrt n\,\Delta^+_{\mathcal T}\rightsquigarrow\sup_t(g_0^\top W)_+$ with
$W\sim N(0,\mathcal I_1^{-1})$. By G-A, $N/n\to\infty$ and $\mathcal I(\hat\theta)/n\to\mathcal I_1$,
$\sqrt n\,\hat\varepsilon_q$ converges in probability to the $q$-quantile of that law, which is a
continuity point for $q$ above the atom at 0. The tail set $\mathcal T$ is defined at $\hat\theta$, and its
lower endpoint converges, so both sups are over asymptotically the same set. (2) On
$\{\Delta^+\le\hat\varepsilon\}$, $p_\infty+\hat\varepsilon\le\alpha$ implies $F_{\theta_0}(T-)\ge1-\alpha$.
Asymptotically $T$ and $\hat\theta$ are independent (their covariance is $G-G\mathcal I^{-1}\mathcal I=0$),
so $\limsup P(F_{\theta_0}(T-)\ge1-\alpha)\le\alpha$. Add the $1-q$ from (1). (3) is Proposition 2.
$\square$

**Withdrawn from the first draft** (adversarial pass 2):
* "$P(p\le\alpha)\le\alpha+E[\hat\varepsilon]+\dots$ for every $n$". This is false. Counterexample:
  $\hat\varepsilon=\delta\,1\{V\in(\alpha,\alpha+\delta]\}$ gives $E\hat\varepsilon=\delta^2$ but size
  $\alpha+\delta$.
* "No extra simulation cost."
* "Finite-sample" for the gap-adjusted rule.

**What this is.** A data-computable, asymptotically calibrated confidence bound on the calibration gap,
whose finite-sample coverage is **measured**. It sits beside Proposition 2's exact but
unknown-parameter bound. It is not a finite-sample theorem. The worst-case construction (a) is
finite-sample exact and measured useless. Whether any finite-sample-exact and useful construction exists
is open.

## 3b. Separation and the Firth fallback at n ≤ 50 (explicit, as the W3 contract requires)
The released `deepgof1()` takes the analyst's own fitted `glm` and does **no** Firth fallback. Its
bootstrap refits use plain `glm`, and only a refit that errors is scored $+\infty$. (Corrected
2026-09-23 against `R/deepgof.R`, lines 236–239.) The Firth fallback exists only in the **simulation
harness** (`fit_pol`), which is what the studies below used.
* **Bootstrap failures are inside the bound by construction.** $F_\theta$ is the law of the score as
  deployed, with failed refits at $+\infty$; `scores_at` implements exactly that. So $\Delta^+$,
  $\hat\varepsilon$ and Proposition 2 all refer to the deployed law, failures included. No extra term is
  needed.
* **Observed-data separation adds one term.** Let $E_{\rm sep}$ be the event that the observed MLE does
  not exist. On $E_{\rm sep}$ the statistic is computed at the Firth fit, which is not covered by the
  expansion G-A. For every constant $e$ (Proposition 2 on $E_{\rm sep}^c$, and a trivial bound on
  $E_{\rm sep}$),
  $$P(p\le\alpha)\le\alpha+e+P(\Delta^+>e)+P_{\theta_0}(E_{\rm sep}).$$
  $P_{\theta}(E_{\rm sep})$ is a property of the design and $\theta$ alone, computable by simulating
  outcomes at $\hat\theta$ (no scores needed). It should be reported with the test. In the rare-event
  regime it is not negligible.
* **Measured, n = 50, covariates as in §4** (`fit_batch` convergence, 4,000 datasets each):

  | intercept (event rate) | P(separation) |
  |---|---|
  | −0.5 (.38) | .0000 |
  | −2.0 (.12) | .0025 |
  | −3.0 (.05) | **.055** |

  The coverage study used intercept −0.5, where separation never occurred (0 of 450 observed fits at
  each n; bootstrap failure rate ≤ .0002). So its coverage numbers do not exercise the Firth path. At an
  event rate of 5%, the separation term alone is .055 and must be added. **At n = 50 with rare events,
  the certified size bound is therefore dominated by separation, not by the calibration gap.** The
  honest recommendation is to report $P_{\hat\theta}(E_{\rm sep})$ beside the p-value.

## 4. Results (`01_gap_bound.py`, 2026-09-23; 300 null datasets and 150 alternatives per n; N = 20,000 per cdf)
Design: $x_1,x_2\sim N(0,1)$, $\eta=-0.5+x_1+0.6x_2$; covariates held fixed per dataset; the true gap
is computed against 20,000 draws at $\theta_0$.

| n | true tail gap: median / q90 / max | $\hat\varepsilon_{.95}$ median | **coverage of $\hat\varepsilon_{.95}$** | $\hat\varepsilon_{.99}$ median | coverage of $\hat\varepsilon_{.99}$ |
|---|---|---|---|---|---|
| 50 | .0016 / .035 / .092 | .041 | **.960** (SE .011) | .058 | .990 (SE .006) |
| 100 | .0033 / .018 / .048 | .022 | **.953** (SE .012) | .031 | .997 (SE .003) |
| 200 | .0023 / .011 / .021 | .013 | **.957** (SE .012) | .019 | .993 (SE .005) |

**W3 exit test: met with $\hat\varepsilon_{.99}$** (coverage .990/.997/.993, all clearly above .95); with
$\hat\varepsilon_{.95}$ the coverage is consistent with .95 but not demonstrably above it.
The worst-case (maximised-Monte-Carlo-type) bound is 3–5 times larger (pilot: .06–.21) and is not
competitive.

**What the bound buys in practice.**
* **As a reported certificate (corrected):** a size bound must include the miss probability. With
  $q=.99$: $\alpha+\hat\varepsilon_{.99}+.01\approx$ .118 at n = 50, .091 at n = 100 and .079 at n = 200
  (medians). The first draft reported .09/.07/.063 without the miss term; that was wrong.
  **Coverage of $\hat\varepsilon_{.95}$ (.953–.960, SE about .012) cannot demonstrate $\ge.95$ (lower 95%
  limit about .93); coverage of $\hat\varepsilon_{.99}$ (.990–.997) clearly exceeds .95.** The 'true gap'
  itself carries Monte Carlo sup-noise of about .005–.01 from N = 20,000, larger than its median, so
  coverage is, if anything, slightly understated.
* **As a test** (reject iff $p+\hat\varepsilon_{.99}\le\alpha$): it is conservative at small n. Its power
  under the quadratic alternative C = 1 is .06 against .28 for the deployed test at n = 50, .55 against
  .69 at n = 100, and .93 against .95 at n = 200. So the gap-adjusted test is a reasonable exact-ish
  option from n ≈ 200, and the plug-in bootstrap with the reported certificate is the right default below.
* The deployed test's measured size in this study was .030, .067 and .063 at n = 50, 100 and 200
  (300 datasets each, SE about .013).
