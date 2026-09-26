# W3 upgrade — the rate of the level error of the deployed test

**Status 2026-09-24.** Proved: conditional on the covariates, the size of the deployed test
(parametric bootstrap, $B=199$, $p=(1+\#\{T^*\ge T\})/(B+1)$) differs from the attainable level
$k/(B+1)$ by at most $L_{B,k}\,\mathbb E[\Delta_n\mid X]$ (exact, no assumptions), and
$\mathbb E[\Delta_n\mid X]=O(n^{-1/2}\log n)$ under explicit regularity (Theorem R). Not proved:
$O(1/n)$. What *is* proved about $1/n$: in the Gaussian surrogate of the map the $n^{-1/2}$ term of
the size has coefficient zero (Proposition R.6), because the residual map is exactly uncorrelated with
the score of the fitted model; whether that cancellation survives the lattice and the piecewise-affine
network is the question the rate study measures. The study (RESULT.md) finds the pooled error within
two of its standard errors ($\le .0025$, SE $.0014$) at every $n\ge100$; single-design cells at
$n\ge100$ exceed $\alpha$ by at most $.0038$, and the three cells outside two SE are all *below* $\alpha$.
So the rate is **not identifiable from data of this size**; that is reported as such, not as "$1/n$".

**Revision 2026-09-24 (adversarial review, ADVERSARIAL_REVIEW.md).** (a) The atom deficit is never
zero at finite $n$ (the score takes finitely many values); it is now bounded, $D_{n,k}\le\frac12
L_{B,k}\sum_t P(T=t)^2$. (b) Lemma R.4 uses the activation-region *partition* ($M$ convex pieces,
error multiplier $M$) instead of inclusion–exclusion over a cover ($2^M-1$), and the explicit
constant of Raič (2019). (c) $c_1^2$ corrected from $16$ to $24(d+1)\kappa^2/\lambda^2$ so that the
tail really is $\le n^{-2}$. (d) The unconditional moment claim of Lemma R.1 is withdrawn (the MLE
fails to exist with positive probability). (e) $m_L=B_nW^{-1/2}\varepsilon$, not $B_n\varepsilon$.
(f) The rate is qualitative: $C_X$ contains $M$ and is not usable at any practical $n$.
(g) The radial-(R3) remark is a conjecture, not a proved statement. (h) The ESM's "25.8" for
$L_{199,10}$ (Section 5, after Lemma A.5) is wrong; the value is 26.84.

Notation follows `manuscript/sc_audit/esm/ESM_audit.tex` (Sections 1, 2, 5). All statements are
conditional on the covariate configuration $X=(x_1,\dots,x_n)$, which is what the parametric bootstrap
conditions on.

---

## 1. What Beran's theory predicts, and why it is not the last word here

Beran (1988, JASA, *Prepivoting test statistics*) gives the parametric bootstrap of a root $R_n$ a level
error of order $n^{-1}$ when $R_n$ is asymptotically pivotal and admits a smooth Edgeworth expansion,
and of order $n^{-1/2}$ when its limit law depends on $\theta$. Here the score $T_n=S((m_n-a)/s)$ is
**not** pivotal: its null limit is $\mathcal L(S((Z(\theta)-a)/s))$ with $Z(\theta)\sim N(0,\Sigma(\theta))$
(display (S.eq:nulllimit)), and $\Sigma(\theta)$ moves with $\theta$. On Beran's classification the
deployed one-stage bootstrap is therefore in the $n^{-1/2}$ class, and only a second stage (bootstrapping
the p-value again) would be in the $n^{-1}$ class. Two facts keep this from being the whole story:

1. **Beran's $n^{-1/2}$ term is a covariance, and here it vanishes at first order.** The level of the
   one-stage bootstrap is $\mathbb E[\Psi_k(V_n)]$ with $V_n=1-F_{X,\hat\theta}(T_n^-)$. Expanding in
   $\hat\theta-\theta_0$, the leading term is
   $\mathbb E[\Psi_k'(V^0_n)\,\partial_\theta F_{X,\theta_0}(T_n)^\top(\hat\theta-\theta_0)]$: the product of a
   function of $T_n$ and the estimation error. If $T_n$ and $\hat\theta$ were independent this would be
   $\mathbb E[\cdot]\,\mathbb E[\hat\theta-\theta_0]=O(1/n)$ (the MLE bias). For this map they are
   *exactly uncorrelated at the linear level* (Lemma R.3(c) below: the influence vectors of the map are
   orthogonal to the score of the fitted logistic model, the Chernoff–Lehmann structure), so the
   $n^{-1/2}$ coefficient is zero in the Gaussian surrogate (Proposition R.6). The bootstrap is then
   "second-order correct in the mean even though the root is not pivotal". This is a heuristic for the
   real test and a theorem for the surrogate.
2. **Edgeworth expansions are not available.** The map is a function of the $K^2$ integer cell counts
   (lattice) and of the $d+1$ continuous sufficient statistics of the logistic fit, pushed through a
   piecewise-affine network. No published expansion covers this object, and I do not construct one.
   Whether the lattice terms, which are of order $n^{-1/2}$ in each of $F_{X,\hat\theta}$ and
   $F_{X,\theta_0}$, cancel in their *difference* (the bootstrap reproduces the same lattice at parameters
   $O(n^{-1/2})$ apart) is plausible and unproved.

So the honest target is a rigorous first-order rate with explicit dependence, plus a precise statement
of where the second-order claim would come from.

---

## 2. Set-up and the exact reduction (no assumptions beyond $H_0$)

Conditional on $X$, for a parameter $\theta$ let $F_{X,\theta}$ be the distribution function of the score
$T(X,y)$ when $y_i\sim\mathrm{Bern}(\pi_\theta(x_i))$ independently. Under $H_0$ with true $\theta_0$,
$T_n\sim F_{X,\theta_0}$ and the bootstrap replicates are i.i.d. $F_{X,\hat\theta_n}$ given the data.
Write $V_n=1-F_{X,\hat\theta_n}(T_n^-)$, $V^0_n=1-F_{X,\theta_0}(T_n^-)$,
$\Delta_n=\sup_t|F_{X,\hat\theta_n}(t)-F_{X,\theta_0}(t)|$, and, for $k\in\{1,\dots,B+1\}$,
$\Psi_k(v)=P(\mathrm{Bin}(B,v)\le k-1)$. Lemma A.3 gives
$P(p_n\le k/(B+1)\mid X)=\mathbb E[\Psi_k(V_n)\mid X]$, and Lemma A.4 gives
$P(V^0_n\le u\mid X)\le u$.

**Lemma R.0 (the Lipschitz constant of $\Psi_k$).** $\Psi_k$ is differentiable with
$-\Psi_k'(v)=\mathrm{Beta}(k,B+1-k)$ density at $v$, so $\Psi_k$ is Lipschitz with constant
$L_{B,k}=\max_v \mathrm{beta}(v;k,B+1-k)$. At $B=199$: $L_{199,10}=26.84$ ($\alpha=.05$) and
$L_{199,20}=19.07$ ($\alpha=.10$). *Proof.* $\Psi_k(v)=\int_v^1\mathrm{beta}(u;k,B+1-k)\,du$ is the
standard binomial–beta identity; differentiate. $\square$

**Proposition R.1 (exact reduction to the mean gap).** Under $H_0$, for every $n$, $B$, $k$,
$$
P\bigl(p_n\le\tfrac{k}{B+1}\bigm|X\bigr)\;\le\;\tfrac{k}{B+1}+L_{B,k}\,\mathbb E[\Delta_n\mid X],
\qquad
P\bigl(p_n\le\tfrac{k}{B+1}\bigm|X\bigr)\;\ge\;\mathbb E[\Psi_k(V^0_n)\mid X]-L_{B,k}\,\mathbb E[\Delta_n\mid X],
$$
and $\mathbb E[\Psi_k(V^0_n)\mid X]=k/(B+1)$ when $F_{X,\theta_0}$ is atomless (in general it lies in
$[k/(B+1)-\,\text{atom deficit},\,k/(B+1)]$). Hence
$$
\Bigl|P\bigl(p_n\le\tfrac{k}{B+1}\bigm|X\bigr)-\tfrac{k}{B+1}\Bigr|\;\le\;L_{B,k}\,\mathbb E[\Delta_n\mid X]
\;+\;\Bigl(\tfrac{k}{B+1}-\mathbb E[\Psi_k(V^0_n)\mid X]\Bigr),
$$
the second term (the atom deficit $D_{n,k}$) being nonnegative and bounded by
$$0\le D_{n,k}\le\tfrac12L_{B,k}\sum_t P(T_n=t\mid X)^2\;\le\;\tfrac12L_{B,k}\max_tP(T_n=t\mid X).$$
At a fixed $X$ the score takes finitely many values (at most $2^n$), so $F_{X,\theta_0}$ is never
atomless and $D_{n,k}>0$ at every $n$; the bound says it is controlled by the collision probability
of the null law of the score. (For continuous covariates distinct response vectors generically give
distinct scores, so this is of the order of $\max_yP_{\theta_0}(y)$, exponentially small, unless the
network collapses many response vectors onto one value; that genericity is not proved here.)

*Proof.* $|V_n-V^0_n|=|F_{X,\hat\theta}(T_n^-)-F_{X,\theta_0}(T_n^-)|\le\Delta_n$ pointwise (left limits
of a difference are limits of the difference), so $|\Psi_k(V_n)-\Psi_k(V^0_n)|\le L_{B,k}\Delta_n$ and
$|\mathbb E\Psi_k(V_n)-\mathbb E\Psi_k(V^0_n)|\le L_{B,k}\mathbb E\Delta_n$. By Lemma A.4 $V^0_n$ is
super-uniform, so by Lemma A.5 with $\varepsilon=0$, $\mathbb E\Psi_k(V^0_n)\le\int_0^1\Psi_k=k/(B+1)$.
Deficit bound: let $T'$ be an independent copy of $T_n$ given $X$ and $W\sim U(0,1)$ independent; the
randomized tail $\tilde V=P(T'>T_n\mid T_n)+W\,P(T'=T_n\mid T_n)$ is exactly $U(0,1)$, and
$V^0_n-\tilde V=(1-W)P(T'=T_n\mid T_n)\ge0$. As $\Psi_k$ is nonincreasing and $L_{B,k}$-Lipschitz,
$0\le\Psi_k(\tilde V)-\Psi_k(V^0_n)\le L_{B,k}(1-W)P(T'=T_n\mid T_n)$; take expectations
($\mathbb E\Psi_k(\tilde V)=k/(B+1)$, $\mathbb E(1-W)=\frac12$, $\mathbb E P(T'=T_n\mid T_n)=\sum_tP(T_n=t)^2$).
$\square$

*The $+\infty$ convention.* If a bootstrap refit fails and is scored $+\infty$ (ESM §1, step 5), the
reference law is stochastically larger; the upper half of R.1 holds verbatim (ESM Remark
rem:fsassumes) and the lower half loses at most $L_{B,k}\,P^*(\text{refit fails})$. In the rate study
there were no failed refits in 114,000 replicates.

*Remark.* The ESM (after Lemma A.5) observes that a bound $k/(B+1)+C\,\mathbb E\Delta_n$ with an
**absolute** $C$ is false. Proposition R.1 does not contradict it: the constant is $L_{B,k}\asymp\sqrt B$,
finite at the deployed $B=199$ (26.84 at $k=10$; the ESM's "25.8" at that spot is a slip and must
be corrected to 26.8). Proposition S.7 avoids the $\sqrt B$ by using a quantile of the gap; R.1
pays it and gets a bound in the mean, which is what a rate statement needs. Both are exact.

**Consequence.** The rate of the size error is at most the rate of $\mathbb E[\Delta_n\mid X]$, with the
explicit multiplier $L_{B,k}$. Everything below is about $\mathbb E[\Delta_n\mid X]$.

---

## 3. Regularity conditions for the rate (conditional on $X$)

Throughout, $\tilde x_i=(1,x_i^\top)^\top$, $\pi_i(\theta)=\Lambda(\tilde x_i^\top\theta)$,
$w_i(\theta)=\pi_i(1-\pi_i)$, $\varepsilon_i=y_i-\pi_i(\theta)$, $U(\theta)=\sum_i\tilde x_i\varepsilon_i$ (score),
$I_n(\theta)=\sum_i w_i(\theta)\tilde x_i\tilde x_i^\top$. The cells $C_{k\ell}$ are the rank cells of the
selected covariate pair $(a,b)$ of $X$; given $X$ and the pair they are fixed index sets. Write
$v_C(\theta)=\sum_{i\in C}w_i(\theta)$ and $a_C(\theta)=\sum_{i\in C}w_i(\theta)\tilde x_i$. Let
$N=\{\theta:\|\theta-\theta_0\|\le r_0\}$.

- **(R1) Bounded, well-conditioned design.** $\max_i\|\tilde x_i\|\le\kappa$;
  $\lambda_{\min}(I_n(\theta))\ge\lambda n$ and $v_C(\theta)\ge v_0 n$ for every cell $C$ and every
  $\theta\in N$; and the covariance $\Sigma_n(\theta)$ of the linearized map (Lemma R.3) has rank
  $K^2-1$ with its nonzero eigenvalues in $[\lambda_m,\lambda_M]$ for $\theta\in N$. (The rank deficiency
  is exact: the $v_C^{1/2}$-weighted sum of the rows is identically zero, since $\sum_C a_C=I_ne_0$.)
- **(R2) The network.** $G=S((\cdot-a)/s)$ is $L$-Lipschitz, and $\mathbb R^{K^2}$ is *partitioned* into at
  most $M$ convex sets on each of which $G$ is affine. (True of the shipped architecture: convolutions,
  ReLU, max-pooling, max/mean readouts and a linear head. Piecewise affineness is Lemma S.1(a)
  (not (b), which is the offset-bias bound); Lipschitz continuity is layerwise. The partition is the set of
  activation regions — fixed on/off pattern of every rectifier, fixed arg-max of every max unit with
  ties broken by index — each an intersection of finitely many open or closed half-spaces, hence convex.
  $M\le 2^{\#\text{ReLU}}\prod(\text{pool size})^{\#\text{max units}}$; it is finite and depends only on the
  frozen weights, and for the shipped network (over 2,000 rectifiers) it is astronomically large.)
- **(R3) Quantitative anti-concentration of the Gaussian surrogate.** With
  $Z_{n,\theta}\sim N(0,\Sigma_n(\theta))$: $P(t<G(Z_{n,\theta})\le t+\delta)\le c_A\delta$ for all $t$,
  $\delta>0$, $\theta\in N$, $n\ge n_0$. This is a uniform, quantitative form of A6. Like A6 it is a
  property of the shipped weights, **assumed, not proved, and not yet checked numerically**; it excludes
  flat regions of $G$ carrying Gaussian mass. It *is* checkable by simulation at a given $(X,\hat\theta)$:
  draw $Z\sim N(0,\Sigma_n(\hat\theta))$, score it, and estimate $\sup_t$ of the density of $G(Z)$ and any
  atom (extending the 199,999-distinct-values check of A6 in ESM §8). Such a check estimates $c_A$; it
  does not prove (R3) uniformly in $\theta\in N$ and $n$.
- **(R4) Axis-selection margin.** With $c_j(\theta)=|\theta_j|\hat\sigma_j$, the top-two set of
  $c(\theta_0)$ is $\{a,b\}$ with margin $\gamma=\min_{j\in\{a,b\},\,l\notin\{a,b\}}(c_j-c_l)>0$; and
  $r_0\le\gamma/(4\max_j\hat\sigma_j)$. (Vacuous when $d=2$: the pair is then always $\{1,2\}$. When
  $d\ge3$ it depends on the unknown $\theta_0$, so it is checkable in a simulation, not by the analyst.)

(R1) and (R4) are conditions on $X$ and $\theta_0$; (R2) is a fact; (R3) is the one assumption about the
network. Constants $C$ below depend only on $(\kappa,\lambda,v_0,\lambda_m,\lambda_M,K,d,L,M,c_A,r_0)$
and never on $n$ or $t$.

---

## 4. Four lemmas

**Lemma R.1 (MLE deviation; exponential tail).** Under (R1), for $\theta\in N$ with
$\|\theta-\theta_0\|\le r_0/2$, $y\sim P_\theta$, and every $r\le r_0/2$,
$$P_\theta\bigl(\|\hat\theta-\theta\|>r\bigm|X\bigr)\;\le\;2(d+1)\exp\Bigl(-\frac{n\lambda^2r^2}{8(d+1)\kappa^2}\Bigr).$$
Consequently $\mathbb E_\theta\bigl[\|\hat\theta-\theta\|^q\,\mathbf 1\{\|\hat\theta-\theta\|\le r_0/2\}\bigr]=O(n^{-q/2})$
for every $q\ge1$. (Unconditional moments are infinite: the MLE fails to exist, under separation, with
positive probability at every finite $n$. Nothing below uses them.)

*Proof.* The log-likelihood $\ell$ is concave with $-\nabla^2\ell(\vartheta)=I_n(\vartheta)\succeq\lambda nI$ on
$N$. For $\|\vartheta-\theta\|=r$ (so $\vartheta\in N$), Taylor's theorem with integral remainder gives
$\ell(\vartheta)-\ell(\theta)\le\|U(\theta)\|r-\tfrac12\lambda nr^2$. On the event
$\{\|U(\theta)\|<\tfrac12\lambda nr\}$ this is $<0$ on the whole sphere, so by concavity the maximizer
lies inside the ball: $\|\hat\theta-\theta\|<r$. Each coordinate $U_j=\sum_i\tilde x_{ij}\varepsilon_i$ is a sum
of independent centred terms bounded by $\kappa$, so Hoeffding gives
$P(|U_j|\ge s)\le2\exp(-s^2/(2n\kappa^2))$; with $\|U\|\le\sqrt{d+1}\max_j|U_j|$ and $s=\lambda nr/(2\sqrt{d+1})$
the display follows. The truncated moment bound is the tail integral
$\mathbb E[W^q\mathbf 1\{W\le r_0/2\}]\le\int_0^{r_0/2}qr^{q-1}P(W>r)\,dr$ with $W=\|\hat\theta-\theta\|$. The
event $\{W>r_0/2\}$ (including non-existence of the MLE) has probability $\le2(d+1)e^{-cn}$ and is
handled separately wherever it occurs below. $\square$

**Lemma R.2 (axis selection is stable).** Under (R1), (R4): if $\|\hat\theta-\theta_0\|\le r_0$ then the
selected pair is $(a,b)$. Hence, for $y\sim P_\theta$ with $\theta\in N$ and $\|\theta-\theta_0\|\le r_0/2$,
$P_\theta(\text{pair}\ne(a,b)\mid X)\le2(d+1)\exp(-n\lambda^2r_0^2/(32(d+1)\kappa^2))$.

*Proof.* $|c_j(\hat\theta)-c_j(\theta_0)|\le\hat\sigma_j\|\hat\theta-\theta_0\|\le\gamma/4$ for every $j$, so
every gap of size $\ge\gamma$ in $c(\theta_0)$ survives with size $\ge\gamma/2$ and the top-two set is
unchanged. Apply Lemma R.1 with $r=r_0/2$ (then $\|\hat\theta-\theta_0\|\le\|\hat\theta-\theta\|+\|\theta-\theta_0\|\le r_0$). $\square$

**Lemma R.3 (linearization of the map, and exact orthogonality).** Fix $\theta\in N$ and $y\sim P_\theta$.
On the event $E=\{\text{pair}=(a,b)\}\cap\{\|\hat\theta-\theta\|\le r_0/2\}$ define the *linear map*
$$m_{L,C}\;=\;v_C(\theta)^{-1/2}\sum_{i=1}^n\bigl(\mathbf 1\{i\in C\}-a_C(\theta)^\top I_n(\theta)^{-1}\tilde x_i\bigr)\,\varepsilon_i
\;=:\;\sum_i\xi_{i,C},\qquad C=1,\dots,K^2 .$$
Then:
(a) $m_L=\sum_i\xi_i$ is a sum of independent centred vectors with $\|\xi_i\|\le C_\xi n^{-1/2}$,
$C_\xi=K(1+\kappa^2/\lambda)v_0^{-1/2}$;
(b) $\Sigma_n(\theta)=\mathrm{Cov}(m_L)$ is the matrix of (R1), and $\Sigma_n(\theta)$ is Lipschitz in
$\theta$ on $N$ in Frobenius norm with an $n$-free constant, as is the factor $B_n(\theta)$ with
$m_L=B_n(\theta)W_\theta^{-1/2}\varepsilon$, $W_\theta=\mathrm{diag}(w_i(\theta))$, so that $\Sigma_n=B_nB_n^\top$
($B_n$ is $K^2\times n$, $\|B_n(\theta)-B_n(\theta')\|_F\le C_B\|\theta-\theta'\|$);
(c) $\mathrm{Cov}(m_L,U(\theta))=0$ exactly;
(d) the deployed map satisfies $m=m_L+R$ with
$\|R\|\le C_R\bigl(\sqrt n\,\|\hat\theta-\theta\|^2+\|\hat\theta-\theta\|\,(1+\|m_L\|)\bigr)$ on $E$.

*Proof.* (a) $|\mathbf 1\{i\in C\}-a_C^\top I_n^{-1}\tilde x_i|\le1+\|a_C\|\,\|I_n^{-1}\|\,\kappa\le1+(n\kappa/4)(\lambda n)^{-1}\kappa$,
and $v_C\ge v_0n$; sum over $K^2$ cells. (b) Each entry of $B_n(\theta)$ is
$v_C^{-1/2}(\mathbf 1_C(i)-a_C^\top I_n^{-1}\tilde x_i)w_i^{1/2}$; $w_i,a_C,I_n,v_C$ are smooth in $\theta$
with derivatives bounded by $\kappa\times$(their size) ($|\Lambda'|\le1/4$, $|\Lambda''|\le1/(6\sqrt3)$),
$\|I_n^{-1}\|\le(\lambda n)^{-1}$, so each entry is Lipschitz with constant $O(n^{-1/2})$; there are
$K^2n$ entries, giving $\|B_n(\theta)-B_n(\theta')\|_F^2\le K^2n\cdot O(n^{-1})\|\theta-\theta'\|^2$. Then
$\Sigma_n=B_nB_n^\top$ inherits Lipschitz continuity since $\|B_n\|_F^2=\mathrm{tr}\Sigma_n\le K^2\cdot O(1)$.
(c) $\mathrm{Cov}(m_{L,C},U)=v_C^{-1/2}\sum_i(\mathbf 1_C(i)-a_C^\top I_n^{-1}\tilde x_i)w_i\tilde x_i^\top
=v_C^{-1/2}(a_C^\top-a_C^\top I_n^{-1}I_n)=0$.
(d) On $E$ the cells are the same index sets for the deployed map and for $m_L$. Write $h=\hat\theta-\theta$.
The score equation $U(\hat\theta)=0$ with $U(\hat\theta)=U(\theta)-I_n(\bar\theta)h$ (mean value, $\bar\theta$ on the
segment) and $\|I_n(\bar\theta)-I_n(\theta)\|\le n\kappa^3\|h\|/(6\sqrt3)$ give
$h=I_n(\theta)^{-1}U(\theta)+\rho_1$, $\|\rho_1\|\le C\|h\|^2$. The numerator of the deployed cell entry is
$\sum_{i\in C}(y_i-\pi_i(\hat\theta))=\sum_{i\in C}\varepsilon_i-a_C^\top h-\rho_{2,C}$ with
$|\rho_{2,C}|\le n\kappa^2\|h\|^2/(12\sqrt3)$ (second-order Taylor of $\Lambda$). Substituting $h$:
$\sum_{i\in C}\varepsilon_i-a_C^\top I_n^{-1}U-a_C^\top\rho_1-\rho_{2,C}=v_C^{1/2}m_{L,C}+O(n\|h\|^2)$. The
denominator: $v_C(\hat\theta)^{1/2}=v_C(\theta)^{1/2}(1+O(\kappa\|h\|))$. Dividing,
$m_C=m_{L,C}(1+O(\|h\|))+O(n\|h\|^2/\sqrt{v_0n})$, which is (d). $\square$

**Lemma R.4 (Gaussian approximation of the linear map on the level sets of $G$).** Under (R1)–(R2),
for $\theta\in N$ and $n\ge n_0$,
$$\sup_t\bigl|P(G(m_L)\le t)-P(G(Z_{n,\theta})\le t)\bigr|\;\le\;C_{BE}\,n^{-1/2},
\qquad C_{BE}=M\,c_K\,C_\xi^3\lambda_m^{-3/2},\quad c_K=42(K^2-1)^{1/4}+16\ (\approx118\text{ at }K=6),$$
the order $(K^2-1)^{1/4}$ being Bentkus (2005) and the explicit constant Raič (2019, Theorem 1.1).

*Proof.* By (R2) $\mathbb R^{K^2}$ is partitioned into $M$ convex regions $Q_1,\dots,Q_M$ with $G=\ell_j$
(affine) on $Q_j$. Then $\{G\le t\}$ is the *disjoint* union of the convex sets $A_j=Q_j\cap\{\ell_j\le t\}$,
so $|P(S\in\{G\le t\})-P(Z\in\{G\le t\})|\le\sum_j|P(S\in A_j)-P(Z\in A_j)|$, $M$ terms. (The first
draft used a closed-polyhedra *cover* and inclusion–exclusion, which costs $2^M-1$ terms; the
partition is sharper and equally rigorous.) Both $m_L$ and $Z_{n,\theta}$ live on the
$(K^2-1)$-dimensional range $\mathcal R$ of $\Sigma_n(\theta)$ (the same subspace, since $m_L=B_n\varepsilon$ and
$Z=B_n\xi$ have the same column space); intersecting with $\mathcal R$ keeps convexity. Standardize by
$\Sigma_n^{-1/2}$ on $\mathcal R$: the summands $\Sigma_n^{-1/2}\xi_i$ are independent, centred, with
$\sum_i\mathrm{Cov}=I_{K^2-1}$ and $\sum_i\mathbb E\|\Sigma_n^{-1/2}\xi_i\|^3\le n\,(\lambda_m^{-1/2}C_\xi n^{-1/2})^3
=\lambda_m^{-3/2}C_\xi^3n^{-1/2}$. Bentkus (2005, *Theory Probab. Appl.* 49, 311–323; main result — the
theorem number is not verified against a copy of the paper): for independent centred
$X_i\in\mathbb R^{d'}$ with $\mathrm{Cov}(\sum X_i)=I$,
$\sup_{C\text{ convex}}|P(\sum_iX_i\in C)-P(N(0,I)\in C)|\le c\,d'^{1/4}\sum_i\mathbb E\|X_i\|^3$, with an
unspecified absolute $c$; Raič (2019, *Bernoulli* 25(4A), Theorem 1.1, verified) gives
$c\,d'^{1/4}\le42d'^{1/4}+16$ for independent, not necessarily identically distributed summands and all
measurable convex sets. The bound holds for lattice summands (no Cramér condition). Apply it to each
of the $M$ convex sets and sum. $\square$

*Remark (the rate is qualitative).* $M$ is astronomically large for the shipped network (over 2,000
rectifiers; even the region-counting bounds for a 36-dimensional input are of order $10^{50}$ or
more), so $C_{BE}n^{-1/2}<1$ would need $n$ far beyond any sample size. The lemma, and Theorem R, are
statements about the **order in $n$** with a constant that is **finite and free of $n$ and $t$** but
has no usable value. The paper must say so, and does (INSERT_ESM, INSERT_MAIN). A route with a
dimension-only constant exists (smooth $G$ at scale $\epsilon$, Lindeberg with third derivatives of
order $\epsilon^{-2}$, then (R3)) but gives only $n^{-1/6}$ up to logarithms.

---

## 5. The theorem

**Theorem R (rate of the calibration gap and of the level error).** Assume $H_0$ with $\theta_0$, and
(R1)–(R4) conditional on $X$. Then there are $n_0$ and $C_X<\infty$, depending only on the constants
listed in Section 3, such that for all $n\ge n_0$ and every configuration $X$ of size $n$ satisfying
(R1)–(R4) with those constants (so the statement is non-asymptotic in form; "$n\to\infty$" means a
sequence of designs with uniform constants)
$$\mathbb E[\Delta_n\mid X]\;\le\;C_X\,\frac{\log n}{\sqrt n},$$
and therefore, by Proposition R.1, for every $k\in\{1,\dots,B+1\}$,
$$\Bigl|P\bigl(p_n\le\tfrac{k}{B+1}\bigm|X\bigr)-\tfrac{k}{B+1}\Bigr|\;\le\;L_{B,k}\,C_X\,\frac{\log n}{\sqrt n}
\;+\;\Bigl(\tfrac{k}{B+1}-\mathbb E[\Psi_k(V^0_n)\mid X]\Bigr),$$
the last term being the atom deficit $D_{n,k}$ of Proposition R.1 (nonnegative, positive at every finite
$n$, at most $\frac12L_{B,k}\sum_tP(T_n=t)^2$, and a one-sided, conservative term). At the deployed $B=199$ the multipliers
are $L_{199,10}=26.84$ at $\alpha=.05$ and $L_{199,20}=19.07$ at $\alpha=.10$.

**Proof.** Fix $\theta\in\{\theta_0,\hat\theta_n\}$; both lie in $N$ with $\|\theta-\theta_0\|\le r_0/2$ on the
event $E_0=\{\|\hat\theta_n-\theta_0\|\le r_0/2\}$, whose complement has probability $\le2(d+1)e^{-c_0n}$ by
Lemma R.1. Off $E_0$ bound $\Delta_n\le1$; the contribution to $\mathbb E\Delta_n$ is exponentially small.
On $E_0$ we bound, for each such $\theta$,
$$\text{(I)}_\theta=\sup_t\bigl|F_{X,\theta}(t)-P(G(Z_{n,\theta})\le t)\bigr|
\quad\text{and}\quad
\text{(II)}=\sup_t\bigl|P(G(Z_{n,\hat\theta})\le t)-P(G(Z_{n,\theta_0})\le t)\bigr|,$$
and use $\Delta_n\le\text{(I)}_{\theta_0}+\text{(I)}_{\hat\theta}+\text{(II)}$. Note that $\hat\theta$ is
data-measurable, so $\text{(I)}_{\hat\theta}$ and (II) are random through $\hat\theta$ only; all
probabilities inside (I) and (II) are over a fresh $y\sim P_\theta$ (the bootstrap draw) or over the
Gaussian surrogate, with $\theta$ held fixed.

*Step 1: (I)$_\theta$, the remainder.* Let $y\sim P_\theta$, $E$ as in Lemma R.3 with $r_n:=c_1\sqrt{\log n/n}$
in place of $r_0/2$, where $c_1$ is chosen so that Lemma R.1 gives $P(\|\hat\theta-\theta\|>r_n)\le n^{-2}$ (that is,
$c_1^2=24(d+1)\kappa^2/\lambda^2$, giving $2(d+1)n^{-3}\le n^{-2}$ for $n\ge2(d+1)$; the first draft's
$16(d+1)\kappa^2/\lambda^2$ gave $2(d+1)n^{-2}$, off by the factor $2(d+1)$), and $n_0$ so that $r_n\le r_0/2$
and $n_0\ge2(d+1)$. Then
$P(E^c)\le n^{-2}+2(d+1)e^{-cn}$ by Lemmas R.1–R.2. On $E$, Lemma R.3(d) gives
$\|R\|\le C_R\bigl(c_1^2\log n/\sqrt n+r_n(1+\|m_L\|)\bigr)$. For any $t$,
$$\bigl|\mathbf 1\{G(m_L+R)\le t\}-\mathbf 1\{G(m_L)\le t\}\bigr|\;\le\;\mathbf 1\{|G(m_L)-t|\le L\|R\|\},$$
because $|G(m_L+R)-G(m_L)|\le L\|R\|$. Split on $\{\|m_L\|\le\rho_n\}$ with $\rho_n=c_2\sqrt{\log n}$:
$m_L$ is a sum of independent vectors bounded by $C_\xi n^{-1/2}$, so by Hoeffding coordinatewise
$P(\|m_L\|>\rho_n)\le2K^2\exp(-\rho_n^2/(2K^2C_\xi^2))\le n^{-2}$ for a suitable $c_2$. On
$E\cap\{\|m_L\|\le\rho_n\}$, $L\|R\|\le\delta_n:=LC_R\bigl(c_1^2\log n/\sqrt n+r_n(1+\rho_n)\bigr)
\le C_3\log n/\sqrt n$. Therefore
$$\text{(I)}_\theta\;\le\;\sup_t P\bigl(|G(m_L)-t|\le\delta_n\bigr)+P(E^c)+P(\|m_L\|>\rho_n)
+\sup_t\bigl|P(G(m_L)\le t)-P(G(Z_{n,\theta})\le t)\bigr| .$$
The first term is a difference of two level-set probabilities of $G(m_L)$, so by Lemma R.4 applied twice
and (R3),
$P(|G(m_L)-t|\le\delta_n)\le P(G(m_L)\le t+\delta_n)-P(G(m_L)\le t-2\delta_n)
\le P(t-2\delta_n<G(Z_{n,\theta})\le t+\delta_n)+2C_{BE}n^{-1/2}\le3c_A\delta_n+2C_{BE}n^{-1/2}$
(the set $\{|G-t|\le\delta\}$ is contained in $\{G\le t+\delta\}\setminus\{G\le t-2\delta\}$, a difference
of two level sets, each covered by Lemma R.4). The last term is Lemma R.4. Altogether
$\text{(I)}_\theta\le3c_AC_3\log n/\sqrt n+3C_{BE}n^{-1/2}+2n^{-2}+2(d+1)e^{-cn}$.

*Step 2: (II), the Gaussian comparison.* Let $\xi\sim N(0,I_n)$ (a fresh Gaussian vector, not the summands
of Lemma R.3) and couple $Z_{n,\theta_0}=B_n(\theta_0)\xi$,
$Z_{n,\hat\theta}=B_n(\hat\theta)\xi$; $\eta:=Z_{n,\hat\theta}-Z_{n,\theta_0}=(B_n(\hat\theta)-B_n(\theta_0))\xi$ is
Gaussian with $\mathbb E\|\eta\|^2=\|B_n(\hat\theta)-B_n(\theta_0)\|_F^2\le C_B^2\|\hat\theta-\theta_0\|^2=:\sigma^2$
(Lemma R.3(b)). By the same indicator inequality,
$\text{(II)}\le\sup_tP(|G(Z_{n,\theta_0})-t|\le L\|\eta\|)\le\sup_tP(|G(Z_{n,\theta_0})-t|\le L\delta)+P(\|\eta\|>\delta)$
for any $\delta$. The Gaussian norm satisfies $P(\|\eta\|>\sigma+s)\le\exp(-s^2/(2\sigma^2))$ (Gaussian
concentration of the Lipschitz functional $\|\cdot\|$, the operator norm of the covariance being $\le\sigma^2$).
Take $\delta=\sigma(1+2\sqrt{\log n})$; then $P(\|\eta\|>\delta)\le n^{-2}$ and, by (R3),
$\text{(II)}\le2c_AL\sigma(1+2\sqrt{\log n})+n^{-2}$. On $E_0\cap\{\|\hat\theta_n-\theta_0\|\le r_n\}$ (probability
$\ge1-n^{-2}-2(d+1)e^{-cn}$) $\sigma\le C_Bc_1\sqrt{\log n/n}$, so $\text{(II)}\le C_4\log n/\sqrt n+n^{-2}$.

*Step 3: assemble.* On the good event, $\Delta_n\le2\max_\theta\text{(I)}_\theta+\text{(II)}\le C_X'\log n/\sqrt n$;
off it, $\Delta_n\le1$ with probability $\le2n^{-2}+4(d+1)e^{-cn}$. Take expectations. The second display
is Proposition R.1. $\square$

**What the logarithm is.** It enters at two places only: decoupling the remainder $R$ from $m_L$ (Step 1)
and decoupling $\eta$ from $Z$ (Step 2), each by a union bound between an anti-concentration event and a
sub-Gaussian tail. *Conjecture (not proved):* under a stronger radial form of (R3),
$P\bigl(|G(Z_{n,\theta})-t|\le\delta(1+\|Z_{n,\theta}\|)\bigr)\le c_A'\delta$, the union bounds become
unnecessary and the rate is $C_Xn^{-1/2}$. This needs (R3)-type control for the non-Gaussian $m_L$ at
scale $n^{-1/2}$ (itself only available up to $C_{BE}n^{-1/2}$) and truncated moments of
$\hat\theta-\theta$ (Lemma R.1); I have not written it out, so the paper states it as expected, not proved.

**What is not claimed.** (i) $O(n^{-1})$. (ii) Any Edgeworth expansion. (iii) A usable numerical value of
$C_X$: it contains $M$ (the number of activation regions) through Lemma R.4; the rate is qualitative.
(iv) Anything about the atom deficit beyond its sign and the collision bound of Proposition R.1.

---

## 6. Why $O(1/n)$ is plausible: the Gaussian surrogate is second-order correct

Replace the map by its Gaussian surrogate and the estimator by its Gaussian limit: $Z\sim N(0,\Sigma_n(\theta_0))$
and $\hat\theta\sim N(\theta_0,I_n(\theta_0)^{-1})$, **independent** (the surrogate of Lemma R.3(c)).
Write $\Phi_\theta(t)=P(G(Z_{n,\theta})\le t)$ and $c_\alpha(\theta)=\Phi_\theta^{-1}(1-\alpha)$. The surrogate
test rejects when $G(Z)>c_\alpha(\hat\theta)$ (the infinite-$B$ bootstrap), so its size is
$s_n=\mathbb E[\phi(\hat\theta)]$ with $\phi(\theta)=1-\Phi_{\theta_0}(c_\alpha(\theta))$ and $\phi(\theta_0)=\alpha$.

**Proposition R.6.** Assume (R1)–(R3), that $\Phi_\theta(t)$ is $C^2$ in $(t,\theta)$ near
$(c_\alpha(\theta_0),\theta_0)$ with a density $\Phi_{\theta_0}'(c_\alpha(\theta_0))>0$, uniformly in $n$. Then
$$s_n=\alpha+\tfrac12\,\mathrm{tr}\bigl(\nabla^2\phi(\theta_0)\,I_n(\theta_0)^{-1}\bigr)+o(n^{-1})=\alpha+O(n^{-1}).$$
The coefficient of $n^{-1/2}$ is zero.

*Proof.* By the implicit function theorem $c_\alpha(\theta)$ is $C^2$ near $\theta_0$, hence so is $\phi$.
Taylor: $\phi(\hat\theta)=\alpha+\nabla\phi(\theta_0)^\top(\hat\theta-\theta_0)+\tfrac12(\hat\theta-\theta_0)^\top
\nabla^2\phi(\theta_0)(\hat\theta-\theta_0)+o(\|\hat\theta-\theta_0\|^2)$. Take expectations: the linear term has
mean zero because $\mathbb E\hat\theta=\theta_0$ (in the surrogate) and because $\hat\theta$ is independent of
$Z$ so that $\nabla\phi(\theta_0)$ is a constant; the quadratic term has mean
$\tfrac12\mathrm{tr}(\nabla^2\phi\,I_n^{-1})=O(1/n)$ by (R1). $\square$

*What this shows and does not show.* In the real test, (a) $\mathbb E(\hat\theta-\theta_0)=O(1/n)$ rather than
$0$ (MLE bias), which contributes $\nabla\phi^\top\mathbb E(\hat\theta-\theta_0)=O(1/n)$, harmless; (b) $T_n$ and
$\hat\theta$ are uncorrelated at the linear level but not independent: their dependence is carried by the
$n^{-1/2}$ non-Gaussian terms, and the resulting covariance contribution
$\mathrm{Cov}(\Psi_k'(V^0_n)\partial_\theta F(T_n),\hat\theta-\theta_0)$ is heuristically
$O(n^{-1/2})\times O(n^{-1/2})=O(1/n)$; (c) the lattice terms of order $n^{-1/2}$ in $F_{X,\theta}$ are the
same function of the same cell counts at the two parameters and are expected to cancel in the difference
to first order. (a) is proved; (b) and (c) are where an Edgeworth expansion would be needed and where I
stop. **So: the $n^{-1/2}$ term has no Gaussian source; if it is present in the deployed test it is a
lattice or network-nonsmoothness effect.** That is a sharp, falsifiable statement, and the rate study is
its test: an error visibly decaying like $n^{-1/2}$ would refute the cancellation; an error at the noise
floor at every $n$ is consistent with it but does not prove it.

---

## 7. Relation to the existing results

| result | assumptions | statement | rate |
|---|---|---|---|
| Lemma S.2 | pivotality (idealized) | size $\le\lfloor\alpha(B+1)\rfloor/(B+1)$ | exact |
| Prop. S.3 | A1–A3, A6 | size $\to k/(B+1)$ | $o(1)$ |
| Prop. S.7 | $H_0$ only | size $\le k/(B+1)+\varepsilon+P(\Delta^+>\varepsilon)$ | exact, gap unknown |
| Thm S.11 | G-A (Edgeworth-type, unverified) | data-computable $\hat\varepsilon_q$, asymptotic coverage | $\Delta^+=O_P(n^{-1/2})$ *under G-A* |
| **Prop. R.1** | $H_0$ only | $\lvert\text{size}-k/(B+1)\rvert\le L_{B,k}\mathbb E\Delta_n$ (+ atom deficit) | exact |
| **Thm R** | (R1)–(R4), no Edgeworth | $\mathbb E\Delta_n\le C_X\log n/\sqrt n$ | $n^{-1/2}\log n$, proved |
| **Prop. R.6** | Gaussian surrogate | size $=\alpha+O(1/n)$ | $1/n$ in the surrogate only |

Theorem R is the first statement in the paper that gives the level error a rate without assuming an
expansion: G-A of Theorem S.11 *assumes* $F_{\theta_0+h/\sqrt n}-F_{\theta_0}=n^{-1/2}g_0^\top h+o(n^{-1/2})$;
Theorem R *proves* $\sup_t|F_{\hat\theta}-F_{\theta_0}|=O(n^{-1/2}\log n)$ in mean from Berry–Esseen and
Lipschitz continuity alone, at the price of the logarithm and of a non-explicit constant.

## 8. References used
- Beran, R. (1988). Prepivoting test statistics: a bootstrap view of asymptotic refinements. *JASA* 83, 687–697.
- Bentkus, V. (2005). A Lyapunov-type bound in $\mathbb R^d$. *Theory of Probability and Its Applications* 49(2), 311–323. doi:10.1137/S0040585X97981123 (Crossref-verified 2026-09-24; no local copy, theorem number unverified, so it is cited without a locator).
- Hoeffding, W. (1963). Probability inequalities for sums of bounded random variables. *JASA* 58(301), 13–30. doi:10.1080/01621459.1963.10500830 (Crossref-verified 2026-09-24).
- Raič, M. (2019). A multivariate Berry–Esseen theorem with explicit constants. *Bernoulli* 25(4A), 2824–2853. doi:10.3150/18-BEJ1072 (Crossref-verified; Theorem 1.1 read in arXiv:1802.06475v4).
- Chernoff, H. and Lehmann, E. L. (1954). The use of maximum likelihood estimates in $\chi^2$ tests for goodness of fit. *Ann. Math. Statist.* 25, 579–586.
- Moore, D. S. and Spruill, M. C. (1975). Unified large-sample theory of general chi-squared statistics for tests of fit. *Ann. Statist.* 3, 599–616.
(All already in the paper's bibliography except Bentkus 2005, Hoeffding 1963 and Raič 2019, which must be added and verified in the compiled .bbl if INSERT_ESM.tex is used.)
