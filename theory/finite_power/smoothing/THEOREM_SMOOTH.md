# W4 — A tight finite-sample power certificate by Gaussian smoothing

Status 2026-09-24: statement and full proof below; the computation is in this folder
(`smooth_core.py`, `01_…05_*.py`, driver `run_all.py`); the numbers are in `RESULT.md` and
`SMOOTH_CERT.json`. Revised 2026-09-24 after `ADVERSARIAL_REVIEW.md`: $n_{\mathrm{cert}}$ rounded up, the
tail level $c^{\circ}$ made explicit, the global constants and their derivation stated, per-rule vs
joint probability stated.

## 0. What is certified and what is not

Everything below is a statement about the **Gaussian model of the map** that the paper's Proposition
S.12 (probabilistic CROWN) already uses: the map is $m=tv+e$ with $v$ a fixed unit direction in
$\mathbb R^{36}$, $e\sim N(0,I_{36})$, $S$ the frozen network, and a signal is converted to a sample size by
$n=4t^{2}$ (a signal $\|\mu\|=0.5$ per $\sqrt n$). It is **not** a guarantee for the finite-$n$ Bernoulli
law of the deployed test (referee point B2): the deployed map is Gaussian only in the limit, its null
covariance is the $\Sigma$ of Theorem S.9 rather than $I$, and its critical value is a bootstrap rank.
The theorem is written for a general, possibly singular, covariance $\Sigma$ so that it also applies to
the local model of Theorem S.9, where the mean $\delta(h)$ lies in $\mathrm{range}(\Sigma)$; the numbers
are computed for $\Sigma=I$ so that they compare with Proposition S.12 cell by cell.

Two rejection rules are certified in the model:

* **Fixed threshold** (the rule of Proposition S.12 and its table): reject iff $S(m)>\bar t$, with
  $\bar t=2.4865$ a $99.9\%$ upper confidence bound on the $.95$ quantile of $S(e)$ under the null
  (from $400{,}000$ null draws, `03_certify.py`).
* **Bootstrap rule** (the deployed one, $B=199$, $\alpha=.05$): given $m$, draw $T^{*1},\dots,T^{*B}$
  i.i.d. from the null law $F_0$ of $S(e)$, independent of $m$, and reject iff
  $\#\{b:T^{*b}\ge S(m)\}\le 9$. Its power is Dufour's formula
  $\beta_{199,10}(\mu)=E[\mathrm{Bin}(9;199,W)]$, $W=1-F_0(S(m)-)$ (Theorem S.9 item 4 gives this as
  the limit of the deployed test's power).

The certificate is a **Monte Carlo certificate**: it holds with probability at least $1-\eta$ over the
Monte Carlo draws used to build it (treating the pseudo-random generator as ideal), where $\eta$ is fixed
in advance and every piece of randomness is charged to it. Nothing in the certified inequality is a point
estimate. **The probability is per direction and per rule** ($1-\eta-.001=.998$ with the critical-value
event); jointly, by the union bound, the five fixed-rule rows hold with probability
$\ge1-5(.001)-.001=.994$ and all twenty cells (two rules, two targets, five directions; the two targets of
one rule and direction share one Clopper–Pearson event) with probability $\ge1-10(.001)-.001-.001=.988$.

Numbering in the ESM insert: Lemma 1+2 = Lemma A.6 (the sixth auxlemma); the fixed-threshold theorem =
Proposition S.12a; the bootstrap-rule theorem = Corollary S.12b (fixed tags, so S.13–S.19 do not shift).

## 1. The smoothing lemma

Throughout, $\Phi$ is the standard normal distribution function, $\Phi^{-1}(0)=-\infty$,
$\Phi^{-1}(1)=+\infty$, and $\Phi(\pm\infty)=1,0$.

**Lemma 1 (Gaussian shift; identity covariance).** Let $h:\mathbb R^{d}\to[0,1]$ be measurable,
$G\sim N(0,I_{d})$, and $H(\mu)=E\,h(\mu+G)$. Then for all $\mu,\mu'\in\mathbb R^{d}$,
$$\Phi\bigl(\Phi^{-1}(H(\mu))-\|\mu'-\mu\|\bigr)\ \le\ H(\mu')\ \le\ \Phi\bigl(\Phi^{-1}(H(\mu))+\|\mu'-\mu\|\bigr).$$
Equivalently, $\mu\mapsto\Phi^{-1}(H(\mu))$ is $1$-Lipschitz for the Euclidean norm wherever it is finite,
and $H$ is identically $0$ or identically $1$ if it takes that value anywhere.

*Proof.* Fix $\mu\neq\mu'$, put $s=\|\mu'-\mu\|>0$ and $u=(\mu'-\mu)/s$. Write $p=H(\mu)$ and let
$\varphi_{\mu}$ be the $N(\mu,I_d)$ density. The likelihood ratio is
$$\frac{\varphi_{\mu'}(x)}{\varphi_{\mu}(x)}=\exp\bigl(s\,u^{\top}(x-\mu)-s^{2}/2\bigr),$$
a strictly increasing function of $u^{\top}(x-\mu)$, and $u^{\top}(G)\sim N(0,1)$.

*Upper bound.* If $p=1$ there is nothing to prove. Otherwise let
$h^{*}(x)=1\{u^{\top}(x-\mu)\ge z\}$ with $z=\Phi^{-1}(1-p)=-\Phi^{-1}(p)$, so that
$E\,h^{*}(\mu+G)=P\{N(0,1)\ge z\}=p$ (for $p=0$, $z=+\infty$ and $h^{*}\equiv0$). Let
$k=\exp(sz-s^{2}/2)$ be the likelihood ratio on the boundary. On $\{h^{*}=1\}$ the ratio is $\ge k$ and
$h^{*}-h\ge0$; on $\{h^{*}=0\}$ the ratio is $<k$ and $h^{*}-h\le0$. Hence
$(h^{*}-h)(\varphi_{\mu'}-k\varphi_{\mu})\ge0$ pointwise, and integrating,
$$H^{*}(\mu')-H(\mu')\ \ge\ k\bigl(H^{*}(\mu)-H(\mu)\bigr)=k(p-p)=0,$$
where $H^{*}(\mu')=E\,h^{*}(\mu'+G)=P\{u^{\top}(\mu'-\mu)+N(0,1)\ge z\}=P\{N(s,1)\ge z\}
=\Phi(s-z)=\Phi(\Phi^{-1}(p)+s)$. So $H(\mu')\le\Phi(\Phi^{-1}(H(\mu))+s)$. (This is the
Neyman–Pearson lemma for the shift family, proved directly.)

*Lower bound.* Apply the upper bound to $1-h$, whose smoothed value is $1-H$:
$1-H(\mu')\le\Phi(\Phi^{-1}(1-H(\mu))+s)=\Phi(-\Phi^{-1}(H(\mu))+s)=1-\Phi(\Phi^{-1}(H(\mu))-s)$.

*Edge values.* If $H(\mu)=0$ then $h=0$ Lebesgue-almost everywhere (the Gaussian density is positive), so
$H\equiv0$; likewise $H(\mu)=1$ forces $H\equiv1$. With the conventions for $\Phi^{-1}(0),\Phi^{-1}(1)$
the display holds in these cases too. $\square$

The constant $1$ is attained by half-space indicators, so the lemma cannot be improved without using
more about $h$ than measurability.

**Lemma 2 (singular covariance).** Let $\Sigma$ be symmetric positive semidefinite of rank $r$, with
spectral decomposition $\Sigma=V\Lambda V^{\top}$ ($V$ is $d\times r$ with orthonormal columns,
$\Lambda$ diagonal and positive), and let $Z\sim N(0,\Sigma)$. For measurable $h:\mathbb R^{d}\to[0,1]$ put
$H(\mu)=E\,h(\mu+Z)$ and $\|\delta\|_{\Sigma^{+}}=(\delta^{\top}\Sigma^{+}\delta)^{1/2}$ with $\Sigma^{+}$
the Moore–Penrose inverse. If $\mu'-\mu\in\mathrm{range}(\Sigma)$ then
$$\Phi\bigl(\Phi^{-1}(H(\mu))-\|\mu'-\mu\|_{\Sigma^{+}}\bigr)\le H(\mu')\le
\Phi\bigl(\Phi^{-1}(H(\mu))+\|\mu'-\mu\|_{\Sigma^{+}}\bigr).$$

*Proof.* $Z\overset{d}{=}V\Lambda^{1/2}G$ with $G\sim N(0,I_{r})$. Since
$\mathrm{range}(V\Lambda^{1/2})=\mathrm{range}(V)=\mathrm{range}(\Sigma)$, there is $w\in\mathbb R^{r}$
with $\mu'-\mu=V\Lambda^{1/2}w$, namely $w=\Lambda^{-1/2}V^{\top}(\mu'-\mu)$, and
$\|w\|^{2}=(\mu'-\mu)^{\top}V\Lambda^{-1}V^{\top}(\mu'-\mu)=\|\mu'-\mu\|_{\Sigma^{+}}^{2}$. Define
$\tilde h(g)=h(\mu+V\Lambda^{1/2}g)$ on $\mathbb R^{r}$. Then $H(\mu)=E\,\tilde h(G)$ and
$H(\mu')=E\,h(\mu+V\Lambda^{1/2}(w+G))=E\,\tilde h(w+G)$. Lemma 1 on $\mathbb R^{r}$ gives the claim. $\square$

If $\mu'-\mu\notin\mathrm{range}(\Sigma)$ the two laws are mutually singular and no such bound can hold;
in the local model of Theorem S.9 this case does not arise, because $\delta(h)\in\mathrm{range}(\Sigma)$ for
every $h$ (Step 7 of its proof), and along a ray $\mu=\lambda u$ with $u\in\mathrm{range}(\Sigma)$ the
distance between $\lambda u$ and $\lambda' u$ is $|\lambda-\lambda'|\,\|u\|_{\Sigma^{+}}$. For the model
computed here, $\Sigma=I$ and $\|v\|_{\Sigma^{+}}=\|v\|=1$.

**Corollary (power along a ray).** For a threshold $c$ put $q_{c}(t)=P\{S(tv+e)>c\}$,
$e\sim N(0,I_{36})$. Then $t\mapsto\Phi^{-1}(q_{c}(t))$ is $1$-Lipschitz: for all $t,t'$,
$q_{c}(t')\ge\Phi(\Phi^{-1}(q_{c}(t))-|t'-t|)$. (Lemma 1 with $h=1\{S>c\}$.)

## 2. The critical value and the bootstrap layer

**Lemma 3 (an upper confidence bound on a null quantile).** Let $S^{*}_{1},\dots,S^{*}_{N_{0}}$ be
i.i.d. from $F_{0}$, $S^{*}_{(k)}$ the $k$-th order statistic, $0<p<1$, and
$c_{p}=\inf\{x:F_{0}(x)\ge p\}$. Let $k$ be the smallest integer with $P\{\mathrm{Bin}(N_{0},p)\ge k\}\le a$.
Then $P\{S^{*}_{(k)}<c_{p}\}\le a$.

*Proof.* $S^{*}_{(k)}<c_{p}$ iff at least $k$ draws are $<c_{p}$. The number of such draws is
$\mathrm{Bin}(N_{0},F_{0}(c_{p}-))$ and $F_{0}(c_{p}-)\le p$ by the definition of $c_{p}$. A binomial with
a smaller success probability is stochastically smaller, so
$P\{\mathrm{Bin}(N_{0},F_{0}(c_{p}-))\ge k\}\le P\{\mathrm{Bin}(N_{0},p)\ge k\}\le a$. $\square$

The paper's $\bar t$ is exactly this bound with $p=.95$, $a=.001$, $N_{0}=400{,}000$. For the
bootstrap rule the levels $w_{1}>\dots>w_{m}$ (here $m=14$, from $.05$ down to $.001$) get bounds
$\hat c_{l}=S^{*}_{(k_{l})}\ge c_{1-w_{l}}$ each at level $a=\eta_{c}/m$ from $N_{0}=2{,}000{,}000$ fresh
null draws (`01_null_quantiles.py`); all $m$ hold together with probability $\ge1-\eta_{c}$.

**Lemma 4 (layer cake for the bootstrap rule).** Let $g(w)=P\{\mathrm{Bin}(B,w)\le k-1\}$, which is
non-increasing in $w$ with $g(1)=0$ (for $k\le B$). Let $1=w_{0}>w_{1}>\dots>w_{m}>0$ and
$\Delta g_{l}=g(w_{l})-g(w_{l-1})\ge0$. For any random $W\in[0,1]$,
$$E\,g(W)\ \ge\ \sum_{l=1}^{m}\Delta g_{l}\,P\{W\le w_{l}\}.$$
If $W=1-F_{0}(S(m)-)$ and $\hat c_{l}\ge c_{1-w_{l}}$, then $\{S(m)>\hat c_{l}\}\subseteq\{W\le w_{l}\}$, so
$$\beta_{B,k}(\mu)=E\,g(W)\ \ge\ \sum_{l=1}^{m}\Delta g_{l}\,P\{S(m)>\hat c_{l}\}.$$

*Proof.* Pointwise, $g(W)\ge\sum_{l}\Delta g_{l}1\{W\le w_{l}\}$: if $W>w_{1}$ the right side is $0\le g(W)$;
if $w_{j+1}<W\le w_{j}$ for some $1\le j\le m$ (with $w_{m+1}=0$) the right side telescopes to
$g(w_{j})-g(w_{0})=g(w_{j})\le g(W)$ because $g$ is non-increasing and $W\le w_{j}$. Take expectations. For the
inclusion: $S(m)>\hat c_{l}\ge c_{1-w_{l}}$ gives $F_{0}(S(m)-)\ge F_{0}(c_{1-w_{l}})\ge1-w_{l}$, i.e.
$W\le w_{l}$. No continuity of $F_{0}$ is used. $\square$

The weights used are $\Delta g=(.461,.132,.131,.113,.084,.050,.022,.006,.001,\dots)$ at
$w=(.05,.045,\dots,.001)$; their sum is $g(.001)=1-4\times10^{-14}$. So the layer cake loses nothing at
large signal, and at moderate signal it loses only the mass of $W$ strictly between consecutive levels.

## 3. The certificate

Fix a direction $v$, a target power $\beta\in(0,1)$, a Monte Carlo budget $\eta$, a tail probability
$\varepsilon\le1-\beta$ with $R=R_{\varepsilon}$ such that $P\{\|e\|>R\}=\varepsilon$ ($e\sim N(0,I_{36})$;
$R=8.74$ for $\varepsilon=10^{-4}$), and a **deterministic** net $t_{1}<\dots<t_{J}$ with draw counts
$M_{1},\dots,M_{J}$, all chosen before any draw.

*Monte Carlo at the net.* At each $t_{j}$ draw $e_{j1},\dots,e_{jM_{j}}$ i.i.d. $N(0,I_{36})$, independent
across $j$, and let $X_{j}=\#\{i:S(t_{j}v+e_{ji})>\bar t\}$. Let $q^{\mathrm{lo}}_{j}$ be the Clopper–Pearson
lower confidence bound for $q_{\bar t}(t_{j})$ at level $1-\eta/J$, i.e. the $\eta/J$ quantile of
$\mathrm{Beta}(X_{j},M_{j}-X_{j}+1)$ ($0$ if $X_{j}=0$), and $\Lambda_{j}=\Phi^{-1}(q^{\mathrm{lo}}_{j})$.

*The tail.* Fix a level $c^{\circ}\ge\bar t$ (in the computation $c^{\circ}=\max(\bar t,\hat c_{1},\dots,\hat c_{m})
=5.8321$, the largest critical value in use, so that one chain serves both rules). (T1) A chain of
intervals $[t_{a},t_{b}]$ covering $[t_{J},t_{g}]$, each with the CROWN lower bound of $S$ over the
Euclidean ball with centre $\tfrac{t_{a}+t_{b}}{2}v$ and radius $\tfrac{t_{b}-t_{a}}{2}+R$ strictly above
$c^{\circ}$. (T2) $t_{g}=(c^{\circ}+c+LR)/\kappa$, where $\kappa=S_{\infty}(v)>0$ is the consistency
certificate of the direction, $c\ge\sup|S-S_{\infty}|$ and $L\ge$ the Lipschitz constant of $S_{\infty}$
are the certified constants of `03_certify.py`: $c=100.6965$ (printed $100.70$; offset range
$[-100.697,82.826]$), $L=1735.636$ (printed $1735.64$). $c$ comes from interval propagation of the offsets
through every layer ($\max(a+x,0)-\max(a,0)\in[\min(l,0),\max(u,0)]$ for every $a$ and $x\in[l,u]$;
pooling and averaging preserve intervals), and $L$ from interval propagation of the gradient through every
activation pattern (rectifier derivatives in $\{0,1\}$, each max-pool routed to one of four inputs), so it
bounds $\|\nabla S\|$ and $\|\nabla S_\infty\|$ alike. Both were checked on 240,000 sampled maps at six
scales and a gradient attack ($c_{emp}=7.51$, $L_{emp}=3.82$). The ESM did not state these constants
before the 2026-09-24 insert; the insert now does. CROWN runs in float64 without outward rounding; the
smallest lower bound in any chain is 5.953 (saddle), a margin of 0.12 over $c^{\circ}$, far above
rounding error.

**Theorem (smoothing certificate, fixed threshold).** With probability at least $1-\eta$ over the net
draws, simultaneously for every $t\ge0$,
$$\Phi^{-1}\bigl(q_{\bar t}(t)\bigr)\ \ge\ \Lambda(t):=\max_{1\le j\le J}\bigl(\Lambda_{j}-|t-t_{j}|\bigr).$$
Consequently, on that event, $q_{\bar t}(t)\ge\beta$ for every $t$ in
$$\mathcal C=\bigcup_{j:\Lambda_{j}\ge z_{\beta}}\bigl[t_{j}-(\Lambda_{j}-z_{\beta}),\ t_{j}+(\Lambda_{j}-z_{\beta})\bigr]
\ \cup\ [t_{J},\infty),\qquad z_{\beta}=\Phi^{-1}(\beta),$$
provided (T1) and (T2) hold; the second set needs no Monte Carlo. Let $t_{\mathrm{cert}}$ be the left end of
the connected component of $\mathcal C$ that contains $+\infty$ and $n_{\mathrm{cert}}=\lceil4t_{\mathrm{cert}}^{2}\rceil$
(rounded **up**: $n\ge n_{\mathrm{cert}}$ implies $t=\sqrt n/2\ge t_{\mathrm{cert}}$; rounding to nearest
would claim sample sizes that are not certified). Then, with probability at least $1-\eta$, the
fixed-threshold test has power at least $\beta$ in the Gaussian map model at every $t\ge t_{\mathrm{cert}}$,
i.e. at every $n\ge n_{\mathrm{cert}}$. On the
additional event, of probability at least $.999$ over the paper's null draws, that $\bar t\ge c_{.95}$,
the same holds for the test that rejects when $S$ exceeds the true $.95$ null quantile $c_{.95}$, whose
size in the model is exactly $.05$ when $F_{0}$ is continuous at $c_{.95}$ and at most $.05$ otherwise.

*Proof.* (i) Each Clopper–Pearson bound fails with probability at most $\eta/J$ (the beta–binomial
identity: $P\{q^{\mathrm{lo}}_{j}>q\}\le\eta/J$ when $X_{j}\sim\mathrm{Bin}(M_{j},q)$). By the union bound
over the $J$ net points, all of $q^{\mathrm{lo}}_{j}\le q_{\bar t}(t_{j})$ hold with probability at least
$1-\eta$. Note that the net is deterministic, so $J$ is not data-dependent. (ii) On that event
$\Lambda_{j}\le\Phi^{-1}(q_{\bar t}(t_{j}))$ for every $j$, and by the Corollary to Lemma 1,
$\Phi^{-1}(q_{\bar t}(t))\ge\Phi^{-1}(q_{\bar t}(t_{j}))-|t-t_{j}|\ge\Lambda_{j}-|t-t_{j}|$ for every $j$
and $t$; take the maximum over $j$. (iii) $\Lambda_{j}-|t-t_{j}|\ge z_{\beta}$ iff
$|t-t_{j}|\le\Lambda_{j}-z_{\beta}$, which gives the intervals. (iv) Tail, (T1): for
$t\in[t_{a},t_{b}]$ and $\|e\|\le R$, $\|tv+e-\tfrac{t_{a}+t_{b}}{2}v\|\le|t-\tfrac{t_{a}+t_{b}}{2}|+\|e\|\le
\tfrac{t_{b}-t_{a}}{2}+R$, so $tv+e$ lies in the ball and $S(tv+e)>\bar t$ by the soundness of CROWN
(the same fact Proposition S.12 rests on); hence $q_{\bar t}(t)\ge P\{\|e\|\le R\}=1-\varepsilon\ge\beta$.
(T2): by Lemma 2 of the paper, $S(m)\ge S_{\infty}(m)-c$; $S_{\infty}$ is positively homogeneous and
$L$-Lipschitz, so $S_{\infty}(tv+e)=tS_{\infty}(v+e/t)\ge t\kappa-L\|e\|$. For $t>t_{g}$ and
$\|e\|\le R$, $S(tv+e)\ge t\kappa-LR-c>c^{\circ}\ge\bar t$, so again $q_{\bar t}(t)\ge1-\varepsilon$. Together (T1)
and (T2) cover $[t_{J},\infty)$. (v) The component of $\mathcal C$ containing $+\infty$ is an interval
$[t_{\mathrm{cert}},\infty)$ (a union of closed intervals and a closed half-line whose connected
component is an interval), and $t\mapsto n=4t^{2}$ is increasing. (vi) $\{S>\bar t\}\subseteq\{S>c_{.95}\}$
when $\bar t\ge c_{.95}$, so the power of the test with the true critical value is at least
$q_{\bar t}$. $\square$

**Theorem (smoothing certificate, bootstrap rule $B=199$, $\alpha=.05$).** Let $\hat c_{1}\ge\dots\ge
\hat c_{m}$ be the level bounds of Lemma 3 (joint probability $\ge1-\eta_{c}$), and at each net point let
$X_{jl}=\#\{i:S(t_{j}v+e_{ji})>\hat c_{l}\}$ (the same draws for all $l$), with Clopper–Pearson lower
bounds $q^{\mathrm{lo}}_{jl}$ at level $1-\eta/(Jm)$ and $\Lambda_{jl}=\Phi^{-1}(q^{\mathrm{lo}}_{jl})$.
Put
$$\underline\beta(t)=\sum_{l=1}^{m}\Delta g_{l}\,\Phi\Bigl(\max_{j}\bigl(\Lambda_{jl}-|t-t_{j}|\bigr)\Bigr).$$
Then with probability at least $1-\eta-\eta_{c}$, $\beta_{199,10}(tv)\ge\underline\beta(t)$ for every
$t\le t_{J}$, and $\beta_{199,10}(tv)\ge(1-\varepsilon)\,g(w_{m})$ for every $t\ge t_{J}$ when (T1)–(T2)
hold with $\bar t$ replaced by $\max_{l}\hat c_{l}$. The function $\underline\beta$ is Lipschitz in $t$ with
constant $\sup\varphi\cdot\sum_{l}\Delta g_{l}\le0.4$, so if $\underline\beta\ge\beta+0.2\,\Delta t$ on a grid of
spacing $\Delta t$ over $[t^{\circ},t_{J}]$, then $\beta_{199,10}(tv)\ge\beta$ for all $t\ge t^{\circ}$.

*Proof.* Lemma 4 with $P\{S(m)>\hat c_{l}\}=q_{\hat c_{l}}(t)$, the union bound over the $Jm$
Clopper–Pearson statements and the $m$ quantile bounds, and the Corollary to Lemma 1 applied to each
$q_{\hat c_{l}}$; the tail as in (iv) above, using $\hat c_{l}\le\max_{l}\hat c_{l}$. The Lipschitz claim:
each inner maximum is $1$-Lipschitz in $t$, $\Phi$ has derivative at most $\varphi(0)<0.4$, and the
$\Delta g_{l}$ are non-negative with sum at most $1$. $\square$

## 4. Reading the theorem

* **Where the tightness comes from.** Proposition S.12 bounds the *value* of $S$ by a linear relaxation
  through 3,200 rectifiers, and its slack is the compounded relaxation error (about 50 score units at the
  certified point). The smoothing certificate never bounds $S$ at all: it measures the power at net
  points and interpolates with a Lipschitz constant that is *exact* (Lemma 1 is sharp). Its only slacks
  are (a) the Clopper–Pearson gap at the net, of order $z_{\eta/J}\sqrt{q(1-q)/M}$ on the probability
  scale, and (b) the interpolation loss, at most half the net spacing on the $\Phi^{-1}$ scale. With
  spacing $0.1$ and $M=10^{5}$ near the certified point both are a few hundredths, so the certified $t$
  exceeds the true one by a few per cent and $n_{\mathrm{cert}}/n_{\mathrm{true}}$ is close to one.
* **What is spent.** $\eta$ is the probability that the certificate is wrong because of its own Monte
  Carlo; $\eta_{c}$ (or the $.001$ carried by $\bar t$) is the probability that a critical value bound
  is wrong; $\varepsilon$ is not a failure probability but a deterministic loss of power in the tail. All
  are fixed in advance and reported.
* **What it is not.** The model is Gaussian with identity covariance and a fixed direction; the
  finite-$n$ law of the deployed test is discrete, has covariance $\Sigma$, and uses a bootstrap
  critical value. Lemma 2 shows that the same argument runs verbatim for $N(\delta,\Sigma)$ along any ray
  in $\mathrm{range}(\Sigma)$, with the Mahalanobis length of the ray in place of $|t-t'|$; the map from
  the local drift to $n$ is then Theorem S.9's, not $n=4t^{2}$. Bridging to the finite-$n$ Bernoulli law
  would need a Berry–Esseen-type bound for the map, which is not attempted here.
* **Cost.** The Monte Carlo is one forward pass per draw (about $6$–$18$ million draws per direction,
  concentrated near the certified point), plus one CROWN pass per tail interval. No retraining and no
  bounding of the network near the signal strengths that matter.
