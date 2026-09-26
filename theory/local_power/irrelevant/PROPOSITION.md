# Proposition P (axis selection against a projection over all directions) — statement and proof

**Status 2026-09-25: restated for the projection statistic as Liu et al. (2024) define it** — the
angle-only weight of Escanciano (2006) on the raw covariates, which is the paper's primary projection
arm (`theory/addon_arms/proj_liu.R`, `proj_A_liu(std = FALSE)`; compiled copy `proj_liu_c.c`). The
previous version analysed the benchmark's first implementation, whose weight carried an extra length
factor $|x_i-x_l||x_j-x_l|$ on standardized covariates; it is kept for the archive in
`PROPOSITION_P_nnweight.md`; it is not in the paper. The rate for the angle-only weight was first derived as Remark (g) of the adversarial review
(`ADVERSARIAL_REVIEW.md`, checked by `07_printed_weight_check.py`); it is now the main statement, with
the proof below. Numbers: `08_proj_liu_raw.py` → `PROJ_LIU_RAW.json` (curve, rate, null quantile),
`PRINTED_WEIGHT_CHECK.log` (drift constant). Re-checked by a second adversarial review on 2026-09-25
(`ADVERSARIAL_REVIEW.md`, "Review 2 (printed weight)"). Notation follows Theorem L (`../THEOREM.md`; Theorem S.9 of
the supplement).

## 0. What is proved, in three lines

1. **Axis selection makes DeepGOF-1 blind to irrelevant covariates.** Once the selection rule picks the
   two relevant axes (A3), $\delta(h)$, $\Sigma$ and the limit power do not depend on the number $k$ of
   added independent covariates that carry no signal. The Hosmer–Lemeshow and Stukel limits, built on
   the linear predictor, are invariant for the same reason.
2. **The projection test is diluted.** Its statistic averages the residual process over all directions
   of the full $p$-dimensional covariate space. Its weight depends on angles only (homogeneous of degree
   zero in the Gram entries), its drift per unit $\tau^2$ falls like $Q^\circ/p^2$ (explicit constant),
   and its null standard deviation falls only like $1/p$ (explicit bounds). Its local signal-to-noise
   index therefore decays like $1/p$.
3. Hence the local sample size the projection test needs, relative to DeepGOF-1, grows at most linearly
   in $p$ (unconditionally), and grows linearly under a condition on the null quantile that is checked
   numerically but not proved.

Every statement is about the limit $n\to\infty$ taken at a fixed $k$; nothing is claimed when $k$ grows
with $n$.

## 1. Setting

**Design ("pure").** $x=(x_R,x_I)$, $x_R\in\mathbb R^2$ relevant, $x_I\in\mathbb R^k$ irrelevant, all
$p=2+k$ coordinates i.i.d. $N(0,1)$. The null model is logistic in $\tilde x=(1,x_R,x_I)$ with true value
$\beta_0=(\beta_{0R},0_k)$, the coefficients of both relevant covariates nonzero (so that the top pair of
A3 is unique): the irrelevant covariates are fitted but absent from the true predictor.
$\pi_0$ and $v=\pi_0(1-\pi_0)$ depend on $x_R$ only. (The niche study's $0.15$ coefficients: §6.)

**Local alternatives.** $\pi_n=\pi_0+n^{-1/2}h$, $h=h(x_R)$, $|h|\le Cv$ (Theorem L). In the niche study
$h=\tau vg$, $g=x_1x_2$ or $x_1^2-1$, $\tau=C\sqrt n$.

**Projected drift.** $\tilde h=h-v\,\tilde x^\top\mathcal I^{-1}b(h)$, $\mathcal I=E[v\tilde x\tilde x^\top]$,
$b(h)=E[\tilde xh]$. Then $E[\tilde h]=0$ and $E[\tilde h\tilde x]=0$ (unweighted). Put
$M=E[\tilde h\,x_Rx_R^\top]\in\mathbb R^{2\times2}$.

**The projection statistic** (Liu et al. 2024, Sec. 2.1, after Escanciano 2006). $T=n^{-2}\sum_{i,j,l}e_ie_jA_{ijl}$,
$e=y-\hat\pi$, computed on the raw covariates, with
$$A_{ijl}=\int_{S^{p-1}}\mathbf 1\{\omega^\top(x_i-x_l)\le0\}\,\mathbf 1\{\omega^\top(x_j-x_l)\le0\}\,d\omega$$
($d\omega$ the uniform probability measure). For nonzero $a=x_i-x_l$, $b=x_j-x_l$ this is
$A_{ijl}=(\pi-\angle(a,b))/(2\pi)$; it is $1/2$ when exactly one of $a,b$ is zero and $1$ when both are
(validated against a brute-force integral over $2\times10^5$ directions in `validate_proj_liu.log`). Hence
$$e^\top Ae=\sum_l E_\omega\Bigl[\Bigl(\sum_ie_i\,\mathbf 1\{\omega^\top(x_i-x_l)\le0\}\Bigr)^2\Bigr]:$$
the statistic averages a residual process over **all directions $\omega$ of $\mathbb R^p$** (and base
points $x_l$). A uniformly random direction puts only a fraction $|\omega_R|^2$, of mean $2/p$, of its
length on the plane of the two relevant covariates. The diagonal is $A_{ii}=1+(n-1)/2$ exactly, so the
diagonal part of $T$ is $(n+1)(2n^2)^{-1}\sum_ie_i^2$, a function of the residual sum of squares alone.

In terms of the Gram entries, $A_{ijl}=F^\circ(\alpha,\beta,\gamma)/(2\pi)$ with
$$F^\circ(\alpha,\beta,\gamma)=\pi-\arccos\frac{\gamma}{\sqrt{\alpha\beta}},\qquad
(\alpha,\beta,\gamma)=(|a|^2,|b|^2,a^\top b),$$
**homogeneous of degree zero**. A constant factor does not change the test, so all quantities below are
stated for the weight $F^\circ$; its population weight function is
$K(x,x')=E_zF^\circ(|x-z|^2,|x'-z|^2,(x-z)^\top(x'-z))$, $z$ an independent copy of $x$ ($A/n\to K/(2\pi)$).

**Limit of the statistic.** As for every residual quadratic form with a fixed bounded weight function
(Escanciano 2006 for this family), $T$ converges under $P_n$, after removal of the deterministic limit of
its diagonal part ($E v/2$), to $(2\pi)^{-1}\sum_j\lambda_j(Z_j+\tau c_j)^2-\text{const}$, $Z_j$ i.i.d.
$N(0,1)$, where $\lambda_j$ are the eigenvalues of $(I-\Pi)\mathcal A(I-\Pi)$ on $L^2(F_X)$,
$\mathcal A(x,x')=\sqrt{v(x)}K(x,x')\sqrt{v(x')}$, $\Pi$ the orthogonal projection onto
$\mathrm{span}\{\sqrt v\,\tilde x_j\}$, and $c_j$ the coordinates of $\tilde h/\sqrt v$. (This
representation is used, not re-derived.) Two scalars summarize it:
$$\sigma=\sum_j\lambda_jc_j^2=E[\tilde h(x)\tilde h(x')K(x,x')]\ \ (\text{drift per }\tau^2),\qquad
\nu^2=2\sum_j\lambda_j^2=2\|(I-\Pi)\mathcal A(I-\Pi)\|_{HS}^2\ \ (\text{null variance}),$$
and $e=\sigma/\nu$ is the **local signal-to-noise index** (free of the constant factor).

**Constants.** $\rho_0=\tfrac12$; $F^\circ$ at $(2,2,1)$, the bulk point divided by $p$:
$$F^\circ_{\alpha\beta}=\frac{1}{12\sqrt3}=0.048113,\qquad F^\circ_{\gamma\gamma}=\frac{1}{3\sqrt3}=0.192450,\qquad
F^\circ_{\alpha\gamma}=-\frac{1}{3\sqrt3},$$
$$Q^\circ=F^\circ_{\alpha\beta}(\operatorname{tr}M)^2+\tfrac12F^\circ_{\gamma\gamma}\|M\|_F^2,\qquad
C^\circ=4(F^\circ_{\alpha\beta})^2+\tfrac12(F^\circ_{\gamma\gamma})^2=\tfrac1{36}.$$
In the niche design with the added coefficients set to zero, $Q^\circ=0.00424$ (interaction), $0.00745$
(quadratic); $E v=0.18368$ (Gauss–Hermite $150\times150$, converged; Review 2 — the earlier $0.00744$ was
$0.2\%$ low; the bands below are unchanged at three decimals).

## 2. Statement

**Proposition P.** In the setting of §1, with every limit taken as $n\to\infty$ at fixed $k$, as
functions of the number $k$ of irrelevant covariates:

**(i) Axis selection gives invariance.** $\delta(h)$, $\Sigma$, $F_0$ and the limit power
$\bar\beta_{B,k}(\delta(h))$ of Theorem L are the same for every $k\ge0$ (A3 holds with the top pair
$(x_1,x_2)$ because the irrelevant population criteria $|\beta_{0j}|\sigma_j$ are zero). The limit powers
of the Hosmer–Lemeshow and Stukel tests are also the same for every $k$.

**(ii) The projection test is diluted at rate $1/p$.** $K$ is positive semidefinite, so every
$\lambda_j\ge0$. With $p=2+k$, as $k\to\infty$,
$$\sigma_k=\frac{Q^\circ}{p^2}\bigl(1+O(p^{-1/2})\bigr),\qquad
\frac{2(Ev)^2C^\circ}{p^2}\,(1+o(1))\ \le\ \nu_k^2\ \le\ \frac{C^\circ}{8p^2}\,(1+o(1)),$$
so that
$$\frac{Q^\circ}{\sqrt{C^\circ/8}}\,(1+o(1))\ \le\ p\,e_k\ \le\ \frac{Q^\circ}{Ev\,\sqrt{2C^\circ}}\,(1+o(1)).$$
Numerically: $p\nu_k\in[0.0433,0.0589]$; $p\,e_k\in[0.072,0.098]$ (interaction), $[0.126,0.172]$
(quadratic).

**Corollary P.1 (local sample size).** Consider the projection test with the exact null quantile
($B=\infty$), and let $q_k$ be the $.95$ quantile of its standardized null limit
$S_k=(\sum_j\lambda_jZ_j^2-\sum_j\lambda_j)/\nu_k$. Let $\tau_{50}^{P}(k)$ and $\tau_{50}^{DG}$ be the
amplitudes at which the projection test and DeepGOF-1 have limit power one half ($\tau_{50}^{DG}$ does
not depend on $k$ by (i)).
(a) Unconditionally, $\tau_{50}^{P}(k)^2\le\theta_1\nu_k/\sigma_k$ with
$\theta_1=\sqrt{19}+\sqrt2+\sqrt{2\sqrt2\sqrt{19}+3}=9.69$; hence $n^*_P(k)/n^*_{DG}=O(p)$.
(b) If $\underline q=\liminf_kq_k>1$ (not proved; checked numerically in §4),
$\tau_{50}^{P}(k)^2\ge\theta_0\nu_k/\sigma_k$ with $\theta_0=\underline q+\sqrt2-\sqrt{2\sqrt2\,\underline q+3}>0$,
so that
$$\frac{n^*_P(k)}{n^*_{DG}}=\Bigl(\frac{\tau_{50}^P(k)}{\tau_{50}^{DG}}\Bigr)^2\asymp p ,$$
and the ratio exceeds one for all large $k$.

**Remark (the mechanism).** The fit removes from $K$ everything that is constant or linear in either
argument. What is left of the departure is the **curvature** of $F^\circ$ in the Gram entries
$(|x-z|^2,|x'-z|^2,(x-z)^\top(x'-z))$, which are of size $p$ while the relevant covariates move them by
$O(1)$; $F^\circ$ is homogeneous of degree zero, so its second derivatives there are of order $1/p^2$, and
so is the drift. The noise is carried by the curvature along the $k$ irrelevant coordinates themselves:
$O(k^2)$ terms of size $O(1/p^2)$ each that survive the fit, a Hilbert–Schmidt norm of order $1/p$.
Drift over noise is of order $1/p$. DeepGOF-1 never forms these quantities; its cells are built on the two
selected axes.

## 3. Proof

### Step 0 — decoupling of the irrelevant block
Because $x_I\perp x_R$, $E[x_I]=0$ and $v,h$ depend on $x_R$ only,
$\mathcal I=\mathrm{diag}(\mathcal I_{RR},(Ev)I_k)$ and $b(h)=(b_R(h),0_k)$; hence
$\tilde h=h-v\,\tilde x_R^\top\mathcal I_{RR}^{-1}b_R$ is a function of $x_R$ alone, **the same for every
$k$**, and $E[\tilde h\,x_I]=E[\tilde h]E[x_I]=0$. Also $E[v\,x_I\,f(x_R)]=0$ for every $f\in L^2$.

### Step 1 — (i)
*DeepGOF-1.* The population cells $C_c$ are products of quantile intervals of $x_1,x_2$, so
$d_c=E[v1_{C_c}]$ and $a_c(h)=E[h1_{C_c}]$ do not involve $x_I$, and the $k$ added columns of $G$ are
$E[v\,x_I\,1_{C_c}]=0$ (Step 0). With the block forms,
$G\mathcal I^{-1}b(h)=G_R\mathcal I_{RR}^{-1}b_R(h)$ and $G\mathcal I^{-1}G^\top=G_R\mathcal I_{RR}^{-1}G_R^\top$,
the $k=0$ quantities; so $\delta(h)$, $\Sigma$, $F_0$ and $\bar\beta_{B,k}$ are those of $k=0$ (Theorem L,
item 4). A3: the population criteria are $|\beta_{0j}|\sigma_j$, zero for the irrelevant coordinates; under
$P_n$, $\hat\beta\to\beta_0$ by contiguity, so the selected pair is $(x_1,x_2)$ with probability tending to
one at each fixed $k$.
*Hosmer–Lemeshow.* The groups are the deciles of $\eta_0$, a function of $x_R$; the group masses, the
weights $w_g$ and the drifts $E[h1_g]$ do not involve $x_I$, and the added columns of the group–score
covariance are $E[v\,x_I1_g]=0$. The limit quadratic form (§S.rivals) is that of $k=0$.
*Stukel.* The added covariates $z(\eta_0)$ depend on $x_R$ only; $E[vz\,x_I^\top]=0$, so $\Omega$ and
$s(h)$ equal their $k=0$ values, and so does the non-centrality $s^\top\Omega^{-1}s$. $\square$

### Step 2 — $K$ is positive semidefinite
For nonzero $a,b$, $F^\circ(|a|^2,|b|^2,a^\top b)/(2\pi)=P(\omega^\top a\le0,\omega^\top b\le0)
=E_\omega[\phi_\omega(a)\phi_\omega(b)]$ with $\phi_\omega(a)=\mathbf 1\{\omega^\top a\le0\}$: an expectation
of products, hence a positive semidefinite function of $(a,b)$. $K(x,x')=E_zF^\circ(\dots)$ is an average of
such functions of $(x-z,x'-z)$ (the event $x=z$ has probability zero), hence positive semidefinite. Then
$(I-\Pi)\mathcal A(I-\Pi)\succeq0$ and $\lambda_j\ge0$. (The same representation shows that $A$ itself is
positive semidefinite at every $n$.)

### Step 3 — Taylor structure of $F^\circ$
$F^\circ$ is real-analytic on $\{\alpha,\beta>0,\ \gamma^2<\alpha\beta\}$ and homogeneous of degree zero:
$F^\circ(t\alpha,t\beta,t\gamma)=F^\circ(\alpha,\beta,\gamma)$. Its $m$-th derivatives are homogeneous of
degree $-m$, so on the box
$$\mathcal B_p=\{|\alpha/p-2|\le\tfrac14,\ |\beta/p-2|\le\tfrac14,\ |\gamma/p-1|\le\tfrac14\}$$
(convex, containing $G_0=(2p,2p,p)$, and inside the domain: there $|\gamma|/\sqrt{\alpha\beta}\le5/7$)
$|\partial^mF^\circ|\le C_mp^{-m}$. With $\rho=\gamma/\sqrt{\alpha\beta}$, $\psi(\rho)=\pi-\arccos\rho$,
$S=\sqrt{\alpha\beta}$: $F^\circ_\gamma=\psi'(\rho)/S$, $F^\circ_{\gamma\gamma}=\psi''(\rho)/S^2$,
$F^\circ_{\alpha\beta}=\{\rho\psi'(\rho)+\rho^2\psi''(\rho)\}/(4\alpha\beta)$,
$F^\circ_{\alpha\gamma}=-\{\psi'(\rho)+\rho\psi''(\rho)\}/(2\alpha S)$; at $(2,2,1)$, $\rho=\tfrac12$,
$\psi'=2/\sqrt3$, $\psi''=4/(3\sqrt3)$, which gives the constants of §1 (checked by finite differences in
`07_printed_weight_check.py`: $0.048113$ and $0.192450$).
Write $G=(\alpha,\beta,\gamma)=(|x-z|^2,|x'-z|^2,(x-z)^\top(x'-z))$, $\Delta=G-G_0$, and
$F^\circ(G)=\sum_{m=0}^{M}T_m(\Delta)+R_{M+1}$ with $T_m$ the Taylor terms at $G_0$. On $\{G\in\mathcal B_p\}$,
$|R_{M+1}|\le C p^{-(M+1)}\|\Delta\|^{M+1}$. The coordinates of $\Delta$ are centred Gaussian chaoses of
order at most two with variance $O(p)$, so $E\|\Delta\|^q=O(p^{q/2})$ and
$P(G\notin\mathcal B_p)\le Ce^{-cp}$ (Hanson–Wright). Off $\mathcal B_p$, $|F^\circ|\le\pi$ and every $T_m$ is
polynomial in $(p,\|\Delta\|)$, so by Cauchy–Schwarz the off-box part of any expectation below is
$O(e^{-cp/2}p^{C})$ and is omitted.

*The first two orders after averaging over $z$.* With $\Delta_\alpha=(|x|^2-p)-2x^\top z+(|z|^2-p)$,
$\Delta_\beta$ likewise, $\Delta_\gamma=x^\top x'-z^\top(x+x')+(|z|^2-p)$:
$E_z\Delta_\alpha=|x|^2-p$, $E_z\Delta_\gamma=x^\top x'$, and
$$E_z[\Delta_\alpha\Delta_\beta]=(|x|^2-p)(|x'|^2-p)+4x^\top x'+2p,\qquad
E_z[\Delta_\gamma^2]=(x^\top x')^2+|x+x'|^2+2p,$$
$$E_z[\Delta_\alpha\Delta_\gamma]=(|x|^2-p)\,x^\top x'+2x^\top(x+x')+2p,\qquad E_z[\Delta_\alpha^2]=\text{a function of }x .$$
Call **removable** ($\mathcal S$) any function of the form $f(x)+f'(x')+\sum_jg_j(x)\tilde x_j'+\sum_j\tilde x_jg_j'(x')$.
Everything above is removable except
$$P_{ab}=(|x|^2-p)(|x'|^2-p),\qquad P_{\gamma}=(x^\top x')^2-|x|^2-|x'|^2+p .$$
Since $F^{\circ\prime\prime}(G_0)=F^{\circ\prime\prime}(2,2,1)/p^2$,
$$E_z[T_0+T_1+T_2]=\mathcal S\text{-part}+J,\qquad J=\frac1{p^2}\Bigl(F^\circ_{\alpha\beta}P_{ab}+\tfrac12F^\circ_{\gamma\gamma}P_\gamma\Bigr).$$
In Hermite terms $P_{ab}=\sum_{a,b\le p}He_2(x_a)He_2(x_b')$ and
$P_\gamma=\sum_{a\ne b}x_ax_bx_a'x_b'+\sum_aHe_2(x_a)He_2(x_a')$, so
$\|P_{ab}\|^2=4p^2$, $\|P_\gamma\|^2=2p(p+1)$, $\langle P_{ab},P_\gamma\rangle=4p$ and
$$p^2\|J\|^2=C^\circ+\bigl(\tfrac12(F^\circ_{\gamma\gamma})^2+4F^\circ_{\alpha\beta}F^\circ_{\gamma\gamma}\bigr)/p .$$

### Step 4 — the drift $\sigma_k$
Removable terms have zero drift: $E[\tilde h\tilde h'f(x)]=E[\tilde hf]E[\tilde h']=0$ and
$E[\tilde h\tilde h'g(x)\tilde x_j']=E[\tilde hg]E[\tilde h\tilde x_j]=0$ (Step 0, §1). For $J$:
$E[\tilde h(|x|^2-p)]=E[\tilde h|x_R|^2]=\operatorname{tr}M$, and
$E[\tilde h\tilde h'(x^\top x')^2]=\sum_{a,b}E[\tilde hx_ax_b]^2=\|M\|_F^2$ (entries with an irrelevant index
vanish by Step 0). Hence $E[\tilde h\tilde h'J]=Q^\circ/p^2$ exactly.
Orders three and four: split each $\Delta$ coordinate as $\Delta^R+\Delta^I$, the part carried by
$(x_R,x_R',z_R)$ and the centred part carried by $(x_I,x_I',z_I)$, independent of the first. A product
survives $E[\tilde h\tilde h'\,\cdot\,]$ only if its relevant part is non-linear in both $x_R$ and $x_R'$,
which needs at least two relevant factors; a single irrelevant factor averages to zero. So
$E[\tilde h\tilde h'T_3]=O(p^{-3})$ (three relevant factors, $|F^{\circ\prime\prime\prime}|=O(p^{-3})$) and
$E[\tilde h\tilde h'T_4]=O(p^{-4}\cdot p)=O(p^{-3})$ (two irrelevant factors contribute $E[\Delta^I\Delta^I]=O(p)$).
The remainder: $|E[\tilde h\tilde h'R_5]|\le Cp^{-5}E[|\tilde h||\tilde h'|\|\Delta\|^5]=O(p^{-5/2})$, since
$|\tilde h|\le C'v(1+|x_R|)$. Therefore $\sigma_k=Q^\circ/p^2+O(p^{-5/2})$. $\square$

### Step 5 — the noise, upper bound
$(I-\Pi)\mathcal A(I-\Pi)$ annihilates every removable term: $\sqrt v f\otimes\sqrt v$ and
$\sqrt v g\otimes\sqrt v\tilde x_j$ have a factor in $\mathrm{range}(\Pi)$ ($\tilde x_0=1$). An orthogonal
projection does not increase the Hilbert–Schmidt norm and $v\le\frac14$, so with
$R=K-(\mathcal S\text{-part})-J=E_z[R_3]$,
$$\nu_k^2\le2\,E\bigl[v\,v'\,(J+R)^2\bigr]\le\tfrac18\bigl(\|J\|+\|R\|\bigr)^2,\qquad
\|R\|\le\|R_3\|_{L^2(x,x',z)}\le Cp^{-3}\bigl(E\|\Delta\|^6\bigr)^{1/2}=O(p^{-3/2})$$
(Jensen for $E_z$). With Step 3, $\nu_k^2\le\frac{C^\circ}{8p^2}(1+O(p^{-1/2}))$. $\square$

### Step 6 — the noise, lower bound
Let $J_I$ be $J$ with both sums restricted to irrelevant coordinates,
$J_I=p^{-2}\{F^\circ_{\alpha\beta}(|x_I|^2-k)(|x_I'|^2-k)+\frac12F^\circ_{\gamma\gamma}P_{\gamma,I}\}$, and
$\Phi=\sqrt{v(x)v(x')}\,J_I$. As a function of $x$, $J_I$ is a sum of $He_2(x_a)$ and $x_ax_b$ ($a\ne b$,
$a,b\in I$), and $E[v\tilde x_cHe_2(x_a)]=E[v\tilde x_cx_ax_b]=0$ for every $c$. So
$\Phi\perp\mathrm{range}(\Pi\otimes I)$ and $\mathrm{range}(I\otimes\Pi)$, and
$\langle(I-\Pi)\mathcal A(I-\Pi),\Phi\rangle=\langle\mathcal A,\Phi\rangle=E[vv'KJ_I]$. In
$K=\mathcal S\text{-part}+J+R$: removable terms are orthogonal to $vv'J_I$, and so are the parts of $J$
that involve a relevant coordinate; so $E[vv'KJ_I]=(Ev)^2\|J_I\|^2+E[vv'RJ_I]$ with
$|E[vv'RJ_I]|\le\frac1{16}\|R\|\|J_I\|=O(p^{-5/2})$, against a main term of order $p^{-2}$. Since
$\|\Phi\|^2=(Ev)^2\|J_I\|^2$, Cauchy–Schwarz gives
$$\nu_k^2\ge2\,\frac{\langle\mathcal A,\Phi\rangle^2}{\|\Phi\|^2}=2(Ev)^2\|J_I\|^2\bigl(1-O(p^{-1/2})\bigr),\qquad
p^2\|J_I\|^2=\frac{4k^2(F^\circ_{\alpha\beta})^2+\frac12(F^\circ_{\gamma\gamma})^2k(k+1)+4kF^\circ_{\alpha\beta}F^\circ_{\gamma\gamma}}{p^2}\to C^\circ. \ \square$$

### Step 7 — the index, and Corollary P.1
Dividing Step 4 by the square roots of Steps 5 and 6 gives (ii).
*Corollary.* In the limit, the standardized statistic is $U=S_k+2\tau X+\theta$ with
$\theta=\tau^2\sigma_k/\nu_k$, $X=\sum_j\lambda_jZ_jc_j/\nu_k$; $E[U]=\theta$ and
$$\operatorname{Var}U=1+\frac{4\tau^2\sum_j\lambda_j^2c_j^2}{\nu_k^2}\le1+\frac{4\tau^2\lambda_{\max}\sigma_k}{\nu_k^2}
=1+\frac{4\theta\lambda_{\max}}{\nu_k}\le1+2\sqrt2\,\theta,$$
using $\lambda_j\ge0$ (Step 2) and $\nu_k^2=2\sum\lambda_j^2\ge2\lambda_{\max}^2$. Power is increasing in $\tau$
(each $(Z_j+\tau c_j)^2$ is stochastically increasing in $\tau\ge0$ and the weights are nonnegative). By
Cantelli's inequality, $q_k\le\sqrt{19}$ always, and $P(U\le q_k)\le\frac12$ once
$\theta-q_k\ge\sqrt{1+2\sqrt2\theta}$, which holds for $\theta\ge\theta_1$: this is (a). Under (b),
$q_k\ge\underline q-\epsilon>1$ for all large $k$; $P(U>q_k)<\frac12$ once $q_k-\theta>\sqrt{1+2\sqrt2\theta}$,
which holds for $\theta<\theta_0(\underline q-\epsilon)$, and $\theta_0(\cdot)$ is continuous and increasing;
let $\epsilon\to0$. So $\theta_0\le\tau_{50}^P(k)^2\sigma_k/\nu_k\le\theta_1$ (up to $1+o(1)$), and (i), (ii)
give the ratio. $\square$

## 4. Numerical check

Pure design unless stated; all projection-test limits use the weight of §1 on raw covariates
(`08_proj_liu_raw.py`; the kernel reproduces `proj_A_liu(std = FALSE)` and `proj_liu_c.c` to $1.3\times10^{-9}$
on a 40-point design with tied rows).

* **Constants.** Finite differences of $F^\circ$ at $(2,2,1)$ reproduce $F^\circ_{\alpha\beta}=0.048113$,
  $F^\circ_{\gamma\gamma}=0.192450$ and $F^\circ_{\alpha\gamma}=-0.192450$. $Q^\circ=0.00424$ (interaction),
  $0.00745$ (quadratic); $C^\circ=1/36$; $Ev=0.1837$. Asymptotic band for $p\nu_k$: $[0.0433,0.0589]$; for
  $p\,e_k$: $[0.072,0.098]$ (interaction), $[0.126,0.172]$ (quadratic).
* **Drift** (`07_printed_weight_check.py` → `PRINTED_WEIGHT_CHECK.log`; Monte Carlo over the irrelevant
  block, 8 samples of 1,500 relevant draws, in-sample projection so that removable terms cancel exactly):
  $p^2\sigma_k=0.00426\pm0.00035$ / $0.00784\pm0.00065$ at $p=20$, $0.00417\pm0.00034$ / $0.00745\pm0.00053$
  at $p=50$, $0.00431\pm0.00028$ / $0.00758\pm0.00042$ at $p=200$ (interaction / quadratic), against
  $Q^\circ=0.00424$ / $0.00745$: **the rate $1/p^2$ and the constant are confirmed within about one SE.** (The
  $M=1000$ reference-sample discretisation cannot resolve the drift at $p=50$ — its values there are of
  either sign — so the drift is checked only by this Monte Carlo.)
* **Noise** ($\nu_k$ from the $1{,}000\times1{,}000$ matrix with its diagonal removed): $p\nu_k=0.048$, $0.051$
  ($p=10$), $0.047$, $0.047$ ($p=20$), $0.065$, $0.065$ ($p=50$); along the study-design curve ($p=2..10$,
  twelve samples each) $0.048$–$0.055$. Inside the band for $p\le20$; at $p=50$ about $10\%$ above the upper
  constant, within the $O(p^{-1/2})$ correction (and in line with `07`'s $0.065$ from a different sample).
  Order $1/p$, as proved; the constant is not located inside the band.
* **Diagonal and base points.** The weight is $1/2$ on the diagonal against about $1/3$ off it, so in an
  $M$-point discretisation the diagonal adds $O(1/M)$ to the drift and to $\sum\lambda_j^2$, whose true values
  are $O(p^{-2})$; it is removed (U-statistic discretisation). Kept, it lowers the curve's ratios by at most
  $0.18$ (interaction, $k=8$: $2.05$ against $1.86$) and raises $p\nu_k$ at $k=8$ from $0.048$ to $0.062$.
  Replacing the reference sample's own points as base points by $1{,}000$ or $4{,}000$ independent draws
  (same reference samples, seeds 1–2, $k=0$ and $8$) changes the ratios by at most $0.025$ and $\nu_k$ by at
  most $1.3\%$.
* **Index along the curve** (study design, twelve reference samples): $p\,e_k$ = $0.045, 0.053, 0.056, 0.068,
  0.071, 0.068$ (interaction) and $0.098, 0.108, 0.115, 0.137, 0.143, 0.143$ (quadratic) at $k=0,1,2,4,6,8$
  ($\pm0.002$–$0.008$): rising up to $k=4$ and flat after, at or just below the lower end of the asymptotic
  bands. The proposition claims nothing at such small $p$.
* **Corollary constants.** Standardized $.95$ null quantile $q_k$: $1.85$–$1.92$ ($p=2$), $1.74$–$1.80$
  ($p=10$, curve), $1.75$ and $1.78$ ($p=10$, pure), $1.71$ and $1.72$ ($p=20$), $1.62$ and $1.63$ ($p=50$): the
  condition $q_k>1$ holds at every $p$ checked, and with $\underline q=1.6$, $\theta_0=0.27$. The numerical
  $\theta^*=\tau_{50}^2\sigma_k/\nu_k$ (deployed $B=199$ rule) is $1.87$–$2.32$ along the curve, inside
  $[\theta_0,\theta_1]=[0.27,9.69]$ and slowly decreasing in $k$. The largest eigenvalue keeps $0.27$–$0.45$ of
  $\sum\lambda_j^2$ at $p\le10$ and $0.13$–$0.17$ at $p=50$, so the standardized null law is not normal; the
  corollary does not need it to be.
* **Invariances.** DeepGOF-1: $\max|\delta_8-\delta_0|=10^{-16}$, $\max|\Sigma_8-\Sigma_0|=7\cdot10^{-16}$
  (`RATE_CHECK.json`). Hosmer–Lemeshow and Stukel in the pure design: ratios $2.24$ and $1.04$
  (interaction), $2.69$ and $1.26$ (quadratic) at every $k\ge1$ (`ARE_IRRELEVANT.json`).
* **Discretisation size.** At $M=2000$ (two reference samples each) the curve's ratios are $0.59$, $0.73$
  ($k=0$) and $2.17$, $1.92$ ($k=8$) for the interaction and $0.54$, $0.68$ and $2.19$, $2.00$ for the
  quadratic, against the twelve-sample $M=1000$ values $0.74\pm0.03$, $0.64\pm0.02$, $2.05\pm0.11$,
  $1.95\pm0.10$: within the reference-sample spread (single-sample ranges $0.58$–$0.85$, $0.52$–$0.80$,
  $1.42$–$2.60$, $1.39$–$2.62$).

## 5. What is *not* claimed

* The rate is asymptotic in $k$. At $p\le10$ the Taylor variable $\|\Delta\|/p$ is not small, so the
  finite-$k$ curve of the paper is computed from the limit formulas at reference samples
  (`08_proj_liu_raw.py`), not from the rate.
* The representation of the projection statistic's limit as a weighted sum of noncentral $\chi^2_1$ with
  the population weight function $K$ is used, not re-derived here.
* The two noise bounds differ by the factor $(1/4)^2/(Ev)^2\approx1.85$; the true constant lies between
  them.
* The corollary's lower bound needs $\liminf q_k>1$; this is checked numerically, not proved. It is
  stated for the exact-quantile test; the deployed $B=199$ layer (Dufour's formula) is included in the
  numbers of `PROJ_LIU_RAW.json` but not in the bounds. The constants $\theta_0,\theta_1$ are crude
  (Cantelli).
* Correlated, non-Gaussian or categorical irrelevant covariates are not covered; with raw covariates of
  unequal scales the angle weight is not scale invariant, and nothing is claimed for that case. Nothing
  is claimed for the all-pairs reading of DeepGOF-1, which Theorem L does not cover.
* The limit drops the diagonal terms $e_i^2$ of the statistic (a deterministic shift as $n\to\infty$). At
  finite $n$ they are not negligible when $k$ is large: in the niche design the null standard deviation of
  the full statistic (fitted model, Bernoulli responses, four covariate samples, 400 replicates each) is
  $1.66$ ($n=200$) and $1.28$ ($n=500$) times its limit at $k=8$, and $0.99$, $1.00$ at $k=0$
  (`10_null_sd_finite_n.py` → `NULL_SD_FINITE_N.log`, Review 2). The excess variance ($1.76$ against $0.64$ times the limit variance) falls roughly like $1/n$. It is
  consistent with the limit over-predicting the simulated power at $C\ge0.5$ (`PRED_VS_SIM_LIU.json`); it
  is not modelled by the proposition.
* Invariance (i) is a large-$n$ statement at fixed $k$: it needs the selection to pick $(x_1,x_2)$, which in
  finite samples fails with a probability that grows with $k$ (the niche study: $97\%$ correct at $n=200$,
  $k=8$; $100\%$ at $n=500$).

## 6. The niche design ($c_{\rm irr}=0.15$)

With nonzero irrelevant coefficients $v$ depends weakly on $u=k^{-1/2}\sum_Ix_j$ ($\operatorname{Var}\eta_0$
rises from $1.64$ to $1.82$ at $k=8$) and Step 0 holds up to an $O(c_{\rm irr}^2k)$ perturbation. (i) then
holds approximately and (ii) is not re-proved; both are computed in the design itself (`ARE_IRRELEVANT.json`
for DeepGOF-1, HL and Stukel; `PROJ_LIU_RAW.json` for the projection test): DeepGOF-1's $\tau_{50}$ moves by
about one per cent over $k=0,\dots,8$, while Hosmer–Lemeshow and Stukel, whose statistics are built on the
fitted linear predictor that now carries the added covariates, drift upward (§3 of `RESULT.md`).

## 7. Credit

The decoupling is elementary. The expansion of a function of pairwise Gram entries around its
high-dimensional bulk value is the device of El Karoui (2010, *Ann. Statist.* 38, 1–50) for kernel random
matrices; here it is applied to the projection statistic's weight function, which is homogeneous of degree
zero. The direction representation is Escanciano's (2006) definition of the weight as an integral over the
sphere. Projection statistics were introduced to temper the curse of dimensionality of smoothing-based
specification tests (Escanciano 2006); a rate for their local signal-to-noise index under added
irrelevant covariates, with explicit constants, was not found in the sources checked and is derived here.
