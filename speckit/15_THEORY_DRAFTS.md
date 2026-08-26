<!-- 15_THEORY_DRAFTS.md — the two proposition drafts from workflow wf_c4eecd88-3e8
(2026-08-21). STATUS: DRAFTS, hostile-referee phase re-running after auth failures; Fable's own
hand-check of Prop 1 passed (see the verification note at the bottom). OPUS MUST TREAT THESE AS
DRAFTS until the attack/repair verdicts land or it re-verifies every step itself. -->

# Theory drafts — Propositions 1 and 2


## DRAFT PROPOSITION 1

### statement_latex

```
%======================================================================
% PROPOSITION 1 -- Validity of the deployed bootstrap calibration
% (statement block; proofs in the supplement)
%======================================================================

\paragraph{Setting.}
Let $\{P_\beta:\beta\in\Theta\}$, $\Theta\subseteq\mathbb{R}^p$, be the null model for a dataset
$D=(y,X)$ with the design $X$ held fixed; under $H_0$, $D\sim P_{\beta_0}$ for some unknown
$\beta_0\in\Theta$. Let $\hat\beta(\cdot)$ denote the maximum-likelihood estimator and let
$S(\cdot)$ be a \emph{frozen} measurable real-valued statistic: a fixed function of the dataset
alone, chosen before $D$ is observed (in our case $S$ contains its own refit). Given $D$, draw
$D^{*}_1,\dots,D^{*}_B$ i.i.d.\ from $P_{\hat\beta(D)}$ (same $X$) and set
\begin{equation}\label{eq:pval}
  p \;=\; \frac{1+\#\{b\le B: S(D^{*}_b)\ge S(D)\}}{B+1}.
\end{equation}
For $t\in\Theta$ let $G_t$ be the law of $S(D^{*})$ under $D^{*}\sim P_t$, write
$\bar G_t(s)=P_t\{S(D^{*})\ge s\}$, and define the \emph{calibration functional}
\begin{equation}\label{eq:pi}
  \pi(D)\;=\;\bar G_{\hat\beta(D)}\bigl(S(D)\bigr),
\end{equation}
the ideal ($B=\infty$) bootstrap $p$-value, a function of $D$ alone. If a required MLE fails to
exist we set $p=1$ (observed data) and $S(D^{*}_b)=+\infty$ (bootstrap data); both conventions can
only increase $p$, so all inequalities below are unaffected; the equality statements assume the
relevant MLEs exist almost surely.

\begin{proposition}[Finite-sample validity of the deployed calibration]\label{prop:validity}
Let $D\sim P_{\beta_0}$ and let $p$ be given by \eqref{eq:pval}.
\begin{enumerate}
\item[(i)] \textup{(Master condition.)} If $\pi(D)$ is \emph{super-uniform} under $P_{\beta_0}$,
i.e.
\[
  \mathbb{P}_{\beta_0}\{\pi(D)\le u\}\;\le\; u \qquad\text{for all } u\in(0,1),
\]
then for every $B\ge 1$ and every $\alpha\in(0,1)$,
\[
  \mathbb{P}_{\beta_0}\{p\le\alpha\}\;\le\;\frac{\lfloor \alpha(B+1)\rfloor}{B+1}\;\le\;\alpha .
\]
If $\pi(D)\sim\mathrm{Unif}(0,1)$, the first inequality is an equality.
\item[(ii)] \textup{(Sufficient conditions.)} Each of the following implies that $\pi(D)$ is
super-uniform:
\begin{enumerate}
\item[(a)] \emph{Pivotality:} $\mathcal{L}_{\beta_0}\bigl(S(D)\bigr)=G$ and
$G_{\hat\beta(D)}=G$ almost surely, for a single law $G$;
\item[(b)] \emph{Conditional pivotality given the fitted model:} for
$\mathcal{L}(\hat\beta)$-almost every $t$,
\[
  \mathcal{L}_{\beta_0}\bigl(S(D)\,\big|\,\hat\beta(D)=t\bigr)\;=\;G_t;
\]
\item[(c)] \emph{Conditional dominance (conservative direction):} for
$\mathcal{L}(\hat\beta)$-almost every $t$,
$\mathcal{L}_{\beta_0}\bigl(S(D)\mid\hat\beta=t\bigr)\preceq_{\mathrm{st}} G_t$.
\end{enumerate}
Under (a) the vector $\bigl(S(D),S(D^{*}_1),\dots,S(D^{*}_B)\bigr)$ is i.i.d., and under (b) it is
conditionally i.i.d.\ given $\hat\beta$; in either case it is exchangeable. Under (a) with $G$
atomless, or (b) with $G_t$ atomless a.s., $\pi(D)\sim\mathrm{Unif}(0,1)$ and the calibration is
exact.
\item[(iii)] \textup{(Unconditional quantitative bound; oracle AUC certificate.)} Without any
condition on $S$, for every $B\ge1$ and $\alpha\in(0,1)$,
\[
  \mathbb{P}_{\beta_0}\{p\le\alpha\}
  \;\le\;\frac{\lfloor \alpha(B+1)\rfloor}{B+1}\;+\;\kappa,
  \qquad
  \kappa:=\sup_{u\in(0,1)}\Bigl(\mathbb{P}_{\beta_0}\{\pi(D)\le u\}-u\Bigr)_{\!+},
\]
and
\[
  \kappa\;\le\; d_{\mathrm K}\bigl(\mathcal{L}(\pi),\mathcal{L}(\pi^{*})\bigr)
  \;\le\;\mathrm{TV}\bigl(\mathcal{L}(\pi),\mathcal{L}(\pi^{*})\bigr)
  \;\le\;2\Bigl(\mathrm{AUC}^{\star}-\tfrac12\Bigr),
\]
where $\pi^{*}:=\bar G_{\hat\beta(D)}\bigl(S(D^{*}_0)\bigr)$ for one additional independent
bootstrap draw $D^{*}_0\mid D\sim P_{\hat\beta(D)}$, $d_{\mathrm K}$ is the Kolmogorov distance,
and $\mathrm{AUC}^{\star}$ is the supremum, over all measurable classifiers, of the AUC for
separating draws of $\pi(D)$ from draws of $\pi^{*}$.
\end{enumerate}
\end{proposition}

\begin{corollary}[What $\mathrm{AUC}=\tfrac12$ does and does not certify]\label{cor:auc}
Let $\mathrm{AUC}^{\star}_S$ denote the optimal AUC of any classifier reading the score alone,
i.e.\ separating $\mathcal{L}_{\beta_0}\bigl(S(D)\bigr)$ from $\mathcal{L}\bigl(S(D^{*}_1)\bigr)$.
\begin{enumerate}
\item[(i)] \textup{(Necessity.)} $\mathrm{AUC}^{\star}_S=\tfrac12$ if and only if
$\mathcal{L}\bigl(S(D)\bigr)=\mathcal{L}\bigl(S(D^{*}_1)\bigr)$. In particular, exchangeability of
$\bigl(S(D),S(D^{*}_1),\dots,S(D^{*}_B)\bigr)$ --- e.g.\ under Proposition~\ref{prop:validity}(ii)(a)
or (ii)(b) --- forces $\mathrm{AUC}^{\star}_S=\tfrac12$; any classifier achieving
$\mathrm{AUC}>\tfrac12$ therefore refutes exchangeability, and quantitatively every classifier $f$
satisfies $\mathrm{AUC}(f)-\tfrac12\le \mathrm{TV}\bigl(\mathcal{L}(S(D)),\mathcal{L}(S(D^{*}))\bigr)$.
\item[(ii)] \textup{(Insufficiency.)} $\mathrm{AUC}^{\star}_S=\tfrac12$ does \emph{not} imply
validity: Example~\ref{ex:marginal} exhibits a null model, estimator and statistic with
$\mathrm{AUC}^{\star}_S=\tfrac12$ exactly, yet rejection probability approaching $\tfrac12$ at
every fixed level. The score-scale AUC is a necessary, not a sufficient, empirical check.
\item[(iii)] \textup{(Oracle sufficiency on the calibrated scale.)} By
Proposition~\ref{prop:validity}(iii), the size distortion at every level and every $B$ is at most
$2(\mathrm{AUC}^{\star}-\tfrac12)$ for the optimal classifier acting on the \emph{calibrated}
scores $\pi$ versus $\pi^{*}$; in particular $\mathrm{AUC}^{\star}=\tfrac12$ on the calibrated
scale certifies validity outright.
\end{enumerate}
\end{corollary}

\begin{example}[Marginal equality without validity]\label{ex:marginal}
Let $\Theta=\{-1,0,1\}$ and $D=(V,Z)\in\{-1,1\}\times\mathbb{R}$. Under $P_0$: $V\sim\mathrm{Unif}\{\pm1\}$
and $Z\mid V\sim N(cV,1)$; under $P_{\pm1}$: $V\sim\mathrm{Unif}\{\pm1\}$ independently of
$Z\sim N(\mp c,1)$. Take $\hat\theta(D)=V$ and $S(D)=Z$, with true value $\theta_0=0$. Then
$S(D)$ and $S(D^{*})$ have the identical marginal law $\tfrac12 N(c,1)+\tfrac12 N(-c,1)$, so
$\mathrm{AUC}^{\star}_S=\tfrac12$ exactly; yet for every $\alpha\ge 1/(B+1)$,
$\liminf_{c\to\infty}\mathbb{P}_{0}\{p\le\alpha\}\ge\tfrac12$. Conditionally on $\hat\theta$ the
observed score sits on the wrong side of its own bootstrap cloud half the time; the failure is
invisible marginally and visible only relative to the fitted model.
\end{example}

\begin{remark}[Relation to exact Monte Carlo tests]\label{rem:mc}
If $\beta_0$ were known and the replicates drawn from $P_{\beta_0}$, the vector
$\bigl(S(D),S(D^{*}_1),\dots\bigr)$ would be i.i.d.\ and \eqref{eq:pval} would be an exact Monte
Carlo test for \emph{any} $S$ (Barnard, 1963; Hope, 1968; Besag and Clifford, 1989). The deployed
procedure is instead a ``local'' Monte Carlo test at a point estimate in the sense of Dufour
(2006), and exchangeability genuinely fails in general: $D\sim P_{\beta_0}$ while
$D^{*}_b\mid D\sim P_{\hat\beta}$, so the replicates are i.i.d.\ only \emph{conditionally} on
$\hat\beta$, and equality of the marginal laws of $S(D)$ and $S(D^{*})$ is strictly weaker than
condition (ii)(b) (Example~\ref{ex:marginal}). Dufour restores exactness by maximizing the
$p$-value over the nuisance space; Proposition~\ref{prop:validity} instead certifies the
un-maximized procedure that is actually run.
\end{remark}

\begin{remark}[Exact similarity is impossible here; approximate validity is the right target]\label{rem:impossible}
In the logistic model with canonical link, $T=X^{\top}y$ is a complete sufficient statistic (full-rank
exponential family with natural parameter space $\mathbb{R}^p$), the likelihood depends on $y$ only
through $T$, and on the event that the MLE exists it is in bijection with $\hat\beta$. If
$\ker(X^{\top})\cap\{-1,0,1\}^n=\{0\}$ --- which holds for Lebesgue-almost-every design with at
least one continuously distributed covariate --- the fibers $\{y:X^{\top}y=t\}$ are singletons, so
conditioning on $\hat\beta$ determines $y$; this is precisely why co-sufficient sampling
(Engen and Lilleg{\aa}rd, 1997; Lindqvist and Taraldsen, 2005) is powerless in this model (Barber
and Janson, 2022). Two consequences. First, condition (ii)(b) can then hold only degenerately, and
no statistic $S$ is exactly pivotal. Second, by Neyman structure (Lehmann and Romano, 2005,
\S4.3), \emph{every} test with exact size $\alpha$ for all $\beta$ satisfies
$\mathbb{E}(\phi\mid T)=\alpha$ a.s., hence $\phi(y)=\alpha$ a.s.: the only exactly similar tests
ignore the data. Exact finite-sample similarity is therefore unattainable for any $S$, learned or
classical; the operative guarantee is part (iii), and an estimated size of $.052$ at nominal $.05$
is what the theory predicts as the achievable target, not a defect of the learned statistic.
\end{remark}

\begin{remark}[Direction of the distortion and the measured facts]\label{rem:direction}
Part (ii)(c) is one-sided: if the bootstrap replicates stochastically dominate the observed score
conditionally on $\hat\beta$, the test remains valid but conservative, paying in power. This is
the regime produced by training $S$ with true nulls as the null class (score-scale AUC $.62$,
deployed size $.006$): the mismatch surfaced as deflation, not inflation. Deployment-matched
training, whose null class reproduces the fit-resample-refit pipeline, gave score-scale AUC
$.5011$ and deployed size $.052$; the monotone relation between the AUC and both deployed size and
power across the fifteen trained rules is the empirical counterpart of the chain in part (iii).
\end{remark}

\begin{remark}[Which classifier to run]\label{rem:whichclassifier}
The three scales are not interchangeable. A classifier on raw scores tests the necessary marginal
condition of Corollary~\ref{cor:auc}(i). A classifier on the pair $(S,\hat\beta)$ is the natural
conditional test, but in the continuous-design logistic model its oracle AUC is degenerately $1$
(Remark~\ref{rem:impossible}: the observed pair lies on a measurable graph over $\hat\beta$),
so it certifies nothing at the oracle level. The scale that matches the deployed test is the
calibrated one: $\pi$ versus $\pi^{*}$, both estimable by an inner (double) bootstrap, and it is
on this scale that part (iii) converts an oracle AUC of $\tfrac12$ into a validity certificate.
An empirical classifier only \emph{lower}-bounds the oracle AUC, so a measured AUC near
$\tfrac12$ is an estimate that the slack $2(\mathrm{AUC}^{\star}-\tfrac12)$ is small, never a proof.
\end{remark}
```

### proof_latex

```
%======================================================================
% SUPPLEMENT: proofs for Proposition 1, Corollary 1, Example 1,
% and the claims of Remark S (impossibility)
%======================================================================

We use two elementary lemmas.

\begin{lemma}[Survival probability integral transform]\label{lem:pit}
Let $W\sim Q$ on $\mathbb{R}$ and $\bar Q(s)=Q([s,\infty))$. Then
$\mathbb{P}\{\bar Q(W)\le u\}\le u$ for all $u\in(0,1)$, with equality for all $u$ if $Q$ is
atomless.
\end{lemma}

\begin{proof}
Fix $u$ and let $A_u=\{s:\bar Q(s)\le u\}$. Since $\bar Q$ is nonincreasing, $A_u$ is an up-set:
$A_u=[s_u,\infty)$ or $(s_u,\infty)$ with $s_u=\inf A_u$ (the case $A_u=\emptyset$ is trivial).
If $s_u\in A_u$ then $\mathbb{P}\{W\in A_u\}\le \mathbb{P}\{W\ge s_u\}=\bar Q(s_u)\le u$. If
$s_u\notin A_u$ then $\mathbb{P}\{W\in A_u\}=\mathbb{P}\{W>s_u\}=\lim_{s\downarrow s_u}\bar Q(s)\le u$
because $\bar Q(s)\le u$ for every $s>s_u$. Since $\{\bar Q(W)\le u\}=\{W\in A_u\}$, the bound
follows. If $Q$ is atomless then $\bar Q(W)=1-F_Q(W)$ a.s.\ with $F_Q(W)\sim\mathrm{Unif}(0,1)$
by the probability integral transform, giving equality.
\end{proof}

\begin{lemma}[AUC and total variation]\label{lem:auc}
Let $P,Q$ be laws on a Polish space, and for a measurable $f$ let
$\mathrm{AUC}(f)=\mathbb{P}\{f(W)>f(W')\}+\tfrac12\mathbb{P}\{f(W)=f(W')\}$ with $W\sim P$,
$W'\sim Q$ independent. Then
\[
  \sup_f \mathrm{AUC}(f)\;\ge\;\tfrac12\bigl(1+\mathrm{TV}(P,Q)\bigr),
  \qquad
  \mathrm{AUC}(f)-\tfrac12\;\le\;\mathrm{TV}(P,Q)\ \ \text{for every } f .
\]
In particular $\sup_f\mathrm{AUC}(f)=\tfrac12$ if and only if $P=Q$.
\end{lemma}

\begin{proof}
\emph{Lower bound.} Choose $A$ with $P(A)-Q(A)=\mathrm{TV}(P,Q)$ and take $f=\mathbf{1}_A$:
\[
  \mathrm{AUC}(\mathbf{1}_A)=P(A)\bigl(1-Q(A)\bigr)
  +\tfrac12\bigl[P(A)Q(A)+(1-P(A))(1-Q(A))\bigr]
  =\tfrac12+\tfrac12\bigl(P(A)-Q(A)\bigr).
\]
\emph{Upper bound.} Let $g(w,w')=\mathbf{1}\{f(w)>f(w')\}-\mathbf{1}\{f(w')>f(w)\}$, so
$\mathrm{AUC}(f)-\tfrac12=\tfrac12\,\mathbb{E}\,g(W,W')$. If $W''\sim Q$ is independent of $W'$,
then $\mathbb{E}\,g(W'',W')=0$ by exchangeability of $(W'',W')$. For each fixed $w'$,
$\bigl|\int g(w,w')\,(P-Q)(dw)\bigr|\le \|g\|_\infty\,\|P-Q\|_{1}\le 2\,\mathrm{TV}(P,Q)$, hence
$\mathbb{E}\,g(W,W')=\mathbb{E}\,g(W,W')-\mathbb{E}\,g(W'',W')\le 2\,\mathrm{TV}(P,Q)$. The
final claim combines the two bounds.
\end{proof}

\begin{proof}[Proof of Proposition~\ref{prop:validity}]
\emph{Reduction to the calibration functional.} Conditionally on $D$, the events
$\{S(D^{*}_b)\ge S(D)\}$, $b=1,\dots,B$, are independent with common probability
$\bar G_{\hat\beta(D)}\bigl(S(D)\bigr)=\pi(D)$, because the $D^{*}_b$ are i.i.d.\
$P_{\hat\beta(D)}$ given $D$ and $S(D)$ is $D$-measurable. Hence
$N:=\#\{b:S(D^{*}_b)\ge S(D)\}$ satisfies $N\mid D\sim\mathrm{Bin}(B,\pi(D))$. Since $N$ is an
integer, $p\le\alpha \iff N+1\le\alpha(B+1)\iff N\le k-1$ with
$k:=\lfloor\alpha(B+1)\rfloor$, and therefore
\begin{equation}\label{eq:reduction}
  \mathbb{P}\{p\le\alpha\}=\mathbb{E}\bigl[h_k\bigl(\pi(D)\bigr)\bigr],
  \qquad
  h_k(x):=\mathbb{P}\{\mathrm{Bin}(B,x)\le k-1\}
  =\sum_{j=0}^{k-1}\binom{B}{j}x^{j}(1-x)^{B-j}.
\end{equation}
The polynomial $h_k$ maps $[0,1]$ to $[0,1]$ and is nonincreasing:
$h_k'(x)=-B\binom{B-1}{k-1}x^{k-1}(1-x)^{B-k}\le0$ (the standard binomial--beta identity; for
$k=0$, $h_0\equiv0$ and all claims are trivial). For $U\sim\mathrm{Unif}(0,1)$,
\[
  \mathbb{E}\,h_k(U)=\sum_{j=0}^{k-1}\binom{B}{j}\!\int_0^1 u^{j}(1-u)^{B-j}\,du
  =\sum_{j=0}^{k-1}\frac{1}{B+1}=\frac{k}{B+1}.
\]

\emph{(i).} Super-uniformity says $\pi(D)$ stochastically dominates $U$. Since $h_k$ is
nonincreasing, $\mathbb{E}\,h_k(\pi)\le\mathbb{E}\,h_k(U)=k/(B+1)\le\alpha$, by
\eqref{eq:reduction}. If $\pi\sim\mathrm{Unif}(0,1)$ the two expectations coincide.

\emph{(ii)(a).} On the a.s.\ event $G_{\hat\beta(D)}=G$ we have $\pi(D)=\bar G\bigl(S(D)\bigr)$
with $S(D)\sim G$; Lemma~\ref{lem:pit} with $Q=G$ gives super-uniformity, and uniformity when $G$
is atomless. For the distributional claim: $\mathcal{L}\bigl(S(D^{*}_{1:B})\mid D\bigr)=G^{\otimes B}$
is constant in $D$, so $S(D^{*}_1),\dots,S(D^{*}_B)$ are i.i.d.\ $G$ and independent of $D$; with
$S(D)\sim G$ the $(B{+}1)$-vector is i.i.d.

\emph{(ii)(b).} Conditionally on $\hat\beta=t$, the observed score has law $G_t$ by hypothesis and
$\pi(D)=\bar G_t\bigl(S(D)\bigr)$; Lemma~\ref{lem:pit}, applied under the conditional law, gives
$\mathbb{P}\{\pi\le u\mid\hat\beta=t\}\le u$ for a.e.\ $t$ (equality if $G_t$ atomless), and
integration over $\mathcal{L}(\hat\beta)$ finishes. For the distributional claim:
$\mathcal{L}\bigl(S(D^{*}_{1:B})\mid D\bigr)=G_{\hat\beta(D)}^{\otimes B}$ is
$\sigma(\hat\beta)$-measurable, hence $S(D^{*}_{1:B})\perp\!\!\!\perp D\mid\hat\beta$; combined
with (ii)(b), the $(B{+}1)$-vector is i.i.d.\ $G_t$ conditionally on $\hat\beta=t$, and a mixture
of i.i.d.\ vectors is exchangeable.

\emph{(ii)(c).} Fix $t$ with
$F_t:=\mathcal{L}(S(D)\mid\hat\beta=t)\preceq_{\mathrm{st}}G_t$ and fix $u$. As in
Lemma~\ref{lem:pit}, $\{\pi\le u\}=\{S(D)\in A^t_u\}$ with $A^t_u=\{s:\bar G_t(s)\le u\}$ an
up-set. Stochastic dominance gives $F_t(A^t_u)\le G_t(A^t_u)$, and the proof of
Lemma~\ref{lem:pit} gives $G_t(A^t_u)\le u$. Integrate over $t$.

\emph{(iii).} \emph{First inequality.} Writing
$h_k(x)=h_k(1)+\int_0^1(-h_k'(s))\,\mathbf{1}\{x\le s\}\,ds$ and taking expectations (Fubini),
\[
  \mathbb{E}\,h_k(\pi)-\mathbb{E}\,h_k(U)
  =\int_0^1\bigl(-h_k'(s)\bigr)\bigl[F_\pi(s)-s\bigr]\,ds
  \;\le\;\kappa\int_0^1\bigl(-h_k'(s)\bigr)ds
  =\kappa\,\bigl(h_k(0)-h_k(1)\bigr)\;\le\;\kappa,
\]
where $F_\pi(s)=\mathbb{P}\{\pi\le s\}$; combine with \eqref{eq:reduction}.
\emph{Second inequality.} Conditionally on $D$, $S(D^{*}_0)\sim G_{\hat\beta(D)}$ exactly, so
Lemma~\ref{lem:pit} applied conditionally gives
$\mathbb{P}\{\pi^{*}\le u\mid D\}\le u$, hence $F_{\pi^{*}}(u)\le u$ for all $u$: $\pi^{*}$ is
super-uniform \emph{by construction}. Therefore, pointwise in $u$,
\[
  F_\pi(u)-u=\bigl[F_\pi(u)-F_{\pi^{*}}(u)\bigr]+\bigl[F_{\pi^{*}}(u)-u\bigr]
  \le d_{\mathrm K}\bigl(\mathcal{L}(\pi),\mathcal{L}(\pi^{*})\bigr),
\]
and taking the supremum gives $\kappa\le d_{\mathrm K}$. \emph{Third inequality.}
$d_{\mathrm K}\le\mathrm{TV}$ always, and
$\mathrm{TV}\bigl(\mathcal{L}(\pi),\mathcal{L}(\pi^{*})\bigr)\le 2(\mathrm{AUC}^{\star}-\tfrac12)$
by the lower bound of Lemma~\ref{lem:auc} applied with $P=\mathcal{L}(\pi)$,
$Q=\mathcal{L}(\pi^{*})$.
\end{proof}

\begin{proof}[Proof of Corollary~\ref{cor:auc}]
(i) If the $(B{+}1)$-vector is exchangeable its coordinates share one marginal law, so
$\mathcal{L}(S(D))=\mathcal{L}(S(D^{*}_1))$; when the two laws coincide, $(f(S(D)),f(S(D^{*}_1)))$
is an exchangeable pair for every $f$, whence $\mathrm{AUC}(f)=\tfrac12$. The converse and the
quantitative bound are Lemma~\ref{lem:auc}. (ii) is Example~\ref{ex:marginal}, verified below.
(iii) is Proposition~\ref{prop:validity}(iii) restated.
\end{proof}

\begin{proof}[Details for Example~\ref{ex:marginal}]
Under $P_0$ write $Z=cV+\varepsilon$, $\varepsilon\sim N(0,1)\perp V$. The bootstrap law of the
score at $t=\pm1$ is $G_t=N(-ct,1)$, so the marginal law of $S(D^{*})$ is
$\tfrac12 N(-c,1)+\tfrac12 N(c,1)$, equal to the marginal law of $S(D)=Z$; by
Lemma~\ref{lem:auc}, $\mathrm{AUC}^{\star}_S=\tfrac12$ exactly. Next,
$\pi(D)=\mathbb{P}\{N(-cV,1)\ge Z\mid D\}=\bar\Phi(Z+cV)=\bar\Phi(\varepsilon+2cV)$. Conditionally
on $V=1$, $\pi=\bar\Phi(\varepsilon+2c)\to0$ almost surely as $c\to\infty$, so by
\eqref{eq:reduction} and dominated convergence, for any $\alpha\ge1/(B+1)$ (i.e.\ $k\ge1$),
\[
  \mathbb{P}_0\{p\le\alpha\}
  =\mathbb{E}\,h_k(\pi)
  \ \ge\ \tfrac12\,\mathbb{E}\bigl[h_k(\pi)\mid V=1\bigr]
  \ \xrightarrow[c\to\infty]{}\ \tfrac12\,h_k(0)=\tfrac12 .
\]
Finally, $\mathcal{L}(S(D)\mid\hat\theta=t)=N(ct,1)\ne N(-ct,1)=G_t$: condition (ii)(b) fails even
though the marginals agree, and a classifier reading the pair $(S,\hat\theta)$ separates the two
classes with AUC tending to $1$.
\end{proof}

\begin{proof}[Proof of the claims in Remark~\ref{rem:impossible}]
The null likelihood is $\exp\{\beta^{\top}X^{\top}y-\psi(\beta)\}$ with respect to counting
measure on $\{0,1\}^n$: a full-rank canonical exponential family in $T=X^{\top}y$ with natural
parameter space $\mathbb{R}^p$, so $T$ is complete sufficient (Lehmann and Romano, 2005,
Theorem~4.3.1). The likelihood, hence existence and uniqueness of the MLE, depends on $y$ only
through $T$; when the MLE exists it solves $X^{\top}\pi_{\mathrm{logit}}(X\beta)=T$, and
$\beta\mapsto X^{\top}\pi_{\mathrm{logit}}(X\beta)$ is injective for full-column-rank $X$ (strict
convexity of $\psi$), so $\hat\beta$ and $T$ generate the same $\sigma$-field on that event. If
$\ker(X^{\top})\cap\{-1,0,1\}^n=\{0\}$, then $X^{\top}y=X^{\top}y'$ forces
$y-y'\in\ker(X^{\top})\cap\{-1,0,1\}^n$, i.e.\ $y=y'$: fibers are singletons. For fixed nonzero
$v\in\{-1,0,1\}^n$, the set of designs with $X^{\top}v=0$ is Lebesgue-null (each continuous column
must be orthogonal to $v$; when an intercept is present, only $v$ with $\mathbf{1}^{\top}v=0$
survive that row, and the remaining constraints on any continuous column are still null), and a
union over the finitely many $v$ remains null. Now let $\phi(y)\in[0,1]$ be the rejection
probability of any (possibly randomized) test with $\mathbb{E}_\beta\phi=\alpha$ for all
$\beta\in\mathbb{R}^p$. By sufficiency $\mathbb{E}(\phi\mid T)$ does not depend on $\beta$, and by
completeness applied to $\mathbb{E}(\phi\mid T)-\alpha$ we get $\mathbb{E}(\phi\mid T)=\alpha$
a.s.\ (Neyman structure; Lehmann and Romano, 2005, Theorem~4.3.2). On singleton fibers
$\mathbb{E}(\phi\mid T)=\phi$, so $\phi(y)=\alpha$ a.s.: the test rejects with probability
$\alpha$ regardless of the data. Applied to the deployed test at any attainable level
$k/(B+1)$, this shows no choice of $S$ makes it exactly similar and nontrivial.
\end{proof}
```

### body_paragraph

```
Proposition 1 replaces the folklore justification of the deployed calibration with exact finite-sample statements, and it separates what is achievable from what is not. The classical route to exactness --- exchangeability of $(S(D),S(D^*_1),\dots,S(D^*_B))$, as in Monte Carlo tests with known null --- holds if $S$ is pivotal or, more sharply, conditionally pivotal given the fitted model; but the replicates are i.i.d.\ only conditionally on $\hat\beta$ while $D$ is drawn at the true $\beta_0$, so equality of the marginal score distributions is strictly weaker than what validity needs (Example 1 attains AUC exactly $\tfrac12$ on the score scale with rejection probability near $\tfrac12$). Worse, in the continuous-design logistic model no nontrivial exactly valid test exists at all: $X^\top y$ is complete sufficient with almost-surely singleton fibers, so by Neyman structure every exactly similar test ignores the data (Remark on impossibility; Barber and Janson, 2022) --- a deployed size of $.052$ at nominal $.05$ is the theoretically achievable target, not a defect. The operative guarantee is therefore quantitative: at every level and every $B$, the size distortion is at most the Kolmogorov distance between the calibrated score $\pi(D)$ and its by-construction super-uniform bootstrap counterpart $\pi^*$, hence at most $2(\mathrm{AUC}^\star-\tfrac12)$ for the best classifier separating the two on the $p$-value scale. This delineates exactly what the measured AUCs certify: the score-scale AUC of $.62$ refutes exchangeability outright, the $.5011$ figure is a necessary check whose empirical value estimates --- but, being a lower bound on the oracle, cannot prove --- that the validity slack is small, and the monotone AUC--size relation across the fifteen trained rules is the empirical counterpart of this chain.
```

### conditions_honest

```
ASSUMED (stated in the setting, none hidden): (1) S is frozen: a fixed measurable function of the dataset chosen before D is observed, trained only on independent simulations -- if S touched D, every statement fails; (2) fixed design X, bootstrap at the same X (random X handled by conditioning throughout); (3) parametric bootstrap with refit inside S, i.i.d. replicates given D; (4) MLE non-existence handled by explicit conventions (observed failure -> p=1; bootstrap failure -> S*=+infinity) that can only inflate p, so all inequalities survive; equality claims assume a.s. existence; (5) validity is marginal over both data and bootstrap randomness; (6) all laws here are discrete (y in {0,1}^n), so exact uniformity of pi never literally holds -- every result is stated in the conservative (super-uniform) direction, which needs no continuity, and atoms only push p upward via the >= convention.

WHAT IS PROVED UNCONDITIONALLY: part (iii) -- size distortion at every alpha and every B is at most kappa = sup_u (P(pi<=u)-u)_+, which is at most the Kolmogorov distance between L(pi) and L(pi*), which is at most 2(AUC*-1/2) for the OPTIMAL classifier on the calibrated (p-value) scale. pi* is super-uniform by construction, with no assumptions.

WHAT AUC = 1/2 CERTIFIES AND WHAT IT DOES NOT: (a) Score-scale AUC (the quantity actually measured: .62 and .5011) tests only marginal equality of L(S(D)) and L(S(D*)). Marginal equality is NECESSARY for exchangeability (so the .62 measurement rigorously refutes it) but NOT SUFFICIENT for validity -- Example 1 has optimal score-scale AUC exactly 1/2 and size tending to 1/2. (b) The pair-scale classifier (S, beta-hat) is the natural conditional diagnostic but is degenerate in the flagship model: its oracle AUC is 1 because beta-hat a.s. determines y, so it certifies nothing. (c) The calibrated scale (pi vs pi*, i.e. classify p-values, estimable by an inner/double bootstrap) is the scale on which the certificate is real: oracle AUC = 1/2 there implies exact validity, and oracle AUC = 1/2 + delta implies size distortion <= 2*delta. (d) In ALL cases an empirical classifier only LOWER-bounds the oracle AUC; a measured AUC near 1/2 is evidence, never proof. The paper must say "estimated certificate", and the monotone AUC-size relation across 15 rules is the empirical support for using it as a surrogate.

THE HONEST WEAKENING FORCED BY THE MODEL: exact conditional pivotality (condition ii-b) is degenerate in logistic regression with continuous covariates -- conditioning on beta-hat (equivalently on the complete sufficient statistic X'y, singleton fibers for a.e. design) determines the data, so ii-b cannot hold nondegenerately, no S is exactly pivotal, and by Neyman structure NO nontrivial test is exactly similar. Consequently: exact size .05 is unattainable by anything; .052 is the right kind of number; conditions ii-a/ii-b are idealizations that can hold only approximately, and the paper's operative theorem is the quantitative part (iii). Part ii-c (conditional stochastic dominance) is the one-sided condition that licenses the .62-AUC rule's observed behavior: valid but conservative (size .006), paying in power.

ALSO NOT CLAIMED: no claim that AUC and size are monotonically related in general (that is reported as a measured fact about these 15 rules, consistent with but not implied by part iii); no claim that pi-superuniformity is necessary for finite-B validity (it is necessary only in the large-B limit, and we do not push that direction in the statement); no asymptotic claims at all -- everything is finite-sample.
```

### literature

```
CITED FOR (verified where noted):

- Barnard, G.A. (1963). Discussion of Bartlett, JRSS-B 25, 294. CITED FOR: origin of the exact Monte Carlo test (rank of observed statistic among simulated replicates).
- Hope, A.C.A. (1968). A simplified Monte Carlo significance test procedure. JRSS-B 30, 582-598. CITED FOR: formalization and power of MC significance tests.
- Dwass, M. (1957). Modified randomization tests for nonparametric hypotheses. Ann. Math. Statist. 28, 181-187. CITED FOR: the finite-resampling rank p-value with the (1+N)/(B+1) form (optional companion to Barnard/Hope).
- Besag, J. and Clifford, P. (1989). Generalized Monte Carlo significance tests. Biometrika 76, 633-642. CITED FOR: exchangeability (not i.i.d.-ness) as the operative condition for exact rank p-values.
- Dufour, J.-M. (2006). Monte Carlo tests with nuisance parameters: a general approach to finite-sample inference and nonstandard asymptotics. J. Econometrics 133(2), 443-477. [VERIFIED: https://ideas.repec.org/a/eee/econom/v133y2006i2p443-477.html] CITED FOR: the deployed procedure = his "local MC test at a point estimate", NOT exact in general; exactness via maximized MC instead. PARTIAL-SCOOP FLAG: Dufour already owns "bootstrap-at-estimate is not finite-sample exact" and its asymptotic validity; our novelty must be pitched as the checkable conditions + the quantitative AUC certificate for the UN-maximized procedure, not the invalidity observation itself.
- Beran, R. (1988). Prepivoting test statistics: a bootstrap view of asymptotic refinements. JASA 83, 687-697. CITED FOR: pi* is precisely a prepivoted root; the classifier on the calibrated scale operationalizes "how uniform is the prepivoted p-value"; deployment-matched training is in effect learning a prepivoted statistic. FLAG: overlap in spirit (uniformity of bootstrap p-values); his results are asymptotic refinements, ours finite-sample bounds -- say so explicitly.
- Hall, P. and Wilson, S.R. (1991). Two guidelines for bootstrap hypothesis testing. Biometrics 47, 757-762. CITED FOR: resample under the null + pivot for accuracy; the classical heuristics our Prop 1 makes exact.
- Barber, R.F. and Janson, L. (2022). Testing goodness-of-fit and conditional independence with approximate co-sufficient sampling. Ann. Statist. 50(5), 2514-2544. [VERIFIED: https://projecteuclid.org/journals/annals-of-statistics/volume-50/issue-5/Testing-goodness-of-fit-and-conditional-independence-with-approximate-co/10.1214/22-AOS2187.full ; arXiv:2007.09851] CITED FOR: (1) the singleton-fiber observation -- CSS is powerless for logistic regression with continuous covariates -- which powers our impossibility remark; (2) the nearest methodological neighbor. SCOOP FLAG (the serious one): their Theorem-1-style guarantee also has the "alpha + distance" form (Type I inflation bounded by an expected TV between an ideal and an implemented conditional law). Differentiate precisely: they MODIFY THE SAMPLER (aCSS) to make the distance small; we KEEP the plain parametric bootstrap and CERTIFY it for a frozen learned S, via the p-value-scale Kolmogorov reduction (their joint-TV-type bound is vacuous here, which is exactly why our reduction to pi is needed) and the classifier-AUC surrogate. Must be cited and compared in the body, not buried. Also check/cite the follow-up line: arXiv:2309.08063 (aCSS with regularization) and arXiv:2511.05281 (conditioning on posterior samples for frequentist GOF) -- active neighboring literature; neither states the AUC certificate, but both must appear in related work.
- Engen, S. and Lillegard, M. (1997). Stochastic simulations conditioned on sufficient statistics. Biometrika 84, 235-240. [VERIFIED] CITED FOR: co-sufficient sampling, the classical exact route our Remark says is unavailable here.
- Lindqvist, B.H. and Taraldsen, G. (2005). Monte Carlo conditioning on a sufficient statistic. Biometrika 92(2), 451-464. [VERIFIED] CITED FOR: same, general theory of conditional Monte Carlo given sufficient statistics.
- Lehmann, E.L. and Romano, J.P. (2005). Testing Statistical Hypotheses, 3rd ed., Springer, Sections 4.3 (Theorems 4.3.1-4.3.2). CITED FOR: exponential-family completeness and Neyman structure, used verbatim in the impossibility remark.
- Lopez-Paz, D. and Oquab, M. (2017). Revisiting classifier two-sample tests. ICLR. CITED FOR: classifier-based two-sample discrepancies; grounds reading a trained classifier's AUC as an estimated distributional distance. (Friedman 2004, Stanford tech report, as the precursor -- optional.)
- Hermans, J., Delaunoy, A., Rozet, F., Wehenkel, A., Begy, V., Louppe, G. (2022). A Trust Crisis in Simulation-Based Inference? Your Posterior Approximations Can Be Unfaithful. TMLR. [VERIFIED: https://arxiv.org/abs/2110.06581] CITED FOR: related-work paragraph only -- classifier/coverage diagnostics for learned inferential procedures in SBI; analogous diagnostic culture, different target (posterior coverage, not test size). Not a scoop.

NOT CITED, CONSIDERED AND REJECTED: conformal-prediction exchangeability lemmas (Vovk et al.) -- same rank lemma, but Besag-Clifford is the right lineage for MC tests and keeps the bibliography on-venue; King (1988)/Andrews-Ploberger (1994) belong to Proposition 2, not here.

SOOP-CHECK SUMMARY: nothing found that states the two load-bearing novelties of Prop 1 -- (i) the reduction of deployed finite-B size distortion to the Kolmogorov distance between calibrated scores kappa <= d_K(L(pi), L(pi*)) with pi* super-uniform by construction, and (ii) the two-sided AUC delineation (score-scale AUC = necessary only, with an explicit counterexample; calibrated-scale oracle AUC = sufficient, with the 2(AUC-1/2) constant). Barber-Janson 2022 is the closest and must be engaged directly; Dufour 2006 owns the invalidity-in-general observation and the maximized alternative; Beran 1988 owns the prepivoting idea asymptotically.
```


## DRAFT PROPOSITION 2

### statement_latex

```
% ---------- Setting (shared notation, to precede Proposition 2) ----------
% M : data -> T \subset R^{KxK} is the fixed residual-map transform (Section 2).
% The training generator defines two laws on the map space:
%   P_0  = law of M(D) when D is drawn from the null class of the training library
%          (null model, design/parameter configurations drawn from the library measure \Pi_0);
%   Q_g  = law of M(D) under departure g, nuisance configurations drawn as in training;
%   \bar{Q} = \int Q_g \, d\Pi(g), where \Pi is the probability measure on departure
%          functions induced by the training library.

\begin{proposition}[$\Pi$-average-power optimality of the idealized rule]
\label{prop:avgpower}
Let $P_0$ and $\bar{Q}=\int Q_g\,d\Pi(g)$ be the null-class and alternative-class laws of
the residual map $M(D)\in\mathcal{T}$ under the training generator, and suppose
$P_0$ and $\bar{Q}$ admit densities $p_0$ and $\bar q$ with respect to a common
$\sigma$-finite dominating measure $\mu$ on $\mathcal{T}$. Define the balanced-class
Bayes posterior
\[
\eta^\star(t) \;=\; \frac{\bar q(t)}{p_0(t)+\bar q(t)},
\qquad t\in\mathcal{T},\ p_0(t)+\bar q(t)>0,
\]
and the mixture likelihood ratio $\Lambda(t)=\bar q(t)/p_0(t)$, with the convention
$\Lambda(t)=\infty$ on $\{p_0=0<\bar q\}$.
\begin{enumerate}
\item[(i)] $\eta^\star$ is the $\mu$-a.e.\ unique minimizer of the balanced-class
population cross-entropy risk over all measurable $f:\mathcal{T}\to[0,1]$, and
$\eta^\star=\Lambda/(1+\Lambda)$ is a strictly increasing transform of $\Lambda$
on $[0,\infty]$.
\item[(ii)] Fix $\alpha\in(0,1)$. Any test of the form
$\varphi^\star(t)=\mathbf{1}\{\eta^\star(t)>c\}+\gamma\,\mathbf{1}\{\eta^\star(t)=c\}$,
with $c\in[0,1]$ and $\gamma\in[0,1]$ chosen so that
$\mathbb{E}_{P_0}\varphi^\star(M)=\alpha$, maximizes the $\Pi$-average power
\[
\int \mathbb{E}_{Q_g}\!\bigl[\varphi(M)\bigr]\,d\Pi(g)
\]
over all (possibly randomized) tests $\varphi:\mathcal{T}\to[0,1]$ satisfying
$\mathbb{E}_{P_0}\varphi(M)\le\alpha$. If moreover $\Lambda(M)$ has a continuous
distribution under $P_0$, the non-randomized test
$\mathbf{1}\{\eta^\star(M)>c_\alpha\}$, with $c_\alpha$ the $(1-\alpha)$-quantile of
$\eta^\star(M)$ under $P_0$, attains the maximum exactly.
\end{enumerate}
\end{proposition}
```

### proof_latex

```
\begin{proof}
(i) With balanced class prior, the population cross-entropy risk of a measurable
$f:\mathcal{T}\to[0,1]$ is
\[
R(f)\;=\;-\tfrac12\,\mathbb{E}_{P_0}\log\bigl(1-f(M)\bigr)
        -\tfrac12\,\mathbb{E}_{\bar Q}\log f(M)
      \;=\;\int \Bigl[-\tfrac12\,p_0(t)\log\bigl(1-f(t)\bigr)
        -\tfrac12\,\bar q(t)\log f(t)\Bigr]\,\mu(dt).
\]
For each $t$ with $p_0(t)+\bar q(t)>0$ the integrand is strictly convex in
$f(t)\in(0,1)$ and is minimized uniquely at
$f(t)=\bar q(t)/\{p_0(t)+\bar q(t)\}=\eta^\star(t)$; hence $\eta^\star$ minimizes
$R$ and any minimizer agrees with $\eta^\star$ $\mu$-a.e.\ on
$\{p_0+\bar q>0\}$ (this is the standard population property of the logarithmic
scoring rule). Since $x\mapsto x/(1+x)$ is strictly increasing on $[0,\infty]$
(with value $1$ at $\infty$), $\eta^\star=\Lambda/(1+\Lambda)$ is a strictly
increasing transform of $\Lambda$, proving (i).

(ii) For any test $\varphi:\mathcal{T}\to[0,1]$, Tonelli's theorem (the integrand
lies in $[0,1]$) gives
\[
\int \mathbb{E}_{Q_g}\!\bigl[\varphi(M)\bigr]\,d\Pi(g)
=\int\!\!\int \varphi\,dQ_g\,d\Pi(g)
=\mathbb{E}_{\bar Q}\bigl[\varphi(M)\bigr].
\]
Thus maximizing $\Pi$-average power subject to
$\mathbb{E}_{P_0}\varphi(M)\le\alpha$ is exactly the problem of testing the simple
hypothesis $P_0$ against the simple alternative $\bar Q$ at level $\alpha$. By the
Neyman--Pearson lemma \citep[Theorem~3.2.1]{LehmannRomano}, a test is most
powerful at level $\alpha$ if (and, up to $(P_0+\bar Q)$-null modification off the
boundary set, only if) it has the form
$\mathbf{1}\{\Lambda(t)>k\}+\gamma'\,\mathbf{1}\{\Lambda(t)=k\}$ with size exactly
$\alpha$. By the strict monotonicity established in (i),
$\{\eta^\star>c\}=\{\Lambda>k\}$ and $\{\eta^\star=c\}=\{\Lambda=k\}$ for
$c=k/(1+k)$, so $\varphi^\star$ is of Neyman--Pearson form with size $\alpha$ and
therefore maximizes $\mathbb{E}_{\bar Q}\varphi(M)$, which equals the
$\Pi$-average power. Finally, if $\Lambda(M)$ is continuously distributed under
$P_0$, then so is $\eta^\star(M)$, the boundary set $\{\eta^\star(M)=c_\alpha\}$
has $P_0$-probability zero, and the non-randomized test
$\mathbf{1}\{\eta^\star(M)>c_\alpha\}$ has size $\alpha$ and is of
Neyman--Pearson form, hence attains the maximum.
\end{proof}

% ---------- Remarks (print immediately after the proposition in the body/supplement) ----------

\begin{remark}[Scope of the optimality claim]\label{rem:scope}
Proposition~\ref{prop:avgpower} is an optimality statement \emph{within the class
of tests measurable with respect to the residual map} $M(D)$; the map is a lossy
compression of $D$, and no claim is made against tests that use the full data.
Moreover, $P_0$ is the training-library mixture null, so the size constraint is an
\emph{average} size over the library's null configurations, in the weighted-power
sense of \citet{Wald1943}; finite-sample size control at the user's actual null
configuration is the business of the bootstrap calibration and
Proposition~\ref{prop:validity}, not of this proposition. The two results are
complementary because the deployed test uses only the \emph{ordering} induced by
the score: replacing the $P_0$-quantile $c_\alpha$ by the data-dependent bootstrap
critical value changes the size guarantee but not the ordering, so
Proposition~\ref{prop:avgpower} identifies what the trained score aims at, while
Proposition~\ref{prop:validity} identifies what the deployed calibration
guarantees.
\end{remark}

\begin{remark}[Idealization gap I: the network estimates the Bayes score]
\label{rem:estimation}
The deployed statistic $S$ is a convolutional network fitted by empirical
cross-entropy minimization on a finite simulated corpus, followed by distillation
into a single frozen network \citep{Hinton2015}; it therefore \emph{estimates}
$\eta^\star$ and inherits three discrepancies: approximation error of the
architecture class, estimation error from the finite corpus, and optimization
error. Because the test is invariant to strictly increasing transforms of the
score, only the \emph{ordering} induced by $S$ matters, and the $\Pi$-average
power of the deployed ordering degrades continuously with its rank disagreement
from $\Lambda$; quantitative transfers from classification risk to testing
performance are available in the two-sample setting
\citep{Kim2021,LopezPaz2017}, and the use of discriminative classifiers as
likelihood-ratio estimators is standard in simulation-based inference
\citep{Cranmer2015}. We make no claim that the finite-corpus network attains the
oracle power of Proposition~\ref{prop:avgpower}.
\end{remark}

\begin{remark}[Idealization gap II: the prior is the training library]
\label{rem:prior}
The optimality criterion is average power with respect to the library-induced
measure $\Pi$, exactly the weighted-power framework of \citet{Wald1943},
\citet{King1987}, and \citet{AndrewsPloberger1994}: choosing the training library
is choosing where power is spent. Against departures outside (or heavily
down-weighted by) the support of $\Pi$, Proposition~\ref{prop:avgpower} confers
no guarantee, and a different library yields a different optimal ordering.
Finally, under the deployment-matched training of Section~\ref{sec:matched} the
null class is the law of the map on bootstrap replicates rather than on true null
draws; Proposition~\ref{prop:avgpower} applies verbatim with $P_0$ redefined as
that law, so matched training changes the reference null the score discriminates
against --- which is precisely what aligns the score with the deployed bootstrap
calibration --- without affecting the form of the optimality statement.
\end{remark}
```

### body_paragraph

```
Proposition~\ref{prop:avgpower} makes precise what the trained statistic is aiming at. The training library plays the role of a prior $\Pi$ over departure functions, and the balanced-class Bayes score for separating the library's null class from its $\Pi$-mixture alternative is a strictly increasing function of the mixture likelihood ratio of the residual map; thresholding it is therefore the Neyman--Pearson test of the mixture null against the mixture alternative, and it maximizes $\Pi$-average power among all tests that are measurable with respect to the map. This is the classical weighted-average-power notion of optimality \citep{Wald1943,King1987,AndrewsPloberger1994}, imported unchanged: training a classifier on simulated departures is not a heuristic but an implicit choice of the weight measure over alternatives. The claim is deliberately modest in two respects, recorded in Remarks~\ref{rem:scope}--\ref{rem:prior}: the optimality is relative to the lossy residual map, with average (not configuration-wise) size under the training null --- finite-sample size at the user's null is delivered by the bootstrap calibration of Proposition~\ref{prop:validity} --- and the deployed network only estimates the Bayes score from a finite corpus, so the oracle power is an idealization, not a finite-sample guarantee.
```

### conditions_honest

```
ASSUMED: (1) P0 and Q-bar (the null-class and Pi-mixture alternative laws of the residual map under the TRAINING generator) have densities w.r.t. a common sigma-finite dominating measure -- harmless here since the map lives on R^{KxK}, but it must be stated. (2) The 'null' in the proposition is the training-library mixture null P0, so the size constraint is AVERAGE size over the library's null configurations (Wald's weighted framework), not size at the user's actual (beta, X). (3) Exact attainment of level alpha needs either randomization on the NP boundary or a continuously distributed likelihood ratio under P0; both variants are in the statement. (4) Optimality holds ONLY among tests measurable w.r.t. the residual map -- the map is lossy compression, and no claim is made against full-data tests. (5) The deployed CNN only ESTIMATES the Bayes score (finite corpus + architecture + optimization + distillation); the proposition covers the idealized rule, and the paper must not let 'optimal' migrate onto the deployed test in the abstract or conclusions. WHAT AUC = .5 CERTIFIES HERE: nothing. The measured AUC-.50 is a calibration diagnostic tied to Proposition 1 (indistinguishability of S-scores between the data and its bootstrap replicates); it carries no information about power or about closeness of the CNN to the Bayes score, and Prop 2 makes no use of it. The one genuine contact point between the measured facts and Prop 2 is deployment-matched training: it redefines the null class (hence P0) as the law of the map on bootstrap replicates, and Prop 2 applies verbatim with that P0 -- so matched training aligns the discrimination target with the deployed calibration without changing the optimality structure. This is stated as Remark 3 rather than smuggled into the theorem. ONE EDITORIAL CORRECTION TO THE TASK MEMO: the King point-optimal-testing citation is 1987 (Econometric Reviews 6(2), 169-218), not 1988.
```

### literature

```
VERIFIED CITATIONS AND WHAT EACH IS CITED FOR: (1) Neyman & Pearson (1933), Phil. Trans. R. Soc. A 231, 289-337 -- the NP lemma itself; cite once, historically. (2) Lehmann & Romano, Testing Statistical Hypotheses (3rd/4th ed., Springer), Theorem 3.2.1 -- the randomized NP lemma actually invoked in the proof. (3) Wald (1943), Trans. Amer. Math. Soc. 54, 426-482 -- weighted average power as the optimality criterion; the framing of the whole proposition. (4) King (1987), 'Towards a theory of point optimal testing', Econometric Reviews 6(2), 169-218, DOI 10.1080/07474938708800129 -- VERIFIED: year is 1987, NOT 1988 as in the task memo; point-optimal/power-envelope tradition. (5) Andrews & Ploberger (1994), Econometrica 62(6), 1383-1414 -- average-power-optimal tests built by integrating over a prior on the alternative; the closest classical template for 'library = prior'. (6) Lopez-Paz & Oquab (2017), 'Revisiting classifier two-sample tests', ICLR -- classifier-as-test-statistic precedent; cite in Remark 2. (7) Kim, Ramdas, Singh & Wasserman (2021), 'Classification accuracy as a proxy for two-sample testing', Ann. Statist. 49(1), 411-434, DOI 10.1214/20-AOS1962 -- VERIFIED; rate/regret transfer from classification to testing; cite in Remark 2. (8) Cranmer, Pavez & Louppe (2015), arXiv:1506.02169 -- calibrated discriminative classifiers approximate likelihood ratios (simulation-based inference); supports the estimation remark. (9) Friedman (2004), 'On multivariate goodness-of-fit and two-sample testing', PHYSTAT2003 proceedings, SLAC-PUB-10325 -- earliest ML-statistic GOF/two-sample precedent; optional one-line nod. (10) Hinton, Vinyals & Dean (2015), arXiv:1503.02531 -- distillation, cited only in Remark 2. SCOOP FLAGS: the mathematical content of Prop 2 (Bayes classifier = monotone in mixture LR = NP-optimal against the mixture) is folklore and is essentially present in Lopez-Paz & Oquab 2017, Kim et al. 2021, and Cranmer et al. 2015 for the TWO-SAMPLE / SBI settings -- we must NOT sell Prop 2 as a mathematical novelty. The defensible novelty lives in the composite-null GOF-for-logistic setting: the residual-map compression, the bootstrap-calibration validity condition (Prop 1), the measured AUC surrogate, and the deployment-matched-training alignment (Remark 3). Also worth a related-work sentence: Kubler et al. (2022) 'AutoML two-sample test' (arXiv:2206.08843) automates learned two-sample statistics -- adjacent, not overlapping our composite-null bootstrap setting. Sources checked: https://www.tandfonline.com/doi/abs/10.1080/07474938708800129 (King 1987), https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-1/Classification-accuracy-as-a-proxy-for-two-sample-testing/10.1214/20-AOS1962.pdf and https://arxiv.org/abs/1602.02210 (Kim et al. 2021).
```


## Fable's hand-verification note (2026-08-21)

Checked by direct calculation before any referee verdict landed:
1. **The counterexample (ex:marginal) is VALID.** Under P0, Z|V ~ N(cV,1); the estimator picks
   theta-hat = V and the bootstrap generates Z* ~ N(-Vc,1) — the SIGN FLIPS. So conditionally the
   observed score sits on the far side of its own bootstrap cloud with probability 1/2 as c grows,
   while the MARGINALS of S(D) and S(D*) are both (1/2)N(c,1)+(1/2)N(-c,1) — identical, AUC exactly
   1/2. Marginal blindness genuinely does not imply validity. This kills the naive version of our
   claim and the draft correctly ships the conditional version instead.
2. **pi* is exactly Uniform(0,1)** when G_t is atomless: conditionally on beta-hat = t,
   S(D*_0) ~ G_t so Gbar_t(S(D*_0)) ~ U(0,1); hence d_K(L(pi), Unif) = d_K(L(pi), L(pi*)) and the
   kappa <= 2(AUC*-1/2) chain follows from the standard TV-AUC relation. Structure sound.
3. **Consequence for the paper's empirics:** our fig-4 AUC (score scale) is the NECESSARY check;
   the SUFFICIENT certificate is the same AUC computed on the CALIBRATED scale (pi vs pi*) — and we
   can compute it from the deposited M2 blocks (pi = the exact-rank p-value itself, pi* = the rank
   of one held-out replicate among the rest). Opus: add this one cheap computation to S3.
