<!-- 16_THEORY_PACKAGE.md - FINAL, from workflow wf_c4eecd88-3e8 (2026-08-21): drafts ->
hostile referees -> repairs. All 14 attack repairs incorporated; claims shown false now appear
only in their weaker TRUE form with the failed hope stated in a remark. Fable's independent
hand-check of the counterexample and the pi*-uniformity argument PASSED (see 15_THEORY_DRAFTS.md).
OPUS: verify each proof line-by-line before submission - this is refereed by agents, not humans. -->

FINAL THEORY PACKAGE — Deep-GOF small-sample paper (Propositions 1 and 2, all attack repairs incorporated)

Repairs applied: P1-R1 (arbitrary measurable estimator; MLE reserved for the impossibility remark), P1-R2 (Remark on direction rewritten; (ii)(c) conceded degenerate in the flagship model), P1-R3 (independent-draws fix in the Corollary proof + empirical-protocol sentence), P1-R4 ("degenerately 1" replaced; pair defined), P1-R5 (π = π* = 1 on the failure event; extended-line lemma), P1-cites (Robins–van der Vaart–Ventura, Bayarri–Berger, Davidson–MacKinnon, Davison–Hinkley + positioning remark), P1-body (both overclaims removed); P2-1 (μ = P₀+Q̄; a.e.-uniqueness weakened), P2-2 (Markov-kernel standing assumption), P2-3 (nonnegativity + boundary cases in the proof), P2-4 (deployed test not map-measurable, stated in full), P2-5 (cost of the redefined P₀ named), P2-6 (both body-paragraph fixes), P2-7 (NPLM + conformal-C2ST cites). Every claim an attack showed false now appears only in its weaker true form, with the failed hope stated in a remark.

==================================================================
PROPOSITION 1 — STATEMENT (body/supplement-ready LaTeX)
==================================================================

```latex
%======================================================================
% PROPOSITION 1 -- Validity of the deployed bootstrap calibration
%======================================================================

\paragraph{Setting.}
Let $\{P_\beta:\beta\in\Theta\}$, $\Theta\subseteq\mathbb{R}^p$, be the null model for a dataset
$D=(y,X)$ with the design $X$ held fixed; under $H_0$, $D\sim P_{\beta_0}$ for some unknown
$\beta_0\in\Theta$. Let $\hat\beta(\cdot)$ be a \emph{fixed measurable estimator} of $\beta$ ---
in our implementation the maximum-likelihood estimator, but no property of maximum likelihood is
used before Remark~\ref{rem:impossible} --- and let $S(\cdot)$ be a \emph{frozen} measurable
statistic with values in $\overline{\mathbb{R}}=\mathbb{R}\cup\{+\infty\}$: a fixed function of
the dataset alone, chosen before $D$ is observed (in our case $S$ contains its own refit). Given
$D$, draw $D^{*}_1,\dots,D^{*}_B$ i.i.d.\ from $P_{\hat\beta(D)}$ (same $X$) and set
\begin{equation}\label{eq:pval}
  p \;=\; \frac{1+\#\{b\le B: S(D^{*}_b)\ge S(D)\}}{B+1}.
\end{equation}
For $t\in\Theta$ let $G_t$ be the law of $S(D^{*})$ under $D^{*}\sim P_t$, write
$\bar G_t(s)=P_t\{S(D^{*})\ge s\}$ (survival on $\overline{\mathbb{R}}$), and define the
\emph{calibration functional}
\begin{equation}\label{eq:pi}
  \pi(D)\;=\;\bar G_{\hat\beta(D)}\bigl(S(D)\bigr),
\end{equation}
the ideal ($B=\infty$) bootstrap $p$-value, a function of $D$ alone.

\emph{Failure conventions.} If the estimator fails to exist at the observed data we set $p=1$,
$\pi(D)=1$ and $\pi^{*}=1$ (see part (iii)); if a bootstrap refit fails we set
$S(D^{*}_b)=+\infty$. Each convention can only increase $p$, keeps
$\mathbb{P}\{p\le\alpha\}=\mathbb{E}\,h_k(\pi)$ intact (proof of
Proposition~\ref{prop:validity}), and preserves the super-uniformity of $\pi^{*}$; all
inequalities below hold with the conventions in force, while the equality statements assume the
relevant fits exist almost surely.

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
$\mathcal{L}_{\beta_0}\bigl(S(D)\,\big|\,\hat\beta(D)=t\bigr)=G_t$;
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
bootstrap draw $D^{*}_0\mid D\sim P_{\hat\beta(D)}$ (with $\pi^{*}:=1$ on the estimator-failure
event), $d_{\mathrm K}$ is the Kolmogorov distance, and $\mathrm{AUC}^{\star}$ is the supremum,
over all measurable classifiers, of the AUC for separating \emph{independent} draws of $\pi(D)$
from \emph{independent} draws of $\pi^{*}$. The law $\mathcal{L}(\pi^{*})$ is super-uniform by
construction, with no assumptions on $S$.
\end{enumerate}
\end{proposition}

\begin{corollary}[What $\mathrm{AUC}=\tfrac12$ does and does not certify]\label{cor:auc}
Let $\mathrm{AUC}^{\star}_S$ denote the optimal AUC of any classifier reading the score alone,
i.e.\ separating independent draws from $\mathcal{L}_{\beta_0}\bigl(S(D)\bigr)$ and from
$\mathcal{L}\bigl(S(D^{*}_1)\bigr)$.
\begin{enumerate}
\item[(i)] \textup{(Necessity.)} $\mathrm{AUC}^{\star}_S=\tfrac12$ if and only if
$\mathcal{L}\bigl(S(D)\bigr)=\mathcal{L}\bigl(S(D^{*}_1)\bigr)$. In particular, exchangeability of
$\bigl(S(D),S(D^{*}_1),\dots,S(D^{*}_B)\bigr)$ --- e.g.\ under
Proposition~\ref{prop:validity}(ii)(a) or (ii)(b) --- forces $\mathrm{AUC}^{\star}_S=\tfrac12$;
any classifier achieving $\mathrm{AUC}>\tfrac12$ therefore refutes exchangeability, and
quantitatively every classifier $f$ satisfies
$\mathrm{AUC}(f)-\tfrac12\le\mathrm{TV}\bigl(\mathcal{L}(S(D)),\mathcal{L}(S(D^{*}))\bigr)$.
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
Let $\Theta=\{-1,0,1\}$ and $D=(V,Z)\in\{-1,1\}\times\mathbb{R}$. Under $P_0$:
$V\sim\mathrm{Unif}\{\pm1\}$ and $Z\mid V\sim N(cV,1)$; under $P_{\pm1}$:
$V\sim\mathrm{Unif}\{\pm1\}$ independently of $Z\sim N(\mp c,1)$. Take the (measurable, non-ML)
estimator $\hat\theta(D)=V$ --- the Setting requires measurability only --- and $S(D)=Z$, with
true value $\theta_0=0$. Then $S(D)$ and $S(D^{*})$ have the identical marginal law
$\tfrac12 N(c,1)+\tfrac12 N(-c,1)$, so $\mathrm{AUC}^{\star}_S=\tfrac12$ exactly; yet for every
$\alpha\ge 1/(B+1)$, $\liminf_{c\to\infty}\mathbb{P}_{0}\{p\le\alpha\}\ge\tfrac12$. Conditionally
on $\hat\theta$ the observed score sits on the wrong side of its own bootstrap cloud half the
time; the failure is invisible marginally and visible only relative to the fitted model.
\end{example}

\begin{remark}[Relation to exact Monte Carlo tests]\label{rem:mc}
If $\beta_0$ were known and the replicates drawn from $P_{\beta_0}$, the vector
$\bigl(S(D),S(D^{*}_1),\dots\bigr)$ would be i.i.d.\ and \eqref{eq:pval} an exact Monte Carlo test
for \emph{any} $S$ (Barnard, 1963; Hope, 1968; Besag and Clifford, 1989; the finite-$B$ rank
$p$-value goes back to Dwass, 1957). The deployed procedure is instead a ``local'' Monte Carlo
test at a point estimate in the sense of Dufour (2006), and exchangeability genuinely fails in
general: $D\sim P_{\beta_0}$ while $D^{*}_b\mid D\sim P_{\hat\beta}$, so the replicates are
i.i.d.\ only \emph{conditionally} on $\hat\beta$, and equality of the marginal laws of $S(D)$ and
$S(D^{*})$ is strictly weaker than condition (ii)(b) (Example~\ref{ex:marginal}). Dufour restores
exactness by maximizing the $p$-value over the nuisance space; Proposition~\ref{prop:validity}
instead certifies the un-maximized procedure that is actually run.
\end{remark}

\begin{remark}[Relation to plug-in $p$-value asymptotics]\label{rem:rvv}
Robins, van der Vaart and Ventura (2000) characterize the asymptotic law of plug-in
(parametric-bootstrap) $p$-values in composite null models: $\pi$ is non-uniform in the limit in
general, and conservative when the asymptotic mean of the statistic depends on the nuisance
parameter; Bayarri and Berger (2000) is the companion analysis. Davidson and MacKinnon (1999)
study the asymptotic order of bootstrap size distortion, and $\pi^{*}$ is precisely a prepivoted
root in the sense of Beran (1988). Proposition~\ref{prop:validity} is the finite-sample,
finite-$B$ counterpart of this line: rather than characterizing the limit law of $\pi$, part
(iii) bounds the size distortion at every $n$ and $B$ by
$d_{\mathrm K}(\mathcal{L}(\pi),\mathcal{L}(\pi^{*}))$, a distance between calibrated $p$-value
laws that is estimable by the inner (double) bootstrap (Davison and Hinkley, 1997, \S4.5).
Deployment-matched training can be read, in their terms, as pushing the learned statistic toward
pivotality; the certificate here does not rely on that reading.
\end{remark}

\begin{remark}[Exact similarity is impossible here; approximate validity is the right target]\label{rem:impossible}
Specialize now to the logistic model with canonical link and take $\hat\beta$ to be the MLE.
Then $T=X^{\top}y$ is a complete sufficient statistic (full-rank exponential family with natural
parameter space $\mathbb{R}^p$), the likelihood depends on $y$ only through $T$, and on the event
that the MLE exists it is in bijection with $\hat\beta$. If
$\ker(X^{\top})\cap\{-1,0,1\}^n=\{0\}$ --- which holds for Lebesgue-almost-every design with at
least one continuously distributed covariate --- the fibers $\{y:X^{\top}y=t\}$ are singletons,
so conditioning on $\hat\beta$ determines $y$; this is precisely why co-sufficient sampling
(Engen and Lilleg{\aa}rd, 1997; Lindqvist and Taraldsen, 2005) is powerless in this model (Barber
and Janson, 2022). Three consequences. First, condition (ii)(b) can then hold only degenerately
and no statistic $S$ is exactly pivotal. Second, the same degeneracy afflicts (ii)(c):
conditionally on $\hat\beta$ the observed score is a point mass at its realized value, so
conditional dominance forces $G_t$ to charge $[s(t),\infty]$ with probability one --- (ii)(b) and
(ii)(c) alike are idealizations in this model, available only in models where the conditioning is
non-degenerate. Third, by Neyman structure (Lehmann and Romano, 2005, \S4.3), \emph{every} test
with exact size $\alpha$ for all $\beta$ satisfies $\mathbb{E}(\phi\mid T)=\alpha$ a.s., hence
$\phi(y)=\alpha$ a.s.: the only exactly similar tests ignore the data. Exact finite-sample
similarity is therefore unattainable for any $S$, learned or classical; the operative guarantee
is part (iii), and an estimated size of $.052$ at nominal $.05$ is what the theory predicts as
the achievable target, not a defect of the learned statistic.
\end{remark}

\begin{remark}[Direction of the distortion and the measured facts]\label{rem:direction}
Part (ii)(c) is one-sided: in models where the conditional law of $S$ given $\hat\beta$ is
non-degenerate, conditional dominance of the bootstrap replicates yields a valid but conservative
test, paying in power. In the flagship logistic model (ii)(c) is degenerate
(Remark~\ref{rem:impossible}), so we do not claim it holds there; what the measured facts
instantiate is the \emph{marginal} conservative direction --- $F_\pi$ super-uniform, $\kappa=0$
in part (iii). Training $S$ with true nulls as the null class produced exactly this regime
(score-scale AUC $.62$, deployed size $.006$): the mismatch surfaced as deflation, not inflation.
Deployment-matched training, whose null class reproduces the fit--resample--refit pipeline, gave
score-scale AUC $.5011$ and deployed size $.052$. The monotone relation between the score-scale
AUC and both deployed size and power across the fifteen trained rules is \emph{consistent with}
the chain in part (iii) but is not implied by it, since the chain lives on the calibrated scale
while the measurement lives on the score scale; we report it as an empirical regularity.
\end{remark}

\begin{remark}[Which classifier to run]\label{rem:whichclassifier}
The three scales are not interchangeable. A classifier on raw scores tests the necessary marginal
condition of Corollary~\ref{cor:auc}(i); because the AUC in Lemma~\ref{lem:auc} is defined for
independent draws from the two laws, the empirical protocol must score the two classes on
statistics computed from \emph{independent simulated datasets} (one score per dataset), never on
within-dataset $(S(D),S(D^{*}))$ couples, whose dependence enters no statement here; our pipeline
complies. A classifier on pairs --- class one $(S(D),\hat\beta(D))$ from true null draws, class
two $(S(D^{*}),\hat\beta(D))$ with $D^{*}$ a bootstrap replicate of the same $D$ --- is the
natural conditional test, but in the continuous-design logistic model its oracle AUC approaches
$1$ up to the bootstrap atom at the observed configuration (the model is discrete, so
$\mathbb{P}\{S(D^{*})=S(D)\mid\hat\beta\}>0$; Remark~\ref{rem:impossible}), so at the oracle
level it certifies nothing. The scale that matches the deployed test is the calibrated one:
$\pi$ versus $\pi^{*}$, both estimable by an inner (double) bootstrap (Davison and Hinkley, 1997,
\S4.5), and it is on this scale that part (iii) converts an oracle AUC of $\tfrac12$ into a
validity certificate. An empirical classifier only \emph{lower}-bounds the oracle AUC, so a
measured AUC near $\tfrac12$ is an estimate that the slack $2(\mathrm{AUC}^{\star}-\tfrac12)$ is
small, never a proof.
\end{remark}

\begin{remark}[Honest accounting: what was hoped, what is true]\label{rem:honest1}
Three natural hopes fail, and we record why. \emph{Hope 1: the measured score-scale AUC near
$\tfrac12$ certifies validity.} False: Example~\ref{ex:marginal} has optimal score-scale AUC
exactly $\tfrac12$ and size near $\tfrac12$; the score-scale AUC is necessary only, and the
sufficiency in Corollary~\ref{cor:auc}(iii) holds on the calibrated scale for the \emph{oracle}
classifier, which any trained classifier merely lower-bounds. The direct evidence that the slack
is small is the deployed size itself, supplemented by an inner-bootstrap classifier on $\pi$
versus $\pi^{*}$. \emph{Hope 2: condition (ii)(c) explains the conservative $.006$ regime.}
In the flagship model (ii)(c) is degenerate (Remark~\ref{rem:impossible}); only the weaker,
marginal statement survives, and that is all Remark~\ref{rem:direction} claims. \emph{Hope 3:
exact validity is achievable by a better statistic.} Impossible by Neyman structure
(Remark~\ref{rem:impossible}); part (iii) is the operative guarantee. Standing assumptions,
none hidden: $S$ frozen and independent of $D$; fixed design, bootstrap at the same $X$;
i.i.d.\ replicates given $D$; failure conventions as stated (they only inflate $p$); all laws
discrete, so results are stated in the conservative super-uniform direction, which needs no
continuity; validity is marginal over data and bootstrap randomness; no asymptotic claims are
made anywhere.
\end{remark}
```

==================================================================
PROPOSITION 1 — SUPPLEMENT PROOFS (corrected LaTeX)
==================================================================

```latex
%======================================================================
% SUPPLEMENT: proofs for Proposition 1, Corollary 1, Example 1,
% and the claims of Remark rem:impossible
%======================================================================

We use two elementary lemmas. Lemma~\ref{lem:pit} is stated on the extended line
$\overline{\mathbb{R}}=\mathbb{R}\cup\{+\infty\}$ because the bootstrap-failure convention gives
$S$ an atom at $+\infty$; the proof is verbatim the same as on $\mathbb{R}$.

\begin{lemma}[Survival probability integral transform, extended line]\label{lem:pit}
Let $W\sim Q$ on $\overline{\mathbb{R}}$ and $\bar Q(s)=Q([s,\infty])$. Then
$\mathbb{P}\{\bar Q(W)\le u\}\le u$ for all $u\in(0,1)$, with equality for all $u$ if $Q$ is
atomless.
\end{lemma}

\begin{proof}
Fix $u$ and let $A_u=\{s:\bar Q(s)\le u\}$. Since $\bar Q$ is nonincreasing, $A_u$ is an up-set in
$\overline{\mathbb{R}}$: $A_u=[s_u,\infty]$ or $(s_u,\infty]$ with $s_u=\inf A_u$ (the case
$A_u=\emptyset$ is trivial). If $s_u\in A_u$ then
$\mathbb{P}\{W\in A_u\}\le\mathbb{P}\{W\ge s_u\}=\bar Q(s_u)\le u$. If $s_u\notin A_u$ then
$\mathbb{P}\{W\in A_u\}=\mathbb{P}\{W>s_u\}=\lim_{s\downarrow s_u}\bar Q(s)\le u$ because
$\bar Q(s)\le u$ for every $s>s_u$. Since $\{\bar Q(W)\le u\}=\{W\in A_u\}$, the bound follows.
If $Q$ is atomless (hence supported on $\mathbb{R}$) then $\bar Q(W)=1-F_Q(W)$ a.s.\ with
$F_Q(W)\sim\mathrm{Unif}(0,1)$ by the probability integral transform, giving equality.
\end{proof}

\begin{lemma}[AUC and total variation]\label{lem:auc}
Let $P,Q$ be laws on a Polish space, and for a measurable $f$ let
$\mathrm{AUC}(f)=\mathbb{P}\{f(W)>f(W')\}+\tfrac12\mathbb{P}\{f(W)=f(W')\}$ with $W\sim P$,
$W'\sim Q$ \emph{independent}. Then
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
then $\mathbb{E}\,g(W'',W')=0$ by exchangeability of the independent identically distributed pair
$(W'',W')$. For each fixed $w'$,
$\bigl|\int g(w,w')\,(P-Q)(dw)\bigr|\le \|g\|_\infty\,\|P-Q\|_{1}\le 2\,\mathrm{TV}(P,Q)$, hence
$\mathbb{E}\,g(W,W')=\mathbb{E}\,g(W,W')-\mathbb{E}\,g(W'',W')\le 2\,\mathrm{TV}(P,Q)$. The final
claim combines the two bounds.
\end{proof}

\begin{proof}[Proof of Proposition~\ref{prop:validity}]
\emph{Reduction to the calibration functional.} On the event that $\hat\beta(D)$ exists, the
events $\{S(D^{*}_b)\ge S(D)\}$, $b=1,\dots,B$, are conditionally independent given $D$ with
common probability $\bar G_{\hat\beta(D)}\bigl(S(D)\bigr)=\pi(D)$, because the $D^{*}_b$ are
i.i.d.\ $P_{\hat\beta(D)}$ given $D$ and $S(D)$ is $D$-measurable. Hence
$N:=\#\{b:S(D^{*}_b)\ge S(D)\}$ satisfies $N\mid D\sim\mathrm{Bin}(B,\pi(D))$. Since $N$ is an
integer, $p\le\alpha \iff N\le k-1$ with $k:=\lfloor\alpha(B+1)\rfloor$, and therefore
\begin{equation}\label{eq:reduction}
  \mathbb{P}\{p\le\alpha\}=\mathbb{E}\bigl[h_k\bigl(\pi(D)\bigr)\bigr],
  \qquad
  h_k(x):=\mathbb{P}\{\mathrm{Bin}(B,x)\le k-1\}
  =\sum_{j=0}^{k-1}\binom{B}{j}x^{j}(1-x)^{B-j}.
\end{equation}
On the estimator-failure event the conventions give $p=1$ and $\pi=1$, and
$\mathbf{1}\{p\le\alpha\}=0=h_k(1)$ for every $\alpha<1$ (since $k\le B$), so
\eqref{eq:reduction} holds on that event as well. The polynomial $h_k$ maps $[0,1]$ to $[0,1]$
and is nonincreasing: $h_k'(x)=-B\binom{B-1}{k-1}x^{k-1}(1-x)^{B-k}\le0$ (the binomial--beta
identity; for $k=0$, $h_0\equiv0$ and all claims are trivial). For $U\sim\mathrm{Unif}(0,1)$,
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
\emph{Second inequality.} On the event that $\hat\beta(D)$ exists, conditionally on $D$ the extra
draw satisfies $S(D^{*}_0)\sim G_{\hat\beta(D)}$ exactly, so Lemma~\ref{lem:pit} applied
conditionally gives $\mathbb{P}\{\pi^{*}\le u\mid D\}\le u$; on the failure event $\pi^{*}=1$, so
$\mathbb{P}\{\pi^{*}\le u\mid D\}=0\le u$ for $u<1$. Hence $F_{\pi^{*}}(u)\le u$ for all
$u\in(0,1)$: $\pi^{*}$ is super-uniform \emph{by construction}. Therefore, pointwise in $u$,
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
$\mathcal{L}(S(D))=\mathcal{L}(S(D^{*}_1))$. Conversely, suppose the two marginal laws coincide
and let $W\sim\mathcal{L}(S(D))$ and $W'\sim\mathcal{L}(S(D^{*}_1))$ be \emph{independent} draws,
as in the definition of $\mathrm{AUC}$ in Lemma~\ref{lem:auc}; then $(W,W')$ is an i.i.d.\ pair,
hence $(f(W),f(W'))$ is exchangeable for every $f$ and $\mathrm{AUC}(f)=\tfrac12$. (The claim
concerns independent draws from the two marginal laws only; a within-dataset couple
$(S(D),S(D^{*}))$ with equal marginals need not be exchangeable, and no statement here is made
about it --- correspondingly, the empirical AUC must be computed from scores on independent
datasets, one per class member.) The quantitative bound is Lemma~\ref{lem:auc}. (ii) is
Example~\ref{ex:marginal}, verified below. (iii) is Proposition~\ref{prop:validity}(iii)
restated.
\end{proof}

\begin{proof}[Details for Example~\ref{ex:marginal}]
Note first that $\hat\theta(D)=V$ is a measurable map of $D$, as the Setting requires; it is not
the MLE of the three-point model, and nothing in Proposition~\ref{prop:validity} or
Corollary~\ref{cor:auc} uses maximum likelihood. Under $P_0$ write $Z=cV+\varepsilon$,
$\varepsilon\sim N(0,1)\perp V$. The bootstrap law of the score at $t=\pm1$ is $G_t=N(-ct,1)$, so
the marginal law of $S(D^{*})$ is $\tfrac12 N(-c,1)+\tfrac12 N(c,1)$, equal to the marginal law of
$S(D)=Z$; by Lemma~\ref{lem:auc}, $\mathrm{AUC}^{\star}_S=\tfrac12$ exactly. Next,
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
classes with AUC tending to $1$ as $c\to\infty$.
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
union over the finitely many $v$ remains null. On singleton fibers, conditioning on $\hat\beta$
determines $y$, so $\mathcal{L}(S(D)\mid\hat\beta=t)$ is a point mass $\delta_{s(t)}$: condition
(ii)(b) requires $G_t=\delta_{s(t)}$ and condition (ii)(c) requires
$G_t\bigl([s(t),\infty]\bigr)=1$ for a.e.\ $t$ --- both degenerate, as claimed. Now let
$\phi(y)\in[0,1]$ be the rejection probability of any (possibly randomized) test with
$\mathbb{E}_\beta\phi=\alpha$ for all $\beta\in\mathbb{R}^p$. By sufficiency
$\mathbb{E}(\phi\mid T)$ does not depend on $\beta$, and by completeness applied to
$\mathbb{E}(\phi\mid T)-\alpha$ we get $\mathbb{E}(\phi\mid T)=\alpha$ a.s.\ (Neyman structure;
Lehmann and Romano, 2005, Theorem~4.3.2). On singleton fibers $\mathbb{E}(\phi\mid T)=\phi$, so
$\phi(y)=\alpha$ a.s.: the test rejects with probability $\alpha$ regardless of the data. Applied
to the deployed test at any attainable level $k/(B+1)$, this shows no choice of $S$ makes it
exactly similar and nontrivial.
\end{proof}
```

==================================================================
PROPOSITION 1 — BODY PARAGRAPH (corrected)
==================================================================

```latex
Proposition~\ref{prop:validity} replaces the folklore justification of the deployed calibration
with exact finite-sample statements, and it separates what is achievable from what is not. The
classical route to exactness --- exchangeability of $(S(D),S(D^*_1),\dots,S(D^*_B))$, as in Monte
Carlo tests with known null --- holds if $S$ is pivotal or, more sharply, conditionally pivotal
given the fitted model; but the replicates are i.i.d.\ only conditionally on $\hat\beta$ while
$D$ is drawn at the true $\beta_0$, so equality of the marginal score distributions is strictly
weaker than what validity needs (Example~\ref{ex:marginal} attains AUC exactly $\tfrac12$ on the
score scale with rejection probability near $\tfrac12$). Worse, in the continuous-design logistic
model no nontrivial exactly valid test exists at all: $X^\top y$ is complete sufficient with
almost-surely singleton fibers, so by Neyman structure every exactly similar test ignores the
data (Remark~\ref{rem:impossible}; Barber and Janson, 2022) --- a deployed size of $.052$ at
nominal $.05$ is the theoretically achievable target, not a defect. The operative guarantee is
therefore quantitative: at every level and every $B$, the size distortion is at most the
Kolmogorov distance between the calibrated score $\pi(D)$ and its by-construction super-uniform
bootstrap counterpart $\pi^*$, hence at most $2(\mathrm{AUC}^\star-\tfrac12)$ for the best
classifier separating the two on the $p$-value scale. This delineates exactly what the measured
AUCs certify: the score-scale AUC of $.62$ refutes exchangeability outright; the $.5011$ figure
passes the necessary marginal check, while the direct evidence that the validity slack is small
is the deployed size of $.052$ itself and, on the theorem's own scale, an inner-bootstrap
classifier on $\pi$ versus $\pi^*$; and the monotone AUC--size relation across the fifteen
trained rules is consistent with, though not implied by, this chain.
```

==================================================================
PROPOSITION 2 — STATEMENT (corrected LaTeX)
==================================================================

```latex
% ---------- Setting (shared notation, to precede Proposition 2) ----------
% M : data -> T \subset R^{KxK} is the fixed residual-map transform (Section 2).

\paragraph{Standing assumption for Proposition~\ref{prop:avgpower}.}
Let $(\mathcal{G},\mathscr{G})$ be a measurable space of departure functions and let
$g\mapsto Q_g$ be a Markov kernel from $(\mathcal{G},\mathscr{G})$ to the map space
$\mathcal{T}\subseteq\mathbb{R}^{K\times K}$, i.e.\ $g\mapsto Q_g(A)$ is
$\mathscr{G}$-measurable for every measurable $A\subseteq\mathcal{T}$; let $\Pi$ be a probability
measure on $(\mathcal{G},\mathscr{G})$, so that $\bar Q:=\int Q_g\,d\Pi(g)$ is a well-defined
probability measure on $\mathcal{T}$. In our setting this is automatic: the training library is a
finitely parameterized simulator, and $Q_g(A)$ is jointly measurable in its parameters. Let $P_0$
be the null-class law of the residual map $M(D)$ under the training generator. Throughout we take
the dominating measure to be $\mu:=P_0+\bar Q$ --- a choice, not an assumption --- and write
$p_0=dP_0/d\mu$, $\bar q=d\bar Q/d\mu$, so that $p_0+\bar q=1$ $\mu$-a.e.

\begin{proposition}[$\Pi$-average-power optimality of the idealized rule]
\label{prop:avgpower}
With the standing assumption above, define the balanced-class Bayes posterior
\[
\eta^\star(t) \;=\; \frac{\bar q(t)}{p_0(t)+\bar q(t)},
\qquad t\in\mathcal{T},
\]
and the mixture likelihood ratio $\Lambda(t)=\bar q(t)/p_0(t)$, with the convention
$\Lambda(t)=\infty$ on $\{p_0=0<\bar q\}$.
\begin{enumerate}
\item[(i)] $\eta^\star$ minimizes the balanced-class population cross-entropy risk over all
measurable $f:\mathcal{T}\to[0,1]$, and any minimizer agrees with $\eta^\star$
$(P_0+\bar Q)$-a.e.; moreover $\eta^\star=\Lambda/(1+\Lambda)$ is a strictly increasing transform
of $\Lambda$ on $[0,\infty]$.
\item[(ii)] Fix $\alpha\in(0,1)$. Any test of the form
$\varphi^\star(t)=\mathbf{1}\{\eta^\star(t)>c\}+\gamma\,\mathbf{1}\{\eta^\star(t)=c\}$,
with $c\in[0,1]$ and $\gamma\in[0,1]$ chosen so that
$\mathbb{E}_{P_0}\varphi^\star(M)=\alpha$, maximizes the $\Pi$-average power
\[
\int \mathbb{E}_{Q_g}\!\bigl[\varphi(M)\bigr]\,d\Pi(g)
\]
over all (possibly randomized) tests $\varphi:\mathcal{T}\to[0,1]$ satisfying
$\mathbb{E}_{P_0}\varphi(M)\le\alpha$. If moreover $\Lambda(M)$ has a continuous distribution
under $P_0$, the non-randomized test $\mathbf{1}\{\eta^\star(M)>c_\alpha\}$, with $c_\alpha$ the
$(1-\alpha)$-quantile of $\eta^\star(M)$ under $P_0$, attains the maximum exactly.
\end{enumerate}
\end{proposition}
```

==================================================================
PROPOSITION 2 — PROOF + REMARKS (corrected LaTeX)
==================================================================

```latex
\begin{proof}
(i) With balanced class prior, the population cross-entropy risk of a measurable
$f:\mathcal{T}\to[0,1]$ is
\[
R(f)\;=\;-\tfrac12\,\mathbb{E}_{P_0}\log\bigl(1-f(M)\bigr)
        -\tfrac12\,\mathbb{E}_{\bar Q}\log f(M)
      \;=\;\int \Bigl[-\tfrac12\,p_0(t)\log\bigl(1-f(t)\bigr)
        -\tfrac12\,\bar q(t)\log f(t)\Bigr]\,\mu(dt),
\]
whose integrand is nonnegative (values in $[0,\infty]$), so $R$ is well defined in $[0,\infty]$
and may be minimized pointwise in $t$. Since $\mu=P_0+\bar Q$, we have $p_0(t)+\bar q(t)=1$ for
$\mu$-a.e.\ $t$, so there is no $\mu$-non-null set where both densities vanish. At a point with
$p_0(t)>0$ and $\bar q(t)>0$ the integrand is strictly convex in $f(t)\in(0,1)$, tends to
$+\infty$ at both endpoints, and is uniquely minimized at
$f(t)=\bar q(t)/\{p_0(t)+\bar q(t)\}=\eta^\star(t)$. At a boundary point the minimizer is at the
boundary and still equals $\eta^\star$: if $\bar q(t)=0<p_0(t)$ the integrand
$-\tfrac12 p_0(t)\log(1-f(t))$ is strictly increasing in $f(t)$, minimized uniquely at
$f(t)=0=\eta^\star(t)$; symmetrically, if $p_0(t)=0<\bar q(t)$ the unique minimizer is
$f(t)=1=\eta^\star(t)$. Hence $\eta^\star$ minimizes $R$ and any minimizer agrees with
$\eta^\star$ $\mu$-a.e., i.e.\ $(P_0+\bar Q)$-a.e.\ (this is the standard population property of
the logarithmic scoring rule). Since $x\mapsto x/(1+x)$ is strictly increasing on $[0,\infty]$
(with value $1$ at $\infty$), $\eta^\star=\Lambda/(1+\Lambda)$ is a strictly increasing transform
of $\Lambda$, proving (i).

(ii) For any test $\varphi:\mathcal{T}\to[0,1]$, the standing kernel assumption makes
$(g,t)\mapsto\varphi(t)$ jointly integrable structure available, and Tonelli's theorem (the
integrand lies in $[0,1]$) gives
\[
\int \mathbb{E}_{Q_g}\!\bigl[\varphi(M)\bigr]\,d\Pi(g)
=\int\!\!\int \varphi\,dQ_g\,d\Pi(g)
=\mathbb{E}_{\bar Q}\bigl[\varphi(M)\bigr].
\]
Thus maximizing $\Pi$-average power subject to $\mathbb{E}_{P_0}\varphi(M)\le\alpha$ is exactly
the problem of testing the simple hypothesis $P_0$ against the simple alternative $\bar Q$ at
level $\alpha$. By the Neyman--Pearson lemma \citep[Theorem~3.2.1]{LehmannRomano}, a test is most
powerful at level $\alpha$ if (and, up to $(P_0+\bar Q)$-null modification off the boundary set,
only if) it has the form $\mathbf{1}\{\Lambda(t)>k\}+\gamma'\,\mathbf{1}\{\Lambda(t)=k\}$ with
size exactly $\alpha$. By the strict monotonicity established in (i),
$\{\eta^\star>c\}=\{\Lambda>k\}$ and $\{\eta^\star=c\}=\{\Lambda=k\}$ for $c=k/(1+k)$, so
$\varphi^\star$ is of Neyman--Pearson form with size $\alpha$ and therefore maximizes
$\mathbb{E}_{\bar Q}\varphi(M)$, which equals the $\Pi$-average power. Finally, if $\Lambda(M)$
is continuously distributed under $P_0$, then so is $\eta^\star(M)$, the boundary set
$\{\eta^\star(M)=c_\alpha\}$ has $P_0$-probability zero, and the non-randomized test
$\mathbf{1}\{\eta^\star(M)>c_\alpha\}$ has size $\alpha$ and is of Neyman--Pearson form, hence
attains the maximum.
\end{proof}

% ---------- Remarks ----------

\begin{remark}[Scope of the optimality claim; the deployed test is outside the class]\label{rem:scope}
Proposition~\ref{prop:avgpower} is an optimality statement \emph{within the class of tests
measurable with respect to the residual map} $M(D)$; the map is a lossy compression of $D$, and
no claim is made against tests that use the full data. Moreover, $P_0$ is the training-library
mixture null, so the size constraint is an \emph{average} size over the library's null
configurations, in the weighted-power sense of \citet{Wald1943}. Most importantly, the
\emph{deployed} test --- reject when $S(D)$ exceeds the bootstrap quantile of
$S(D^{*}_1),\dots,S(D^{*}_B)$ --- is \emph{not} map-measurable: its critical value depends on $D$
through $\hat\beta$ and the resampling, not through $M(D)$ alone, so the deployed test lies
outside the class over which Proposition~\ref{prop:avgpower} optimizes, and no formal optimality
transfers to it. What the proposition certifies is that the \emph{ordering} of the map space
induced by $\eta^\star$ is the $\Pi$-average-power-optimal ordering; the deployed procedure
composes that ordering with a different, conditionally calibrated cutoff, whose finite-sample
size properties are the business of Proposition~\ref{prop:validity}, not of this proposition.
\end{remark}

\begin{remark}[Idealization gap I: the network estimates the Bayes score]\label{rem:estimation}
The deployed statistic $S$ is a convolutional network fitted by empirical cross-entropy
minimization on a finite simulated corpus, followed by distillation into a single frozen network
\citep{Hinton2015}; it therefore \emph{estimates} $\eta^\star$ and inherits three discrepancies:
approximation error of the architecture class, estimation error from the finite corpus, and
optimization error. Because the test is invariant to strictly increasing transforms of the score,
only the \emph{ordering} induced by $S$ matters, and the $\Pi$-average power of the deployed
ordering degrades continuously with its rank disagreement from $\Lambda$; quantitative transfers
from classification risk to testing performance are available in the two-sample setting
\citep{Kim2021,LopezPaz2017}, the use of discriminative classifiers as likelihood-ratio
estimators is standard in simulation-based inference \citep{Cranmer2015}, and the decoupling of
finite-sample validity from classifier quality --- validity by construction, power degrading with
classifier error --- has a two-sample analogue in conformal classifier tests
\citep[arXiv:2507.17026]{ConformalC2ST}; learned Neyman--Pearson statistics for goodness of fit
also appear in the physics literature \citep[arXiv:2305.14137]{NPLM}. We make no claim that the
finite-corpus network attains the oracle power of Proposition~\ref{prop:avgpower}.
\end{remark}

\begin{remark}[Idealization gap II: the prior is the training library]\label{rem:prior}
The optimality criterion is average power with respect to the library-induced measure $\Pi$,
exactly the weighted-power framework of \citet{Wald1943}, \citet{King1987}, and
\citet{AndrewsPloberger1994}: choosing the training library is choosing where power is spent.
Against departures outside (or heavily down-weighted by) the support of $\Pi$,
Proposition~\ref{prop:avgpower} confers no guarantee, and a different library yields a different
optimal ordering. Finally, under the deployment-matched training of Section~\ref{sec:matched} the
null class is the law of the map on bootstrap replicates rather than on true null draws;
Proposition~\ref{prop:avgpower} applies verbatim with $P_0$ redefined as that law --- with the
caveat that ``size'' in the proposition then refers to that reference law, which is not the null
hypothesis of any user's data-generating configuration; actual size control is delivered
separately, by the bootstrap calibration under the conditions of
Proposition~\ref{prop:validity}. Matched training thus changes the reference null the score
discriminates against --- which is precisely what aligns the score with the deployed bootstrap
calibration --- without affecting the form of the optimality statement.
\end{remark}

\begin{remark}[Honest accounting: what was hoped, what is true]\label{rem:honest2}
Two natural strengthenings fail. \emph{Hope 1: $\eta^\star$ is the $\mu$-a.e.\ unique minimizer
for an arbitrary dominating measure.} False: off the support of $P_0+\bar Q$ a minimizer is
arbitrary; uniqueness holds exactly $(P_0+\bar Q)$-a.e., which is how (i) is stated (and taking
$\mu=P_0+\bar Q$ makes the density assumption a choice rather than a hypothesis). \emph{Hope 2:
the optimality of Proposition~\ref{prop:avgpower} transfers to the deployed test.} It does not:
the deployed test is not map-measurable (Remark~\ref{rem:scope}), and the proposition certifies
only the ordering that the trained score aims at. Separately, the measured deployment-matched
AUC $\approx.50$ is a calibration diagnostic belonging to Proposition~\ref{prop:validity}; it
carries no information about power or about closeness of the network to the Bayes score, and
this proposition makes no use of it. The mathematical content of
Proposition~\ref{prop:avgpower} --- Bayes posterior monotone in the mixture likelihood ratio,
Neyman--Pearson against the mixture, Wald's weighted power --- is classical; the contribution of
this paper is its combination with the composite-null bootstrap calibration and the residual-map
compression, not the proposition in isolation.
\end{remark}
```

==================================================================
PROPOSITION 2 — BODY PARAGRAPH (corrected)
==================================================================

```latex
Proposition~\ref{prop:avgpower} makes precise what the trained statistic is aiming at. The
training library plays the role of a prior $\Pi$ over departure functions, and the balanced-class
Bayes score for separating the library's null class from its $\Pi$-mixture alternative is a
strictly increasing function of the mixture likelihood ratio of the residual map; thresholding it
is therefore the Neyman--Pearson test of the mixture null against the mixture alternative, and it
maximizes $\Pi$-average power among all tests that are measurable with respect to the map. This
is the classical weighted-average-power notion of optimality
\citep{Wald1943,King1987,AndrewsPloberger1994}, imported unchanged: training a classifier on
simulated departures amounts, in the ideal limit, to an implicit choice of the weight measure
over alternatives. The claim is deliberately modest in three respects, recorded in
Remarks~\ref{rem:scope}--\ref{rem:prior}: the optimality is relative to the lossy residual map,
with average (not configuration-wise) size under the training null; the deployed test itself lies
outside the optimality class, since its bootstrap critical value is not map-measurable --- the
proposition certifies the ordering the score aims at, while finite-sample size at the user's null
is the subject of Proposition~\ref{prop:validity} and is guaranteed there under its conditions;
and the deployed network only estimates the Bayes score from a finite corpus, so the oracle power
is an idealization, not a finite-sample guarantee.
```

==================================================================
BIBLIOGRAPHY ADDITIONS (beyond the drafts' verified lists)
==================================================================

- Robins, J.M., van der Vaart, A., Ventura, V. (2000). Asymptotic distribution of p values in composite null models. JASA 95(452), 1143–1156. [Cited in Remark rem:rvv — asymptotic ancestor; must be engaged for JASA.]
- Bayarri, M.J., Berger, J.O. (2000). P values for composite null models. JASA 95(452), 1127–1142. [Companion.]
- Davidson, R., MacKinnon, J.G. (1999). The size distortion of bootstrap tests. Econometric Theory 15(3), 361–376.
- Davison, A.C., Hinkley, D.V. (1997). Bootstrap Methods and Their Application, CUP, §4.5 (double bootstrap). [Operational basis of the calibrated-scale classifier.]
- Grosso et al., Goodness of fit by Neyman–Pearson testing (NPLM line), arXiv:2305.14137. [Prop 2 related work, adjacent not overlapping.]
- Conformal C2ST, arXiv:2507.17026. [Prop 2 Remark estimation — validity decoupled from classifier quality in the two-sample case.]
- King (1987), not 1988 — Econometric Reviews 6(2), 169–218 (already corrected in the draft; keep).

==================================================================
NOTE TO THE EXECUTOR (Opus) — verify by hand before submission
==================================================================

1. Recompute Example 1 end-to-end (G_t = N(−ct,1); marginal mixture equality; π = Φ̄(ε+2cV); limit ½) and confirm the text now calls θ̂ = V a "measurable, non-ML estimator" consistent with the generalized Setting everywhere (search for stray "MLE" outside Remark rem:impossible).
2. Verify the binomial–beta identity h_k'(x) = −B·C(B−1,k−1)x^{k−1}(1−x)^{B−k} and the k = ⌊α(B+1)⌋ edge case for α(B+1) integer, e.g. B = 999, α = .05.
3. Check the failure-convention bookkeeping: on the estimator-failure event p = π = π* = 1, h_k(1) = 0 for k ≤ B, and F_{π*}(u) ≤ u still holds — the proof asserts all three; confirm no other statement quantifies over π without the convention.
4. Confirm the empirical AUC pipeline (both the .62 and the .5011 measurements) scored the two classes on scores from independent simulated datasets, one score per dataset — the Corollary proof now requires this; if any within-dataset couples were used, the numbers must be recomputed before quoting them next to Corollary cor:auc.
5. Verify RVV (2000) and Bayarri–Berger (2000) volume/issue/pages against JASA 95(452), and King 1987 (not 1988) in the .bib; run the citation-audit rules in the memory file (do not "fix" the known false positives).
6. Grep abstract, introduction, and conclusion for "optimal", "valid", "exact": deployed test must never be called optimal (Prop 2 covers the idealized ordering only) and validity must always carry its condition or the quantitative bound.
7. Confirm every \label/\ref pair compiles (rem:rvv, rem:honest1, rem:honest2, ConformalC2ST/NPLM bib keys are new) under xelatex and under the JASA template's theorem environments.
8. Check Remark rem:direction and the body paragraph both say "consistent with" (not "counterpart of") for the 15-rule monotonicity, and that the .5011 sentence no longer claims to estimate the slack.
9. If feasible before submission, actually run the inner-bootstrap classifier on π vs π* at n = 100 and report its estimated calibrated-scale AUC in the empirical section — it is the certificate the theory says to measure, and a referee will ask.
10. Cross-check Prop 2's standing assumption against Section 2's simulator description (finitely parameterized ⇒ Markov kernel automatic) and make sure Section sec:matched exists and matches Remark rem:prior's reference.