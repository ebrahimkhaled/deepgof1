# W2 — An optimality anchor for DeepGOF-1

**Status 2026-09-23: statement and proof drafted; measurement in `02_np_gap.py`; literature check and
adversarial pass pending.**

## 1. Why this is needed
Both rejections asked why *this* network should be a good statistic. The paper's answer so far was
empirical. The natural theoretical anchor for a statistic learned by classification is the
Neyman–Pearson lemma. The network was trained to separate null maps from alternative maps drawn from a
stated prior, so the best it can do is the most powerful test of "null" against "prior-averaged
alternative", and the question is how close it gets.

## 2. Setting
The map $m\in\mathbb R^{36}$. $P_0$ is its law under the null. The null class is bootstrap replicates, as
in training. $\Pi$ is the shipped departure prior (`release/deepgof1/training/01_corpus.R`), and
$P_1=\int P_h\,d\Pi(h)$ is the prior-averaged alternative law. Suppose $P_1\ll P_0$ with density ratio
$f=dP_1/dP_0$. A **map-based test** is any measurable $\psi(m)\in[0,1]$. Its level is $E_0\psi$ and its
prior-averaged power is $E_1\psi$.

$\varphi^*=1\{f>k\}$ (randomised on $\{f=k\}$) with $E_0\varphi^*=\alpha$ is the Neyman–Pearson test.
$\varphi_S=1\{S>c\}$ (randomised on ties) with **exactly** $E_0\varphi_S=\alpha$ is the score test at the
same level. Randomising on $\{S=c\}$ always achieves exact level. **Exact equality is needed.** With only
$E_0\varphi_S\le\alpha$ the identity gains the term $e^{c-b}(\alpha-E_0\varphi_S)$. Counterexample without
it: $S=\log f+b$ exactly, thresholded at a smaller level $\alpha'<\alpha$, gives $E_0|f-\hat f|=0$ but a
positive gap. (Adversarial review 2026-09-23.)

**Scope.** $\varphi_S$ is the *oracle* score test, thresholded at the exact null quantile of $S$ under
$P_0$. It is not the deployed bootstrap test at $B=199$; Theorem L and Proposition 1 connect the two.
The upper bound $E_1\varphi^*$ holds for every test whose level is at most $\alpha$ under each null law
in the mixture, because such a test also has level at most $\alpha$ under the mixture.

## 3. The theorem
**Theorem NP.** For every map-based test $\psi$ with $E_0\psi\le\alpha$,
$$E_1\psi\;\le\;E_1\varphi^*.$$
And for every constant $b\in\mathbb R$,
$$0\;\le\;E_1\varphi^*-E_1\varphi_S\;\le\;E_0\bigl|f-e^{S-b}\bigr|.$$
More sharply, the gap equals
$E_0\bigl[(\varphi^*-\varphi_S)(f-\hat f)\bigr]+E_0\bigl[(\varphi^*-\varphi_S)(\hat f-\hat k)\bigr]$
with $\hat f=e^{S-b}$ and $\hat k=e^{c-b}$, and the second term is $\le0$.

*Proof.* The first claim is the Neyman–Pearson lemma. For the second: $E_1\psi=E_0[\psi f]$ for any
test. Since $E_0\varphi^*=E_0\varphi_S=\alpha$, for any constant $\hat k$,
$$E_1\varphi^*-E_1\varphi_S=E_0[(\varphi^*-\varphi_S)f]
=E_0[(\varphi^*-\varphi_S)(f-\hat f)]+E_0[(\varphi^*-\varphi_S)(\hat f-\hat k)].$$
Take $\hat k=e^{c-b}$. Then $\varphi_S=1\{\hat f>\hat k\}$ up to ties, and pointwise
$(\varphi^*-\varphi_S)(\hat f-\hat k)\le0$: where $\varphi_S=1$, $\hat f\ge\hat k$ and
$\varphi^*-\varphi_S\le0$; where $\varphi_S=0$, $\hat f\le\hat k$ and $\varphi^*-\varphi_S\ge0$; and where
$\varphi_S\in(0,1)$ (the tie $S=c$), $\hat f=\hat k$ and the product is 0 whatever the randomisation. So
the second term is $\le0$. The first is at most $E_0|f-\hat f|$ because $|\varphi^*-\varphi_S|\le1$.
Nonnegativity of the gap is the first claim. $\square$

**Remarks.**
1. The bound involves only the null law and the two density ratios, and it holds for **every** $b$,
   so it can be minimised over $b$. A network trained by cross-entropy on balanced classes estimates
   $\log f$ up to the constant $b$ (the log prior odds), which is why this is the natural quantity.
2. The bound is not needed when $f$ is computable: then the gap itself can be computed. §4 does that
   under a Gaussian surrogate for the map law. (Corrected after adversarial review.) With $b$ chosen so
   that $E_0e^{S-b}=1$, $E_0|f-\hat f|=2\,\mathrm{TV}(P_1,Q)\le\sqrt{2\,\mathrm{KL}(P_1\|Q)}$ for
   $Q=\hat f\,dP_0$. That KL is **not** the classifier's excess cross-entropy, which averages Bernoulli KLs
   between posteriors. Converting one into the other costs an unbounded factor $1/((1-\eta)(1-\hat\eta))$.
   No claim linking the bound to the training loss is made.
3. **Honest scope.** Optimality is **prior-averaged** and **map-based**. It says nothing about a
   single named alternative, where a test aimed at that alternative beats any omnibus test. It says
   nothing about tests that use more than the map either; Corollary 2 (the representation limit) is
   that boundary.

## 4. Measurement (`02_np_gap.py` → `NP_GAP.json`, 2026-09-23)
Held-out maps from the shipped prior: 300 per family, for the 100 training families (fseed 60101+f) and
for 100 fresh families; 30,000 bootstrap-null maps, all on fresh dataset seeds.

**Surrogate world** ($P_0=N(\mu_0,\Sigma_0)$, $P_1$ = the equal-weight mixture of $N(m_f,\Sigma_0)$; the
exact NP test is available; 200,000 draws; SE .001):

| families | NP power | DeepGOF-1 power | gap (SE) | share of NP | theorem bound $\min_bE_0\lvert f-e^{S-b}\rvert$ |
|---|---|---|---|---|---|
| training | .273 | .202 | .071 (.001) | **73.9%** | .702 |
| fresh | .289 | .204 | .085 (.001) | **70.5%** | .711 |

**Real held-out maps.** The fair out-of-sample comparison builds both tests from the training families
and scores them on the fresh families' maps. The **closed-form Gaussian-mixture LR test gets .246 (SE .004)
and DeepGOF-1 gets .225 (SE .003)**.

**Reading.**
* The anchor is a number, **within the Gaussian surrogate**: DeepGOF-1 attains 73.9% (training
  families) and 70.5% (fresh families) of the most powerful map-based test for its prior. The two
  differ by far more than their SE (about .001), so there is a small generalisation loss to new
  families; it is not "the same". (Corrected after adversarial pass 2.) On real maps the comparator is
  the surrogate-LR test, not the NP test, so no "% of optimum" is claimed there.
* The generic bound is vacuous here (.70), because the network's score is not calibrated as a log
  likelihood ratio. The measured gap is the usable quantity, which is why Remark 2 matters.
* **A closed-form statistic beats the network out of sample by .021** (about 4 SE). That statistic is the
  log of a Gaussian mixture over the prior's family means, which is exactly what Theorem NP says to
  approximate. It is a direct lead for DeepGOF-2: the analytic prior-LR statistic, calibrated by the same
  bootstrap. The level of that analytic test under the bootstrap has not been checked, and nothing is
  claimed for it yet.

## 5. Prior work and credit (literature check, 2026-09-23)
* **The inequality is not new mathematics.** It follows in one line from the exact identity for the
  excess type-II error of plug-in NP classifiers in Tong (2013), JMLR 14, 3011–3040 (proof of Proposition
  1, p. 3019), and Zhao, Feng, Wang & Tong (2016), JMLR 17 (proof of Proposition 2.4, eq. A.15; its
  level-difference term vanishes at exact level). Present it as the Neyman–Pearson analogue of
  Devroye, Györfi & Lugosi (1996), Corollary 6.1. Those papers bound the error in sup norm with margin
  conditions; the $L^1(P_0)$ form with a free constant $b$ is a trivial variant.
* **What is new is the use:** a measured optimality gap for a frozen, pretrained goodness-of-fit
  statistic against the prior-averaged NP test for its own published training prior. The closest works
  are:
  * Grosso, Letizia, Pierini & Wulzer (2024), SciPost Phys. 16, 123, Appendix A: per named alternative,
    not prior-averaged.
  * Kim, Ramdas, Singh & Wasserman (2021), Ann. Statist. 49(1), 411–434: relative efficiency of
    classifier-accuracy tests.
  * Cranmer, Pavez & Louppe, arXiv:1506.02169.
* Remark 2 (Pinsker): set $b=\log E_0[e^S]$ so that $\hat f$ is a density ratio.
* Related learned tests, correctly attributed: Lopez-Paz & Oquab (2017), ICLR; Liu, Xu, Lu, Zhang,
  Gretton & Sutherland (2020), ICML, PMLR 119, 6316–6326; **Minwoo Kim, Genton, Huser & Castruccio
  (2025), Stat. Comput. 35(1), art. 22, doi:10.1007/s11222-024-10551-0** (not Ilmun Kim).

## 5b. Original list
Plug-in Neyman–Pearson classification (Rigollet & Tong 2011, JMLR; Tong 2013, JMLR; Tong, Feng & Li 2018,
Sci. Adv.); classifier two-sample tests (Lopez-Paz & Oquab 2017, ICLR; Kim, Ramdas, Singh & Wasserman
2021, Ann. Statist.); learned test statistics (Liu et al. 2020, ICML; Kim et al. 2025, Stat. Comput.).
The inequality is elementary and is very likely stated in the plug-in NP literature. The claim would
then be its **use**: a measured optimality gap for a frozen goodness-of-fit statistic against its own
published prior.
