# Hostile-referee transcripts (wf_c4eecd88-3e8)

VERDICT: REPAIRABLE. The core mathematics of Proposition 1(i)-(iii) and its proofs is sound — I verified every step and could not construct a counterexample to any numbered claim. All defects are in the surrounding apparatus (one internally inconsistent remark, one false sentence in a proof, one false word in a remark, one overclaiming sentence in the body paragraph, two missing must-cites). Exact repairs listed below.

## (a) Is the statement true as written?

Yes. I checked each piece by direct computation:

- **Reduction (eq. reduction):** N|D ~ Bin(B, π(D)) is correct (D*_b i.i.d. from P_β̂(D) given D; S(D) is D-measurable; the ≥ tie convention matches the closed survival function in π). The k = ⌊α(B+1)⌋ bookkeeping is right, including the integer-α(B+1) edge case. E h_k(U) = k/(B+1) and the monotonicity of h_k are correct.
- **(i), (ii)(a), (ii)(b), (ii)(c):** all verified. The conditional-independence step in (ii)(b) — S(D*_{1:B}) ⊥ D | β̂ because the conditional law is σ(β̂)-measurable — is stated and is exactly what makes the mixture-of-i.i.d. exchangeability argument go through. This is the place where "D*_b i.i.d. only conditionally on β̂ while D is drawn at β₀" usually kills drafts; this draft handles it correctly by making conditional pivotality the hypothesis rather than a consequence.
- **(iii):** the layer-cake identity for h_k, κ ≤ d_K via the triangle split with F_π*(u) ≤ u, and TV ≤ 2(AUC★−½) via the indicator classifier in Lemma auc are all correct. The genuinely nice observation — π* is super-uniform *exactly*, with no assumptions, because conditionally on D the extra draw is distributed exactly G_β̂(D) — is proved correctly.
- **Example 1:** verified numerically in structure: G_t = N(−ct,1), marginals of S(D) and S(D*) both ½N(c,1)+½N(−c,1), π = Φ̄(ε+2cV), rejection → ½. The counterexample does what it claims.
- **Impossibility remark:** the exponential-family completeness, the singleton-fiber argument (y−y′ ∈ ker(Xᵀ)∩{−1,0,1}ⁿ), the Lebesgue-null design argument, and the Neyman-structure conclusion are all correct.

## (b) Proof gaps and internal inconsistencies (the repairs)

**R1 — The Example violates the Setting's own MLE stipulation.** The Setting fixes β̂ = the MLE. In Example 1, θ̂(D) = V is *not* the MLE of that three-point model (the MLE compares |z−cv|, |z+c|, |z−c| and picks θ̂ = 0 with high probability under P₀). So as written, the counterexample lives outside the stated framework. Trivial repair with a real payoff: nothing in any proof of Prop 1 uses MLE-ness — only that β̂ is a measurable function of D. State the Setting for an arbitrary measurable estimator (which is also more faithful to "S contains its own refit"), and invoke the MLE only in the impossibility remark, where completeness/sufficiency actually need it.

**R2 — Remark "direction" contradicts Remark "impossible".** Remark direction attributes the size-.006 regime to condition (ii)(c). But by the draft's own singleton-fiber argument, conditionally on β̂ the observed score is a point mass δ_{s(t)}, so (ii)(c) requires G_t to be supported on [s(t), ∞) a.s. — (ii)(c) is exactly as degenerate in the flagship model as (ii)(b), yet the draft only concedes this for (a)/(b). A referee will put these two remarks side by side. Repair: say the .006 regime instantiates the *marginal* conservative direction (F_π super-uniform, κ = 0), and that (ii)(c) is the clean sufficient mechanism in models where the conditioning is non-degenerate; do not claim (ii)(c) holds in the deployed model.

**R3 — False sentence in the Corollary proof.** "When the two laws coincide, (f(S(D)), f(S(D*₁))) is an exchangeable pair for every f" — false for the *coupled* pair, which is what that notation denotes: equal marginals plus dependence do not give exchangeability (take W uniform on {1,2,3}, W′ = W+1 mod 3: equal marginals, P(W′>W) = 2/3). The statement is true only for *independent* draws W ~ L(S(D)), W′ ~ L(S(D*)), which is how Lemma auc defines AUC. Repair: rewrite the sentence with independent draws, and add one sentence stating that the empirical AUC protocol must use score samples from independent datasets for the two classes (which the deployed pipeline presumably does — say so).

**R4 — "Oracle AUC is degenerately 1" is false, not just loose.** The model is discrete (y* ∈ {0,1}ⁿ), so P{S(D*) = S(D) | β̂} > 0 and the pair-scale oracle AUC is 1 − O(atom mass), not 1. Also, the remark never says whether the bootstrap-class pair is (S(D*), β̂(D)) or (S(D*), β̂(D*)). Repair: "approaches 1 up to the bootstrap atom at the observed configuration" and fix the pair definition.

**R5 — π and π* are undefined on the MLE-failure event in part (iii).** The p = 1 / S* = +∞ conventions cover the deployed p-value, but κ and the d_K chain quantify over L(π), which is undefined when β̂(D) does not exist. One-line repair: set π := 1 and π* := 1 on the failure event; then h_k(1) = 0 = 1{p ≤ α} for all α < 1, F_π*(u) ≤ u is preserved, and the whole of (iii) goes through. Also note Lemma pit is being applied on the extended real line (atom at +∞ from the bootstrap-failure convention) — it holds there, but say so.

## (c) Scoop check

Not scooped on the two load-bearing claims. My searches found no prior statement of (i) the finite-B reduction κ ≤ d_K(L(π), L(π*)) with π* exactly super-uniform by construction, or (ii) the two-sided AUC delineation with the explicit marginal-equality counterexample. The high-dimensional bootstrap literature proves Kolmogorov bounds between *statistic* laws, not between calibrated p-value laws with a by-construction super-uniform benchmark. But two must-engage omissions, one serious for this venue:

- **Robins, van der Vaart & Ventura (2000), "Asymptotic Distribution of P Values in Composite Null Models," JASA 95:1143–1156.** This is *the* JASA paper on exactly this phenomenon: plug-in (parametric-bootstrap) p-values in composite nulls are non-uniform — conservative when the statistic's asymptotic mean depends on the nuisance parameter — with Bayarri & Berger (2000, JASA) as the companion. It is the asymptotic ancestor of Prop 1 and of the "deployment-matched training ≈ making the statistic asymptotically pivotal" story. Submitting to JASA without citing it invites a desk-level "have they read our own journal" note. It does not scoop the finite-sample reduction, but the paper must position against it explicitly: RVV characterize the asymptotic law of π; Prop 1(iii) bounds the finite-B size distortion by a distance that is estimable by the inner bootstrap.
- **Davidson & MacKinnon, "The Size Distortion of Bootstrap Tests" (Econometric Theory 1999)** — owns the phrase and the asymptotic-order analysis of bootstrap size distortion; cite in the related-work sentence. Also cite **Davison & Hinkley (1997, §4.5)** for the inner/double bootstrap that estimates L(π*), since Remark whichclassifier relies on it operationally.

The draft's own flags (Dufour 2006 owning invalidity-in-general; Barber–Janson 2022 as nearest neighbor; Beran 1988 for prepivoting — note π* *is* Beran's prepivoted root, and the draft says so) are accurate and adequately differentiated.

## (d) Body-paragraph overclaim

One sentence overclaims: "the .5011 figure is a necessary check whose empirical value estimates … that the validity slack is small." False by the paper's own Example 1: the score-scale AUC can be exactly ½ while the slack is ~½; a score-scale AUC therefore does not estimate the calibrated-scale slack at all — it estimates the marginal score-scale TV, which the theorem never connects to size. Repair: "the .5011 figure passes the necessary marginal check; the direct evidence that the slack is small is the deployed size of .052 (and, on the theorem's own scale, an inner-bootstrap classifier on π vs π*)." Similarly weaken "the monotone AUC–size relation … is the empirical counterpart of this chain" (the chain is calibrated-scale, the measurement score-scale) to "consistent with" — the conditions_honest field already has the correct phrasing; the body paragraph regressed from it. The same regression appears in Remark direction ("empirical counterpart of the chain in part (iii)") — align both with conditions_honest.

## Summary of required repairs

1. Restate the Setting for an arbitrary measurable estimator; reserve MLE for the impossibility remark (fixes the Example/Setting inconsistency).
2. Rewrite Remark direction: (ii)(c) is degenerate in the flagship model too; the .006 regime is the marginal conservative direction.
3. Fix the exchangeable-pair sentence in the Corollary proof (independent draws), and state the independent-datasets requirement for the empirical AUC.
4. "Degenerately 1" → "1 up to the bootstrap atom"; specify the pair.
5. Define π = π* = 1 on the MLE-failure event; note Lemma pit on the extended line.
6. Add Robins–van der Vaart–Ventura (2000) + Bayarri & Berger (2000) with an explicit positioning paragraph; add Davidson & MacKinnon (1999) and Davison & Hinkley (1997).
7. Fix the one overclaiming body-paragraph sentence and the "empirical counterpart" phrasing (two places).

With these, the proposition, corollary, example, and impossibility remark are correct, tight, and honestly delimited; the finite-B κ ≤ d_K reduction and the AUC delineation stand as the novel content.

Sources: [Robins, van der Vaart & Ventura (2000), JASA](https://www.tandfonline.com/doi/abs/10.1080/01621459.2000.10474310) · [RVV full text](https://statweb.rutgers.edu/ztan/material/robins-vandervaart-ventura.pdf) · [Davidson & MacKinnon, The Size Distortion of Bootstrap Tests](https://www.academia.edu/32046851/The_Size_Distortion_of_Bootstrap_Tests) · [Beran, impact of the bootstrap (prepivoting overview)](https://www.stat.ucdavis.edu/~beran/boot.pdf) · [Rank-transformed subsampling / prepivoting discussion](https://arxiv.org/pdf/2301.02739) · [Bootstrapping ℓp-statistics in high dimensions (Kolmogorov-distance bootstrap bounds, statistic scale)](https://arxiv.org/pdf/2006.13099)

==== ATTACK ====

# AE Report — Proposition 2 (Π-average-power optimality)

## (a) Is the statement true as written?

The core of (ii) is true: Tonelli collapses Π-average power into power against the mixture Q̄, and the randomized NP lemma does the rest. I could not break (ii). But **part (i) is false in the letter**, by a trivial counterexample:

**Counterexample to the uniqueness clause of (i).** Take T = R, μ = Lebesgue on [0,2], with p₀ and q̄ supported on [0,1]. Any measurable f that equals η* on [0,1] but is arbitrary on (1,2] attains the same risk R(f) = R(η*). So η* is **not** "the μ-a.e. unique minimizer" — uniqueness holds only μ-a.e. **on {p₀+q̄ > 0}**, equivalently (P₀+Q̄)-a.e. The proof already says this correctly ("agrees with η* μ-a.e. on {p₀+q̄>0}"); the statement over-promises relative to its own proof. Fix in one move: take μ = P₀+Q̄ as the dominating measure — then the density assumption is automatic (delete it as an "assumption"; state it as a choice) and {p₀+q̄>0} is μ-full.

Two more letter-level defects:

1. **Q̄ is not well-defined as written.** "Π is a probability measure on departure functions" presupposes a σ-algebra on the space of g's and requires g ↦ Q_g(A) to be measurable for every measurable A ⊆ T (a Markov kernel). Without this, both the definition Q̄ = ∫Q_g dΠ(g) and the Tonelli step in the proof are vacuous. This must be an explicit standing assumption, not folklore. (In this paper it is harmless — the library is a finitely-parameterized simulator, so the kernel property is free — say so in one sentence.)

2. **Boundary of the cross-entropy risk.** R(f) can be +∞ (f ≡ 0 on a Q̄-positive set). The pointwise-minimization argument is valid because the integrand is nonnegative, but the proof never says so; add the clause "the integrand is nonnegative, so R is well-defined in [0,∞] and may be minimized pointwise." Also, at points with q̄(t) = 0 or p₀(t) = 0 the per-t minimizer is at the boundary of [0,1], not in the open interval where the proof asserts strict convexity; one sentence covers it.

## (b) Proof gaps — conditioning, exchangeability, estimated parameters

Here is the good news and the trap. **Proposition 2 genuinely does not need any of the bootstrap/conditioning machinery**: it lives entirely under the training generator, D* never appears, so the "D*_b i.i.d. only conditionally on β̂ while D is drawn at true β" subtlety is Proposition 1's problem, not this one's. I tried to make Prop 2 secretly depend on it and failed — the statement is honestly self-contained. The remaining gaps are the three in (a) plus:

- **Remark 1's bridge sentence is the one soft spot.** "Replacing the P₀-quantile c_α by the data-dependent bootstrap critical value changes the size guarantee but not the ordering" is true, but it hides the sharper point a referee will make: the deployed test 1{S(D) exceeds the bootstrap quantile of S(D*)} is **not map-measurable** — its critical value depends on D through β̂ and the resampling, not through M(D) — so the deployed test is not even in the class over which Proposition 2 optimizes, and no formal optimality transfers to it. The remark half-says this; say it in full, in exactly that sentence ("the deployed test lies outside the class of Proposition 2; what the proposition certifies is that the ordering of T induced by η* is the Π-average-power-optimal ordering, and the deployed procedure composes that ordering with a different, conditionally calibrated cutoff").
- **Remark 3 is correct but should name its cost:** with P₀ redefined as the bootstrap-replicate law, "size" in the proposition is size under a law that is not any null hypothesis a user cares about. One clause suffices.
- NP citation (Lehmann–Romano Thm 3.2.1), the c = k/(1+k) change of scale, {η*>c} = {Λ>k}, and the continuous-Λ non-randomized variant are all correct as written.

## (c) Scooping

The mathematical content — Bayes posterior is a monotone transform of the mixture LR; NP against the mixture = weighted-average-power optimal — is **folklore**, present in Wald (1943), essentially a textbook exercise in Lehmann–Romano, and instantiated for classifier statistics in [Lopez-Paz & Oquab (2017)](https://arxiv.org/pdf/1610.06545), Kim et al. (2021), and Cranmer et al. (2015). The draft already concedes this in its literature field, which is the right posture. My searches found no paper proving this proposition in the composite-null, bootstrap-calibrated, GOF-for-logistic setting — the claimed novelty perimeter (residual-map compression + Prop 1 validity condition + matched training) survives. Two adjacent items the related-work paragraph should add so a referee cannot "discover" them: [Goodness of fit by Neyman–Pearson testing (arXiv:2305.14137)](https://arxiv.org/pdf/2305.14137) — the physics NPLM line, learned NP statistics used explicitly for GOF — and [Conformal C2ST (arXiv:2507.17026)](https://arxiv.org/html/2507.17026), which gives exact finite-sample p-values from arbitrary trained classifiers with power degrading in classifier error; the latter directly reinforces Remark 2's rank-disagreement sentence and predates you on "validity decoupled from classifier quality" in the two-sample case. Neither scoops Prop 2; both must be cited. King 1987 (not 1988) is confirmed correct in the draft.

## (d) Does the body paragraph overclaim?

Mostly disciplined. Two sentences need surgery:

1. "training a classifier on simulated departures is **not a heuristic but** an implicit choice of the weight measure" — the theorem covers the *idealized* rule only; the trained CNN's implicit weight measure is distorted by architecture, corpus, and optimization. Weaken to "training a classifier … amounts, in the ideal limit, to an implicit choice of the weight measure."
2. "finite-sample size at the user's null **is delivered by** the bootstrap calibration of Proposition 1" — Prop 1 delivers size only *under its condition* (the one the AUC merely surrogates). Write "is the subject of Proposition 1" or "is delivered under the condition of Proposition 1." Otherwise Prop 2's body paragraph launders Prop 1's caveat.

The `conditions_honest` field is accurate, including the correct statement that AUC ≈ .5 certifies nothing for Prop 2.

## VERDICT: REPAIRABLE

All repairs are local; none touches the architecture of the result.

1. Set μ = P₀+Q̄ (drop the dominating-measure assumption as an assumption); restate (i)'s uniqueness as (P₀+Q̄)-a.e. — this kills the counterexample in (a).
2. Add the standing assumption that g ↦ Q_g is a Markov kernel on a measurable space of departures (one sentence; note it is automatic for the finitely-parameterized library).
3. In the proof of (i): note the integrand is nonnegative (so R ∈ [0,∞] and pointwise minimization is legitimate) and handle the boundary cases q̄(t)=0 or p₀(t)=0 in one clause.
4. Rewrite Remark 1's bridge: state explicitly that the deployed test is not map-measurable and lies outside Prop 2's optimality class; the proposition certifies the optimal *ordering*, nothing about the deployed *test*.
5. Add one clause to Remark 3 naming what "size under bootstrap-replicate P₀" is not.
6. Body paragraph: the two wording fixes in (d).
7. Related work: add arXiv:2305.14137 (NPLM/learned NP GOF) and arXiv:2507.17026 (conformal C2ST) as adjacent-not-overlapping.

Sources: [Lopez-Paz & Oquab, Revisiting Classifier Two-Sample Tests](https://arxiv.org/pdf/1610.06545) · [Goodness of fit by Neyman–Pearson testing](https://arxiv.org/pdf/2305.14137) · [Conformal C2ST](https://arxiv.org/html/2507.17026) · [Neyman–Pearson lemma (mixture-alternative form)](https://en.wikipedia.org/wiki/Neyman%E2%80%93Pearson_lemma) · [Focused test statistics / weighted-average alternatives](https://arxiv.org/html/2507.17831v2)