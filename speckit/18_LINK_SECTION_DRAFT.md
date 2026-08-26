# 18 — THE LINK-FAMILY SECTION: DRAFTED MANUSCRIPT TEXT (2026-08-21)

Drop-in prose for the manuscript, drafted from `results/probit/VERDICT.csv` (gate PB, 9 blocks,
81,000 datasets). Every number here is reproducible by `python code/probit_verdict.py`.

**The single most important instruction for the executor:** this material is TWO findings of
DIFFERENT KINDS, and the paper must carry both. Reporting only the first would be selective
reporting of exactly the sort that produced the 2026-08-17 retraction.

| | probit | complementary log-log |
|---|---|---|
| what it is | a **shared blind spot** — a statement about the PROBLEM | a **genuine boundary** — a statement about OUR METHOD |
| arms above chance (Bonferroni over the 11, 6000 alts/cell, SE .0028, threshold .0579) | n=50 **none** · n=100 only Tsiatis .0620 + EF .0613 · n=200 **none** | n=50 two · n=100 three · n=200 **six**, led by DEF.poly2 at .2438 |
| honest verb | "no method detects it" | "the directional tests detect it and ours does not" |

---

## §X. Link misspecification: a shared blind spot and a real boundary

*(Discussion / boundary section. ~½ page. No new section heading of its own if space is tight —
fold under the existing limitations discussion.)*

Two results follow, and they differ in kind.

**A shared blind spot.** Against a probit link — the alternative closest to the logistic, differing
from it only in tail weight — no procedure we examined detects misspecification at these sample
sizes. Across n = 50, 100 and 200, with 6000 alternative datasets per cell and every test placed at
exactly the .05 level by a randomized matched-level device, no test's power was distinguishable
from chance under a Bonferroni correction across the eleven, the single exception being a marginal
excess at n = 100 (Tsiatis .062 and EF .061 against a .058 threshold). The pretrained test read
.042, .051 and .046. This is a statement about the problem rather than about any method: at
n <= 200 the logit-probit distinction carries too little information for any of twelve procedures to
recover, and an analyst should not expect a goodness-of-fit test to adjudicate it.

**A real boundary.** Against a complementary log-log link the arms separate, and not in our favour.
The directional tests succeed where the pretrained test does not: at n = 200 a polynomial-in-index
test reaches .244 against the pretrained test's .106, and the gap widens with sample size
(-.047, -.028 and -.138 at n = 50, 100 and 200). The reason is structural. A change of link is
monotone in the linear predictor, which is the axis on which directional tests are built; the
residual map is indexed by covariate ranks and detects departures that occupy a region of covariate
space. A change of link occupies no region — it shifts the whole surface — so the map has no
location to find.

What separates the tests that succeed from those that fail is not which of them was designed for
link misspecification, but which axis each one resolves. The complementary log-log link is
asymmetric and departs from the logistic most in the upper probability tail. Tests that give that
region its own resolution find it: directed polynomial terms in the linear predictor reach .244 and
.144, and Hosmer-Lemeshow with equal-width bins reaches .140. Tests that average the region away do
not: the same Hosmer-Lemeshow statistic with equal-count deciles reaches .047, and a test directed
specifically at Stukel's generalized-link family reaches .068. **The same statistic under two
binning rules differs threefold, which locates the effect in the resolution of the fitted-
probability axis rather than in any test's design intent.**

This is the argument of Figure 3 seen from its other side. The classical battery reads the
fitted-probability axis, either by grouping on the fitted risk or by testing directed departures in
the linear predictor; that is why a departure localized in a region of covariate space is averaged
away by it, and equally why a change of link — which lives entirely on that axis — is not.

We therefore state the method's scope rather than defend it. The pretrained test detects departures
that are localized in covariate space, where the classical battery is weakest (Figure 3), and does
not detect departures that are monotone in the linear predictor, where the directional battery is
strongest. The two are complementary, and an analyst whose concern is link form should use a
test that resolves the fitted-probability axis — the directed variants of the Ebrahim-Farrington
test were the strongest arm in every link cell we examined. Combining the two does not resolve this: the Cauchy combination of the eleven
reaches .104 at n = 200 against its own best component's .244, and adding the pretrained statistic
gives .115. The combination buys robustness by diluting a single informative component.

## §Y. Sensitivity to the training corpus (limitations)

*(This paragraph is MANDATORY. It is the honest counterweight to every learned-power number in the
paper, and it is the measurement that vindicates the earlier retraction.)*

The power of a pretrained test depends on the simulated corpus it was trained on, and that
dependence is not small. Training three rules on three independent corpora and evaluating them on
the same held-out data, power against a complementary log-log alternative at n = 200 read .023,
.125 and .172 — a range of .149. On the identical datasets the best classical test read .249, .254
and .230, moving only by Monte Carlo error, because an analytic test never sees a training corpus.
We therefore report every learned-power figure in this paper as a mean across three corpora with
its range attached, and we recommend the same of any pretrained test. A single-corpus power figure
is not a property of the method.

## §Z. Small-sample fragility of the battery (Table 1 companion)

At n = 50 the Pigeon-Heyse statistic returns a missing value on 2.50% of datasets; separation, by
contrast, is negligible on the designs studied here (0.03% at n = 50, zero beyond). Both the missing
values and separation vanish by n = 200. Because the Cauchy combination
averages a transform of all eleven p-values, one missing component propagates to the combination
unless it is computed with the missing entries skipped. We treat a test that fails as a test that
cannot reject, and report the failure rate alongside the power.

---

## Figure disposition

**Figure 5 (`code/figs/fig5_probit.py`) goes to the SUPPLEMENT, not the body.** PB3 failed on
cloglog, so the across-seed band is wide and the figure would over-promise stability in a body
position. Caption, as a declarative sentence per house style:

> **Figure S-k.** Against a probit link no test of twelve exceeds chance at n <= 200, while against a
> complementary log-log link the directional tests reach .244 at n = 200 and the pretrained test
> .106. Grey dots are the eleven classical tests, one per training corpus; the band on the
> pretrained test is its range across the three corpora. All arms are randomized to exactly
> alpha = .05 on held-out nulls, so vertical position is power alone. Neither link family appears in
> the training corpus.

## Forbidden in this section

Claiming any probit or link-family advantage. Quoting a single-corpus number (PB3 failed). Calling
the probit result a demonstration of the method. Presenting the `probit_extreme` cell as power — it
is a conservatism measurement in a near-deterministic regime. Reporting the blind-spot paragraph
without the boundary paragraph. **Claiming that only tests purpose-built for link families beat the
pretrained test — this is FALSE and a referee will check it in minutes: `HL-equalwidth` (.1398, not
purpose-built for anything) beats it, and `DEF.stukel` (.0683, built precisely for link families)
loses to it.** The correct axis of explanation is which tests resolve the fitted-probability scale,
never which tests were designed for what.
