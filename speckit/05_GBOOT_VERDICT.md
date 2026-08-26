# 05 — G-BOOT VERDICT: the deployed-calibration study (2026-08-19)

**Run:** `code/c_boot_deployed.R`, n = 200, B = 999, 3600 datasets, ~11 h on 12 workers.
**Deposit:** `results/CBOOT_pvalues_n200.csv` — every per-replicate p-value, both rules.
**Raw verdict:** `results/CBOOT_VERDICT_n200.txt`.
Criteria were frozen in `04_MASTER_PLAN.md` §2 before the run. Reported here as they fell.

---


> ## ⚠️ CORRECTION (2026-08-19, same day): the diagnosis below is HALF WRONG
>
> This document concluded that the learned score is not pivotal **while the classical Cauchy
> combination is**. The second half is refuted by `08_UNIVERSAL_VERDICT.md`. That conclusion rested
> on only two designs (2-cov and 4-cov, n=200, 50% events) — too narrow to see the classical rule's
> own miscalibration. Across ten held-out design classes, CCT-full using its analytic null reaches
> **.0810 at p=12 and .0977 in the 5%-event regime**, with spread .0482 versus Route B's .0285.
> **Neither rule is pivotal**, and the learned rule with amortized calibration is the better
> calibrated of the two. Everything else in this document stands.


## The result

**Size is valid. The power advantage is not.**

### Size — the feared failure did NOT happen

| rule | design | realized size | 95% Clopper–Pearson |
|---|---|---|---|
| ensemble | all | **.0417** | [.0311, .0546] ✅ covers .05 |
| CCT-full | all | .0517 | [.0398, .0657] ✅ |
| ensemble | 2-cov | **.0250** | [.0141, **.0409**] ⛔ **excludes .05** |
| ensemble | 4-cov | .0583 | [.0410, .0802] ✅ |
| CCT-full | 2-cov | .0500 | [.0340, .0706] ✅ |
| CCT-full | 4-cov | .0533 | [.0368, .0745] ✅ |

The bootstrap is **not inflated** — the fear that drove this gate was unfounded, and the earlier
.083 alarm (120 reps, B = 99, a reduced 4-feature rule, CI [.034, .132]) was noise, as suspected.

But the split reveals something the pooled number hides. **The ensemble is conservative on the
2-covariate design (.025, half nominal) and near-nominal on the 4-covariate design (.058); the
pooled .0417 is an average of two wrong numbers.** CCT-full is correctly calibrated in *both*
(.050 / .053).

**Uniformity (KS vs U(0,1)):** ensemble **p = .0001 — FAILS** the pre-declared α = .01;
CCT-full p = .236, passes. The ensemble's null p-values mean .537 rather than .5 — shifted
conservative.

### Power — the deployed gain is negative

Both rules got their p-value from the *same* B refits, so this is the fairest comparison in the
project: size-matched by construction, no threshold chosen by anyone.

| family | ensemble | CCT-full | diff |
|---|---|---|---|
| quad | .721 | .584 | **+.137** |
| quad2 | .370 | .360 | +.010 |
| cubic | .020 | .040 | −.020 |
| cloglog | .105 | .080 | +.025 |
| probit | .040 | .060 | −.020 |
| bump | .050 | .155 | −.105 |
| sine | .286 | .578 | **−.292** |
| step | .055 | .125 | −.070 |
| interaction | .455 | .555 | −.100 |
| mixed | .610 | .550 | +.060 |
| multivar_quad | .503 | .271 | **+.232** |
| multivar_int | .330 | .365 | −.035 |

**MEAN12: ensemble .2954 | CCT-full .3103 | DEPLOYED GAIN −.0149.**
Sign test 5 wins / 7 losses, two-sided p = .774 — indistinguishable from a coin flip.

### The comparison that matters

| | idealized null reference | deployed bootstrap |
|---|---|---|
| ensemble mean-12 | .4093 | **.2954** (−.114) |
| CCT-full mean-12 | .2912 | **.3103** (+.019) |
| gain | +.118 (size-matched ≈ +.090) | **−.015** |

**The entire advantage was an artifact of a null reference no analyst can construct.** CCT-full
barely moves between the two regimes; the ensemble loses .114.

## Verdict against the frozen criteria

| # | Criterion | Result |
|---|---|---|
| a | size ∈ [.035,.065], CI covers .05 | **PARTIAL** — holds pooled, **fails on 2-cov** |
| b | null p-values uniform, KS at α = .01 | **FAIL** (p = .0001) |
| c | deployed mean-12 gain positive | **FAIL** (−.0149) |

**G-BOOT FAILS.** The pre-declared fallback applies verbatim: *"Size fine, gain gone → the paper is
a negative result about deployment and still publishable at JCGS/CSDA, but the amortization pitch
is retired."*

## Diagnosis — why, and what it means

The ensemble score is **not pivotal**: its bootstrap null distribution depends on the design in a
way the Cauchy combination's does not. The bootstrap resamples y from the *fitted* model at the
observed X, so it holds β̂ fixed and never reproduces the sampling variability of β̂ itself. A
smooth monotone combination of p-values absorbs that; a 20-network nonlinear score does not — its
observed value sits systematically low in its own bootstrap distribution on 2-covariate designs.

This is precisely the object **Proposition 2** was written to address (`04_MASTER_PLAN.md` §2,
G-THEORY). We now have evidence that its hypothesis **fails for the learned statistic and holds for
the classical one**. That is not a footnote — it is the finding.

## What the paper becomes

**The honest, and genuinely more interesting, thesis:**

> A learned test statistic can be substantially more powerful than any classical combination when
> calibrated against a null reference drawn from the true model — and lose that entire advantage,
> along with exact level, under the only calibration an analyst can actually perform. The gap is
> not a tuning problem; it is a failure of pivotality that the classical combination does not share.

This is a **deployment-gap result for learned inference**, and it generalizes past this test: any
paper reporting power for a learned statistic calibrated against an oracle null is reporting a
number its readers cannot reproduce in practice. Given how fast simulation-based and amortized
testing is growing (§ `02_C3_LITERATURE.md`), that warning has an audience.

**What is now forbidden.** Do not report +.118, +.090, or any idealized-reference gain as the
paper's result. Do not write "validity does not rest on the learner" — measured false at α = .01.
Do not retry calibrations until one passes; the criteria were frozen, and this is the outcome.

## Next steps, in order

1. **Confirm the mechanism.** Is the 2-cov conservatism a pivotality failure or a rank-granularity
   artifact? Cheap check with the deposit already on disk: compare the observed-score percentile
   distribution across designs. **No new simulation needed.**
2. **Does it survive a double bootstrap / β̂-resampling correction?** If a calibration that accounts
   for β̂ variability restores both uniformity and the gain, *that* is a method contribution and the
   paper is positive again. If it does not, the negative result stands and is the paper.
3. **Re-run at n = 500** with `code/fast_battery.R` (4.4× cheaper, verified identical) to check
   whether the gap narrows with sample size — the pivotality story predicts it should.
4. Then re-decide the venue. JCGS remains correct for either outcome.

## Two process notes

- **The D2 fix paid for itself within minutes.** A bug in this script (`mean()` without
  `na.rm = TRUE` for the CCT statistic) printed `NA` for the mean-12 gain on 5 of 3600 replicates.
  Because every per-replicate p-value was deposited, the whole power table was recovered by
  re-scoring the CSV — instead of repeating an 11-hour run. This is exactly the standing rule.
- **The pre-declaration did its job.** The fallback for this outcome was written before the run, so
  there was nothing to renegotiate when the result arrived.
