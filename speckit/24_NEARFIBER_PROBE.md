# 24 — GATE CF: THE NEAR-FIBER PROBE (pre-declaration, frozen 2026-08-22 before any data)

**Frozen.** Authorised by the author 2026-08-22 ("option 1 and 2"). One afternoon, one number.

## 1. The question

Every learned statistic in this project pays a ~.07 "price of validity" because the deployed
reference — the parametric bootstrap at β̂ — is not the true null: a classifier tells them apart at
AUC .62 (phase 12), and any score sensitive enough to be powerful keys on that gap.

But the composite null has an exact escape that no gate has tried: **under H0, conditional on the
sufficient statistic T = X′y, the distribution of y is EXACTLY UNIFORM on the fiber
{y* ∈ {0,1}ⁿ : X′y* = T}** — the exp(β′X′y) factor is constant on the fiber, so β drops out
entirely. A reference drawn from the fiber is exact by algebra, not asymptotics (Cox's conditional
logistic; the exact-test literature). The known obstacle: with continuous covariates the exact
fiber nearly degenerates — hence NEAR-fiber: uniform on {y*: ‖X′y* − T‖ ≤ ε} with standardized
columns, sampled by Metropolis bit-flips/swaps with indicator acceptance.

**If near-fiber replicates are indistinguishable from true nulls (AUC → .50), the .07 tax is
refunded at the root** — sensitive statistics could keep their sensitivity with no pivoting, no
deployment-matching, no price. That reopens the deep route with the nuisance removed, not paid for.

## 2. Design (fixed)

- n = 200, D1 (textbook p=2, the best-understood cell); 500 true-null datasets.
- Per dataset: 199 near-fiber replicates at each tolerance ε ∈ {0.5, 1.0, 2.0} (in standardized-X′y
  units), via swap/flip MCMC, 20n burn-in, n steps between draws; acceptance = indicator of the
  ε-ball. Chain diagnostics (acceptance rate, ‖X′y*−T‖ distribution) deposited.
- **CF-AUC:** train the standard classifier (the phase-12/13 protocol, same CNN on K=6 maps) to
  distinguish near-fiber-replicate maps from FRESH true-null maps. Report per ε.
- **CF-size:** score the maps with a NAIVELY-trained CNN (true nulls vs departures — the exact
  construction that collapsed to size .006 under the parametric bootstrap) and calibrate it against
  the near-fiber reference. Its deployed size is the practical test of the refund.
- Comparators printed alongside: parametric bootstrap AUC (~.62 expected) and its naive-CNN size
  (~.006 expected) on the same datasets — the probe's own baseline, regenerated, not quoted.

## 3. Criteria (frozen now)

| # | Criterion | Reading |
|---|---|---|
| **CF1** | best-ε AUC ≤ **.53** | the route OPENS: nuisance killed by conditioning |
| **CF2** | naive-CNN size under near-fiber calibration ∈ [.035,.065] | the refund is REAL, not just an AUC |
| — | AUC ∈ (.53,.58) | partial: report, no route claim |
| — | AUC ≥ .58 | route stays CLOSED; fiber degeneracy wins; one sentence in the paper |

Both CF1 AND CF2 must hold to claim the route opens. CF1 without CF2 is reported as a curiosity.

## 4. Honest risks, stated in advance

The MCMC may mix badly (tiny ε → frozen chain whose "replicates" are near-copies of y, trivially
passing AUC while being useless — **guarded**: replicates must differ from y in ≥ 10% of positions,
reported); large ε → just an unconditional bootstrap in disguise (the ε-sweep exposes the
tradeoff). Either failure mode is a legitimate CLOSED verdict, reported as such.

## 5. Forbidden

Tuning ε beyond the three declared values. Claiming the route opens on CF1 alone. Re-running with
a different statistic if CF2 fails. Any power claim from this probe — it measures validity
machinery only; power under near-fiber calibration is a FUTURE gate if CF1+CF2 pass.

## 6. Execution

`code/cf_probe.R` (fiber MCMC + maps) → `code/cf_probe_score.py` (AUC + size) → verdict appended
here. Deposit `results/cf/`. Runs AFTER S3 finishes (the machine is S3's until then), BEFORE SEL.
