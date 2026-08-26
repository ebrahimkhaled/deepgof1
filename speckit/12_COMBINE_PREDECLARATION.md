# 12 — COMBINATION TEST: PRE-DECLARATION (frozen 2026-08-20, before the confirmatory run)

## What is claimed and what is not

**Exploratory finding (NOT a claim):** combining the deployed p-values of the residual map and the
classical Cauchy rule, themselves computed from the SAME bootstrap refits, gave power **.5460** at a
common realised size against **.5330** for the classical rule alone — **+.0130**, with SE **.0158**.
That is inside noise. It is the first positive number against the classical benchmark and nothing
more until confirmed.

## The frozen rule

p_combined = Cauchy-combine(p_map, p_cct), equal weights. The Cauchy combination is valid under
ARBITRARY dependence (Liu & Xie, JASA 2020), so no independence assumption is needed and no further
calibration is applied. p_map uses the λ=200 pivoted CNN on the variance-standardised 12×12 map;
p_cct is the classical battery. Both come from the same B=199 bootstrap refits.

## The confirmatory design

n = 400, p = 2 (as before). **4000 fresh alternative datasets** (SE .0079, half the current) and
**2000 fresh nulls**, seed base **71000** — never used. Alternative families drawn from the same
random library but a disjoint seed stream. Nothing is re-tuned: λ, the architecture, the
standardisation and the combination weights are all fixed by this document.

## Success criterion (one number, declared now)

**The combined rule's power at a common realised size exceeds the classical rule's by ≥ +.025**,
i.e. more than 3 SE at the new precision. Size of the combined rule must lie in [.035, .065].

## Fallbacks, declared now

- **Criterion met** → the first construction in this project to beat the classical benchmark under
  deployment. It still requires gates M4/M5 (other designs, other n) before any general claim.
- **Gain positive but < +.025** → report as "no advantage established"; the paper stays the
  price-of-validity result. **Do not re-tune and re-test.**
- **Gain ≤ 0** → the combination adds nothing; record and close the line.

## Why the gain may be capped, stated in advance

Kendall's τ between p_map and p_cct on alternatives is **.617** — the two rules agree substantially,
so they are not detecting largely disjoint departures. A combination gains most when components are
near-independent. This is a reason to expect a modest gain at best, and it is recorded before the
result so it cannot be produced afterwards as an excuse.
