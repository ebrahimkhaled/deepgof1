# 01 — C1 PRE-DECLARATION (written BEFORE the run; do not edit after)

**Timestamp of declaration:** 2026-08-15, before any confirmatory data was generated.
**This document is frozen.** Results are reported against it verbatim. The LEGofT program
paid a 78% winner's-curse shrinkage and a referee timestamp audit for doing this wrong once;
this is the discipline applied from the start.

## The rule under confirmation (fixed, no further selection)

The **deep ensemble of 20 nnets** (R `nnet`, size = 8, decay = 0.05, maxit = 300,
init seeds 5001..5020), trained on:
- features: 11 battery z-scores + grouped residuals dual-coded (risk-ordered g1..g10,
  sorted s1..s10, hl, mabs) — 33 features
- corpus: 100 random departure families × 100 replicates + 2000 nulls per n
  (`make_family`/`gen_random` as committed in `code/exp_diverse_training.R`, unchanged)

The CNN alternatives were evaluated and LOST (0.7546 and 0.7479 vs 0.7955 at n=1000,
`results/curves/torch_summary_n1000.csv`); C2 is closed. No architecture changes after this
declaration.

## Everything fresh (nothing reused from the exploratory phase)

| Component | Exploratory used | C1 uses |
|---|---|---|
| Corpus family seeds | SEED0 = 1000, 3000, 7000 | **SEED0 = 21000** (never used) |
| Corpus replicate seeds | fseed·10⁵+i | same scheme, new fseeds |
| Null training seeds | 900000+i | **1900000+i** |
| Evaluation cells | stored Stage-1 cells (reused many times) | **freshly generated**: 12 confirmatory families × 1000 reps + conformal null 2000 + test null 2000 per n, seed base "C1" via digest |
| Sample sizes | 200/500/1000 | the same three, declared here |

## Success criteria (all three must hold; declared before the run)

1. **Power:** ensemble mean12 − CCT-full mean12 ≥ **+0.015 at every n** ∈ {200, 500, 1000},
   both computed on the SAME fresh evaluation cells. (+0.015 = roughly half the exploratory
   gains, allowing for selection shrinkage; MC SE of a 12-family mean ≈ .005.)
2. **Size:** ensemble rejection rate on the fresh test-null within **[0.035, 0.065]** at
   every n.
3. **No design leakage:** the gain at each n is positive on the 2-covariate subset alone.

**If any criterion fails, that is the reported result.** No re-runs with new seeds, no
criterion adjustment, no swapping the rule. A failure kills the paper's headline and the
honest fallback is the exploratory-methods paper (failure catalogue + recipe + the need for
larger confirmatory studies).

## What C1 does NOT test (stated now to prevent scope drift later)

- No kind-holdout in C1 (already run, exploratory, reported as such).
- No per-family claims from C1 (single-family trades have corpus-seed variance ±.1;
  only the three means and the three sizes are confirmatory quantities).
- n = 200 remains a small-sample setting where the 10-group features are coarse; C1's
  criterion is the mean, not any single family.

## Execution

Script: `code/c1_confirmatory.R` (written after this declaration, implementing it exactly).
Runtime estimate: ~2.5 h on 12 workers (corpus + fresh evaluation cells are the cost).
Output: `results/C1_confirmatory_n{200,500,1000}.csv` + `results/C1_VERDICT.txt`.
