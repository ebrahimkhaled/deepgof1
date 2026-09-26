# Feasibility: an EXACTLY valid DeepGOF-1 by confidence-set maximised Monte Carlo (Dufour 2006)

*One-day study, 2026-09-16. Design D1 (p = 2, independent covariates, 50 % events), n = 100, B = 199,
γ = .01, α = .05. Paired: the deployed plug-in test and the exact variant see the same datasets and
the same common random numbers. 117 null + 117 alternative datasets (the run's 205-minute queuing
deadline stopped the 300 + 300 plan; the balanced paired prefix is used). Worker scripts
`10_core.R` … `15_report.R`; per-dataset rows in `per_dataset.csv`; the numbers below are the
verbatim output of `15_report.R` (`_report_stdout.txt`).*

## The construction
- Confidence set `C_γ` for θ = (b0, β1, β2): Bernstein's inequality coordinate-wise on the score
  S_j(θ) = Σ_i x_ij (y_i − plogis(x_i'θ)) (independent, mean zero, |term| ≤ |x_ij|, variance
  Σ x_ij² π_i(1−π_i) at θ), union bound over the 3 coordinates → coverage ≥ 1 − γ at every n,
  conditional on X. Verified: θ0 ∈ C_γ on 117/117 null datasets (and 400/400 in the bench).
- p_B(θ) = (1 + #{b : T*_b(θ) ≥ T_obs})/(B+1) with y*_b(θ) = 1{U_b ≤ plogis(Xθ)}, one U matrix per
  dataset (common random numbers), each replicate refit + mapped + scored by the shipped pipeline
  (verified to 2.3e-15 against `.deepgof_map()` + `deepgof_score()`).
- p_exact = sup_{θ ∈ C_γ} p_B(θ): 11³ grid over the bounding box of C_γ, then 3 rounds of
  Nelder–Mead from the best 3 grid points; 480 s budget per dataset. Level ≤ α + γ at every n;
  reject at p_exact ≤ .04 for a guaranteed .05.

## Results (verbatim)
```
--- 1. the Bernstein confidence set
  bounding-box half-width / Wald SE, median over 234 datasets: 20.0  35.7  31.3 ; q90 63.1 104.5 102.0
  truncation at 200 Wald SE binds (C_gamma numerically unbounded): 82 / 234 = 0.350
  theta0 inside C_gamma (null arm only): 117 / 117
  11^3 grid points inside C_gamma: median 40, q90 148, max 366 ; zero for 1 dataset
  radial extent of C_gamma in Wald-SE units: min-direction 3.8, median 5.4, max-direction 29.4

--- 2. size (null arm, 117 datasets)
  p_plugin <= .05 (deployed)           2/117 = 0.0171   [0.0021, 0.0604]
  p_exact  <= .05                      0/117 = 0.0000   [0.0000, 0.0310]
  p_exact  <= .04 (= alpha - gamma)    0/117 = 0.0000   [0.0000, 0.0310]
  p_exact >= p_plugin on all null datasets: TRUE ; median (p_exact - p_plugin) 0.105

--- 3. power (alternative arm, 117 datasets)
  p_plugin <= .05 (deployed)          41/117 = 0.3504   [0.2645, 0.4441]
  p_exact  <= .05                     23/117 = 0.1966   [0.1289, 0.2802]
  p_exact  <= .04 (guaranteed .05)    20/117 = 0.1709   [0.1077, 0.2516]
  paired loss plugin.05 - exact.04    +0.1795   [+0.1099, +0.2490]   (discordant 21 vs 0 of 117)

--- 4. where the supremum sat
  sup at theta-hat: null 13.7 %, alt 20.5 % ; when elsewhere, median distance 19.8 (null) / 8.2 (alt) Wald SEs
  Nelder-Mead improved on the grid: 132 / 234 = 0.564 ; median gain 0.025
  search hit the 480 s budget (p_exact then a strict LOWER bound, exactness void): 35 / 234 = 0.150

--- 5. compute
  theta evaluations per dataset: median 115, q90 239, max 442
  refits per dataset: median 17,269, q90 46,722, max 87,958   (deployed test: 199 refits ~ 1.7 s)
  seconds per dataset: median 192, q90 491, max 516, mean 235
  cost ratio exact/plug-in, by refits: median 87 x, max 442 x
```

## Reading (numbers only; the author decides what is published)
1. **Exact validity is attainable.** The construction is sound (θ0 covered 117/117; size 0/117 at
   both thresholds; p_exact ≥ p_plugin on every null dataset).
2. **The price is roughly half the power** at n = 100 on this design: .171 vs .350, paired
   difference +.18 with CI [+.11, +.25] and 21 discordant pairs against 0. This exceeds every gain
   the DeepGOF-2 programme is pursuing.
3. **Compute:** median 87× and up to 442× the deployed test's refits; 3–8 minutes per dataset.
4. **Two structural obstacles:** the Bernstein score set is numerically unbounded on 35 % of
   datasets (Dufour's compactness must be forced by a cap), and the supremum search hit its budget on
   15 %, where exactness is void. A rigorous supremum beyond p = 3 is not practical.
5. The exact variant is therefore a *reference construction*, not a shippable test: it shows what
   exactness costs and why the plug-in bootstrap with the finite-sample bound (Proposition 2) is the
   trade the paper makes.
