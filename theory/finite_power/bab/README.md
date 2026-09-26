# Probabilistic branch-and-bound for the finite-sample power certificate: prototype and negative result

Setting: map m = t v + e, e ~ N(0, I_36), frozen network S (`../crown_core.py`, NET), t_bar = 2.4865.
Goal: the smallest t with P{S(t v + e) > t_bar} >= .95 PROVED. Baseline `../07_pcrown.py`: stripe t = 40.06.
Monte Carlo truth: 5.52. The 10x-in-n target is t <= 17.4.

## Headline

| stripe_x1 | certified t | n ratio to truth | regions |
|---|---|---|---|
| 07_pcrown (baseline) | 40.06 | 52.7x | 1 |
| usage-aware charge, z optimised (`07_certify.py`) | 38.34 | 48.3x | 1 |
| + stratified BaB, k = 1 (`03_bab.py`) | **38.00** | 47.4x | 33 (303 s) |
| target | 17.4 | 10x | - |

All five directions (usage-aware, one region): step 30.78 (was 32.07), saddle 44.84 (47.33), diffuse
32.83 (34.37), random 58.99 (62.47). Each certificate is 4-6% smaller in t, 9-12% smaller in n. All were
validated on 200,000 draws (`CERT_UA.log`, `VALID_run_stripe_x1_grad_k1_t38.json`).

**BaB does not reach the target, and the diagnostics show that no practical amount of branching can.**

## Why: the bound falls off a cliff as z grows (01_zsens.py, 06/07 logs)
Every intermediate band is a CROWN linear form plus or minus z Gaussian sd. The relaxation offsets
compound through 10 layers, so the final intercept a(z) falls steeply as z grows. At t = 38 it falls
about 2 units for every +0.05 in z (a = 17.4 at z = 4.0; 6.5 at 4.25; -3.9 at 4.5). The final sd is only
about 2.5. The certified t is simply where a(z*) crosses t_bar, with z* ≈ 4.2-4.5 set by the band budget.
Near the certified point the true S is 60 or more, but the lower function is about 50 units below it
(validation "min gap" 30-53). This slack is compounded relaxation, not noise.

## What was tried, measured
1. **Input-space stratification (the brief's option 1), `03_bab.py`.** Pin k orthonormal directions U
   (top eigenvectors of E[∇S∇S']; 96% of the gradient energy lies in one direction). Cells are
   Gaussian-quantile boxes in s = U'e, so P(C) is exact. The remaining directions are handled
   probabilistically, conditional on the cell. k = 1 certifies 38.0 with 33 cells. k = 2 is worse at
   equal cost (bound .44 after 101 evaluations at t = 38).
   *Limit (`02_pinning_limit.py`)*: point cells are the infinite-refinement limit, and even they FAIL at
   t = 17.4 unless at least 24 of the 36 directions are pinned. The centre cell with k = 24 still has a
   bound of .20; only k = 36, a fully deterministic point, succeeds. Covering 95% of the Gaussian mass with
   near-zero-width boxes in 24 or more dimensions is exponential, beyond any budget.
2. **Pre-activation splitting, BaB-prob style (Wang & Tsiotras 2025, arXiv 2509.25647), `04_relu_bab.py`.**
   Split units are made exact in the final back-substitution. The leaf region is over-approximated on the
   band event E by the halfspaces {U_j(e) >= 0} or {L_j(e) <= 0}, and the leaf probability is bounded by
   the exact bivariate normal (Owen's T), taking the min over ancestors' functions and split constraints.
   **The bound increases with splitting**: at t = 39 the root is .0596; after 300 leaves it is .0695
   (first-layer splits only) and 1.5 (any layer). The relaxed forms of deep units are so loose that the
   two halfspaces of a split overlap heavily, so the children's bounds sum to more than the parent's. A
   sound fix needs tighter intermediate forms, which is the same wall as below.
3. **Sharper accounting of P(E^c).**
   - *Oracle ceiling (`05_accounting_ceiling.py`)*: charging the true P(E^c), estimated from 40k draws,
     instead of any bound still certifies only down to about t = 35. The true P(E^c) hardly depends on t:
     .025 at z = 4.0 against the union bound's .20.
   - *Hunter-Worsley (`hw.py`, `06_pcrown_hw.py`)*: the maximum-spanning-tree bivariate bound cuts the
     union charge by only about 30%.
   - *Usage-aware charge (`bab_core.cell_pcrown_ua`)*: a relaxation is charged only for the statements it
     actually uses. This roughly halves the charge. It is the main source of the 40.06 → 38.34 gain.

## Soundness argument (usage-aware, stratified; k = 0 is the one-cell case)
1. s = U'e ~ N(0, I_k) and e_perp = (I - UU')e are independent. A cell C is a box in s with
   P(C) = Π(ph_i - pl_i), exactly.
2. For any linear form f = λ'e + c: on C, f >= c + min_box(U'λ)'s + G, where G = λ_perp'e_perp ~ N(0, |λ_perp|²)
   is independent of s. So the event {G < -z|λ_perp|} has conditional probability Φ(-z), exactly.
3. CROWN is run layer by layer; each ReLU band is [a_lo - z s_lo, -a_up + z s_up] (the box minimum is in a).
   Every band, the z_C and the set of "used" units are deterministic functions of (t, U, cell, z).
4. The only relaxations used are: stable-active x -> x (needs x >= 0); stable-inactive x -> 0 (needs x <= 0);
   unstable with a triangle upper line and a nonzero negative coefficient anywhere (needs l <= x <= u);
   unstable met only with nonnegative coefficients: lower line 0 or x, valid for EVERY x (needs nothing).
5. Each needed fact follows from one Gaussian statement about that unit's own CROWN form. That form is
   valid whenever the earlier needed facts hold (induction over layers). The statements fail with
   conditional probability Φ(-a_lo/s_lo), Φ(-a_up/s_up), or Φ(-z) per side. Their sum is the charge ch_C,
   by the union bound.
6. On the complement, S(t v + e) >= F_C(e) = a_C + min_box g's + G_fin (CROWN soundness). Hence
   P(S <= t_bar | C) <= ch_C + Φ((t_bar - a_C)/sd_C).
7. Summing over the partition: P(S <= t_bar) <= Σ_C P(C) min(1, ch_C + tail_C). Cells with an infinite
   box end count as bound 1.
8. Certified when the total is <= .05. Nothing is estimated from samples. Sampling is used only for
   validation.
9. Validation (200k draws): the share of draws with any used relaxation invalid was .0049-.0055 against
   charges of .045-.048 (k = 0), and .0025 against .024 (k = 1). There were no lower-function breaches,
   and MC power was 1.000.

## Can BaB reach 10x? Assessment
No, not with CROWN-type linear bounds on this network. The obstacle is not the final Gaussian tail and not
a few influential ReLUs. It is the compounding of about 2,000 undecided units across 10 layers, driven by
bands that must be about 4.2 sd wide to hold jointly with probability .95. Three measured ceilings bound
the remaining room:
- perfect band accounting: about 35;
- perfect intermediate bounds (`../10_oracle_bounds.py`): about 17.6;
- input stratification: needs 24 or more pinned dimensions at 17.4.

Pre-activation splitting made the bound worse. Reaching 17.4 would take something structurally different:
- bounds that do not compound, such as exact per-layer MILP/SDP intermediate bounds, which the oracle
  diagnostic says would at best just touch 17.6;
- or a certificate that is not a lower linear function at all.

Either costs orders of magnitude more than this prototype. The honest number for the paper stays at the
21-57x bracket (now 19-53x with the usage-aware charge, if adopted).

## Files
`bab_core.py` (cell probabilistic CROWN, usage-aware charge), `03_bab.py` (stratified BaB),
`04_relu_bab.py` (pre-activation BaB), `hw.py` + `06_pcrown_hw.py` (Hunter-Worsley), `07_certify.py`
(certify + validate, all directions), `08_validate_strat.py` (validate a stratified run). Diagnostics only,
not certificates: `00_explore.py`, `01_zsens.py`, `02_pinning_limit.py`, `05_accounting_ceiling.py`.
`run_stripe_x1_grad_k1_t35.json` predates the usage-aware charge. Results: `RESULT.json`, `CERT_UA*.json`,
`VALID_*.json`, `*.log`. Run with OMP_NUM_THREADS=4.
