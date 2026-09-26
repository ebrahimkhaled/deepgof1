# Calibration-gap pilot -- DeepGOF-1, design D1, n in {50, 200}

*Pre-declared measurement of the calibration gap of Theorem 1 of `AE_FINITE_SAMPLE_BOUND.md`.*
*Measurement only. No interpretation beyond the pre-declared band check.*

Generated 2026-09-16 08:46:51

## 1. What was measured

For one dataset (X, y) drawn under a known null parameter theta0 = (b0, beta):

* **L0** -- law of the deployed score T(X, y') with y' ~ Bern(pi_theta0(x_i)): the *true null law*.
* **L1** -- law of T(X, y*) with y* ~ Bern(pi_thetahat(x_i)), thetahat the fit on the observed y:
  the *bootstrap reference law* the released test actually uses.

Both at the same observed X, both through the whole deployed pipeline (separation policy, axis
choice, 6x6 rank-binned residual map, shipped scaler, shipped frozen network). From M scores of
each law, by the two-sample empirical DFs F0, F1:

```
Delta_n    = sup_t |F1(t) - F0(t)|                                   (two-sided Kolmogorov)
Delta_n^+  = sup_t (F1(t) - F0(t))_+                                 (the over-rejecting direction)
eps_tau    = sup { (F1(t) - F0(t))_+ : 1 - F0(t^-) <= tau },  tau = 0.30
DKW        = 2*sqrt(log(4/gamma)/(2M)),  gamma = 0.01
V  = P_L1(T* >= T_obs)   (realised infinite-B bootstrap tail)
V0 = P_L0(T* >= T_obs)   (the tail an analyst who knew theta0 would report)
```

Generator: `gen_design_null()` of `code/universal_design.R` with D1's configuration
(p = 2, independent covariates, 50% events), in a local copy that also returns theta0.
Separation policy as in the paper: Firth refit (`brglm2`, `method="brglmFit"`, `type="AS_mean"`)
whenever any |beta_j| > 8 or the IRLS did not converge -- applied to the observed fit and inside
both simulated laws. A refit that errors outright scores +Inf, the shipped convention.

## 2. Verification of the scoring path

| check | result |
|---|---|
| shipped base-R forward pass vs `rdeploy/test_output.csv` (7 maps) | max abs difference **4.641e-11** |
| batched forward pass vs `rdeploy/test_output.csv` | max abs difference **4.641e-11** |
| batched forward pass vs shipped `deepgof_score()`, 600 maps | max abs difference **3.997e-15** |
| fast per-replicate path vs shipped `.deepgof_map()` + `deepgof_score()`, 480 real D1 fits | max abs difference **1.332e-15**, Firth-flag mismatches **0** |

The pilot does not re-implement the test. `.deepgof_map()` is pulled verbatim out of the shipped
`code/rdeploy/deepgof.R` and is the reference against which the fast path is checked; the scaler
and the network weights are the shipped `scaler.csv` / `weights.csv`. Two speed-ups were made and
verified rather than assumed:

1. **Fixed bins.** With p = 2 the shipped axis rule ("two largest |beta_j|*sd(x_j), returned in
   model order") selects both covariates and returns c(x1, x2) for *every* fit, so the rank bins
   depend only on X, which is held fixed across the M replicates of a dataset.
2. **Batched network.** The shipped scalar forward pass costs about 30 ms per map in base R; the
   batched re-expression (`01_batched_net.R`) does the identical arithmetic on a chunk of maps at
   once and agrees with it to 4e-15.

## 3. Timing, and the choice of M

Measured single-core cost of one score (simulate y, fit with the separation rule, map, scale,
forward pass), R 4.4.1 on this machine:

| | shipped scalar path | fast path, batched network |
|---|---|---|
| n = 50 | ~31 ms | ~2.0 ms (fit ~1.3, network ~0.7) |
| n = 200 | ~37 ms | ~3.9 ms (fit ~3.0, network ~0.7) |

The run was sized from that timing. **M = 100,000** scores per law per dataset and **D = 20**
datasets per cell, so a cell costs D * 2 * M = 4,000,000 scores. The DKW term at gamma = 0.01 is
**2*sqrt(log(4/gamma)/(2M)) = 0.01095**.

* cell n = 50: 65.5 min wall on 8 workers.
* cell n = 200: 69.8 min wall on 8 workers.

## 4. Results

Per-dataset rows: `per_dataset.csv` (one row per dataset, all quantities above).
M scores of each law for dataset rep 1 of each cell: `scores_n<N>_rep1.csv.gz`.

### Cell n = 50  (D1: p = 2, independent covariates, 50% events)

D = 20 datasets completed (protocol D = 20), M = 100,000 scores per law per dataset,
DKW term = 0.01095.

| quantity | median | q90 | q95 |
|---|---|---|---|
| eps_tau (raw) | 0.0051 | 0.0297 | 0.0438 |
| eps_tau + DKW | 0.0161 | 0.0406 | 0.0548 |
| Delta_n^+ (raw) | 0.0114 | 0.0395 | 0.0463 |
| Delta_n^+ + DKW | 0.0223 | 0.0505 | 0.0573 |
| Delta_n (two-sided, raw) | 0.0135 | 0.0395 | 0.0463 |

The q90 and q95 columns are the epsilon of Theorem 1 at eta = .10 and eta = .05.

| | value |
|---|---|
| mean(V - V0) | -0.0017 |
| empirical size at .05 under the **bootstrap** reference (mean 1{V <= .05}) | 0.0500 |
| empirical size at .05 under the **true-null** reference (mean 1{V0 <= .05}) | 0.0500 |
| certified size alpha + eps_tau(q90), raw | **0.0797** |
| certified size alpha + eps_tau(q90) + DKW | **0.0906** |
| Firth refit rate, observed fit | 0.0000 |
| Firth refit rate inside L0 | 0.0000 |
| Firth refit rate inside L1 | 0.0001 |
| outright refit failures (+Inf scores), L0 / L1 | 0.0000 / 0.0000 |
| observed event rate, mean over datasets | 0.4990 |

### Cell n = 200  (D1: p = 2, independent covariates, 50% events)

D = 20 datasets completed (protocol D = 20), M = 100,000 scores per law per dataset,
DKW term = 0.01095.

| quantity | median | q90 | q95 |
|---|---|---|---|
| eps_tau (raw) | 0.0013 | 0.0066 | 0.0129 |
| eps_tau + DKW | 0.0123 | 0.0176 | 0.0238 |
| Delta_n^+ (raw) | 0.0033 | 0.0085 | 0.0167 |
| Delta_n^+ + DKW | 0.0143 | 0.0195 | 0.0276 |
| Delta_n (two-sided, raw) | 0.0052 | 0.0122 | 0.0167 |

The q90 and q95 columns are the epsilon of Theorem 1 at eta = .10 and eta = .05.

| | value |
|---|---|
| mean(V - V0) | 0.0006 |
| empirical size at .05 under the **bootstrap** reference (mean 1{V <= .05}) | 0.1000 |
| empirical size at .05 under the **true-null** reference (mean 1{V0 <= .05}) | 0.1000 |
| certified size alpha + eps_tau(q90), raw | **0.0566** |
| certified size alpha + eps_tau(q90) + DKW | **0.0676** |
| Firth refit rate, observed fit | 0.0000 |
| Firth refit rate inside L0 | 0.0000 |
| Firth refit rate inside L1 | 0.0000 |
| outright refit failures (+Inf scores), L0 / L1 | 0.0000 / 0.0000 |
| observed event rate, mean over datasets | 0.4978 |

## 5. Reading

On 20 datasets per cell with M = 100,000 scores per law, at n = 50 the certified size alpha + eps_tau(q90) is 0.0797 before the DKW inflation and 0.0906 after it, which is outside the pre-declared band [.035, .065]; and at n = 200 the certified size alpha + eps_tau(q90) is 0.0566 before the DKW inflation and 0.0676 after it, which is inside before inflation and outside after it the pre-declared band [.035, .065]. The measured empirical sizes and the
Firth rates are in the table above; the author decides what, if any, of this is published.

## 6. Files

| file | what |
|---|---|
| `00_verify_and_bench.R` | forward-pass verification, fast-path verification, timing |
| `01_batched_net.R` | batched re-expression of the shipped forward pass + its verification |
| `02_pilot.R` | the pilot; writes one row file per dataset as it finishes, resumable |
| `03_profile.R` | per-component timing |
| `04_report.R` | renders this file from whatever rows exist on disk |
| `rows/n<N>_rep<i>.csv` | per-dataset checkpoint, written the moment a dataset finishes |
| `per_dataset.csv` | all per-dataset rows, rebuilt from `rows/` |
| `cell_summary.csv` | the per-cell table above, machine readable |
| `scores_n<N>_rep1.csv.gz` | the M scores of each law for one dataset per cell, for plotting |
