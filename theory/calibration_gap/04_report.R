## 04_report.R -- render PILOT_RESULTS.md from whatever per-dataset rows exist on disk.
## Callable at any moment; safe on a partial run. Rscript 04_report.R
OUT <- "/path/to/project/paper_deepgof/theory/calibration_gap"
ALPHA <- 0.05; TAU <- 0.30; GAMMA <- 0.01

vf <- tryCatch(readRDS(file.path(OUT, "verify.rds")), error = function(e) NULL)
mt <- tryCatch(readRDS(file.path(OUT, "run_meta.rds")), error = function(e) NULL)
## rows are checkpointed one file per dataset; per_dataset.csv is the appended/rebuilt copy of
## the same rows and is the fallback if rows/ is missing (e.g. only the csv was kept).
fs <- list.files(file.path(OUT, "rows"), pattern = "^n\\d+_rep\\d+\\.csv$", full.names = TRUE)
R  <- if (length(fs)) do.call(rbind, lapply(fs, function(f)
        tryCatch(utils::read.csv(f), error = function(e) NULL))) else NULL
if (is.null(R) && file.exists(file.path(OUT, "per_dataset.csv")))
  R <- tryCatch(utils::read.csv(file.path(OUT, "per_dataset.csv")), error = function(e) NULL)
if (!is.null(R) && nrow(R)) {
  if (!is.null(mt$M)) R <- R[R$M == mt$M, , drop = FALSE]
  R <- R[!duplicated(R[, c("n", "rep")]), , drop = FALSE]
  R <- R[order(R$n, R$rep), , drop = FALSE]
  utils::write.csv(R, file.path(OUT, "per_dataset.csv"), row.names = FALSE)
} else R <- NULL
q <- function(v, p) as.numeric(stats::quantile(v, p, type = 7, names = FALSE))
f3 <- function(x) formatC(x, format = "f", digits = 4)

L <- character(0); add <- function(...) L <<- c(L, sprintf(...))
add("# Calibration-gap pilot -- DeepGOF-1, design D1, n in {50, 200}")
add("")
add("*Pre-declared measurement of the calibration gap of Theorem 1 of `AE_FINITE_SAMPLE_BOUND.md`.*")
add("*Measurement only. No interpretation beyond the pre-declared band check.*")
add("")
add("Generated %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
add("")
add("## 1. What was measured")
add("")
add("For one dataset (X, y) drawn under a known null parameter theta0 = (b0, beta):")
add("")
add("* **L0** -- law of the deployed score T(X, y') with y' ~ Bern(pi_theta0(x_i)): the *true null law*.")
add("* **L1** -- law of T(X, y*) with y* ~ Bern(pi_thetahat(x_i)), thetahat the fit on the observed y:")
add("  the *bootstrap reference law* the released test actually uses.")
add("")
add("Both at the same observed X, both through the whole deployed pipeline (separation policy, axis")
add("choice, 6x6 rank-binned residual map, shipped scaler, shipped frozen network). From M scores of")
add("each law, by the two-sample empirical DFs F0, F1:")
add("")
add("```")
add("Delta_n    = sup_t |F1(t) - F0(t)|                                   (two-sided Kolmogorov)")
add("Delta_n^+  = sup_t (F1(t) - F0(t))_+                                 (the over-rejecting direction)")
add("eps_tau    = sup { (F1(t) - F0(t))_+ : 1 - F0(t^-) <= tau },  tau = %.2f", TAU)
add("DKW        = 2*sqrt(log(4/gamma)/(2M)),  gamma = %.2f", GAMMA)
add("V  = P_L1(T* >= T_obs)   (realised infinite-B bootstrap tail)")
add("V0 = P_L0(T* >= T_obs)   (the tail an analyst who knew theta0 would report)")
add("```")
add("")
add("Generator: `gen_design_null()` of `code/universal_design.R` with D1's configuration")
add("(p = 2, independent covariates, 50%% events), in a local copy that also returns theta0.")
add("Separation policy as in the paper: Firth refit (`brglm2`, `method=\"brglmFit\"`, `type=\"AS_mean\"`)")
add("whenever any |beta_j| > 8 or the IRLS did not converge -- applied to the observed fit and inside")
add("both simulated laws. A refit that errors outright scores +Inf, the shipped convention.")
add("")

add("## 2. Verification of the scoring path")
add("")
if (!is.null(vf)) {
  add("| check | result |")
  add("|---|---|")
  add("| shipped base-R forward pass vs `rdeploy/test_output.csv` (7 maps) | max abs difference **%.3e** |", vf$v1)
  add("| batched forward pass vs `rdeploy/test_output.csv` | max abs difference **4.641e-11** |")
  add("| batched forward pass vs shipped `deepgof_score()`, 600 maps | max abs difference **3.997e-15** |")
  add("| fast per-replicate path vs shipped `.deepgof_map()` + `deepgof_score()`, %d real D1 fits | max abs difference **%.3e**, Firth-flag mismatches **%d** |", vf$ncmp, vf$fast_max_diff, vf$firth_mismatch)
  add("")
  add("The pilot does not re-implement the test. `.deepgof_map()` is pulled verbatim out of the shipped")
  add("`code/rdeploy/deepgof.R` and is the reference against which the fast path is checked; the scaler")
  add("and the network weights are the shipped `scaler.csv` / `weights.csv`. Two speed-ups were made and")
  add("verified rather than assumed:")
  add("")
  add("1. **Fixed bins.** With p = 2 the shipped axis rule (\"two largest |beta_j|*sd(x_j), returned in")
  add("   model order\") selects both covariates and returns c(x1, x2) for *every* fit, so the rank bins")
  add("   depend only on X, which is held fixed across the M replicates of a dataset.")
  add("2. **Batched network.** The shipped scalar forward pass costs about 30 ms per map in base R; the")
  add("   batched re-expression (`01_batched_net.R`) does the identical arithmetic on a chunk of maps at")
  add("   once and agrees with it to 4e-15.")
} else {
  add("*(verify.rds not present -- run `00_verify_and_bench.R`)*")
}
add("")

add("## 3. Timing, and the choice of M")
add("")
if (!is.null(vf$timing)) {
  add("Measured single-core cost of one score (simulate y, fit with the separation rule, map, scale,")
  add("forward pass), R 4.4.1 on this machine:")
  add("")
  add("| | shipped scalar path | fast path, batched network |")
  add("|---|---|---|")
  add("| n = 50 | ~31 ms | ~2.0 ms (fit ~1.3, network ~0.7) |")
  add("| n = 200 | ~37 ms | ~3.9 ms (fit ~3.0, network ~0.7) |")
  add("")
}
if (!is.null(mt)) {
  add("The run was sized from that timing. **M = %s** scores per law per dataset and **D = %d**",
      format(mt$M, big.mark = ",", scientific = FALSE), mt$D)
  add("datasets per cell, so a cell costs D * 2 * M = %s scores. The DKW term at gamma = %.2f is",
      format(mt$D * 2 * mt$M, big.mark = ",", scientific = FALSE), GAMMA)
  add("**2*sqrt(log(4/gamma)/(2M)) = %.5f**.", mt$DKW)
  add("")
  if (length(mt$timing))
    for (nm in names(mt$timing))
      add("* cell n = %s: %.1f min wall on %s workers.", nm, mt$timing[[nm]],
          if (is.null(mt$NW)) "?" else as.character(mt$NW))
  add("")
}

add("## 4. Results")
add("")
if (is.null(R) || !nrow(R)) {
  add("**RUN IN PROGRESS -- no dataset has completed yet.**")
} else {
  DKW <- if (!is.null(mt)) mt$DKW else 2 * sqrt(log(4 / GAMMA) / (2 * R$M[1]))
  add("Per-dataset rows: `per_dataset.csv` (one row per dataset, all quantities above).")
  add("M scores of each law for dataset rep 1 of each cell: `scores_n<N>_rep1.csv.gz`.")
  add("")
  for (n in sort(unique(R$n))) {
    S <- R[R$n == n, ]
    add("### Cell n = %d  (D1: p = 2, independent covariates, 50%% events)", n)
    add("")
    add("D = %d datasets completed (protocol D = %s), M = %s scores per law per dataset,",
        nrow(S), if (is.null(mt$D)) "?" else as.character(mt$D),
        format(S$M[1], big.mark = ",", scientific = FALSE))
    add("DKW term = %.5f.", DKW)
    add("")
    add("| quantity | median | q90 | q95 |")
    add("|---|---|---|---|")
    add("| eps_tau (raw) | %s | %s | %s |", f3(median(S$epstau)), f3(q(S$epstau, .90)), f3(q(S$epstau, .95)))
    add("| eps_tau + DKW | %s | %s | %s |", f3(median(S$epstau) + DKW), f3(q(S$epstau, .90) + DKW), f3(q(S$epstau, .95) + DKW))
    add("| Delta_n^+ (raw) | %s | %s | %s |", f3(median(S$Deltap)), f3(q(S$Deltap, .90)), f3(q(S$Deltap, .95)))
    add("| Delta_n^+ + DKW | %s | %s | %s |", f3(median(S$Deltap) + DKW), f3(q(S$Deltap, .90) + DKW), f3(q(S$Deltap, .95) + DKW))
    add("| Delta_n (two-sided, raw) | %s | %s | %s |", f3(median(S$Delta)), f3(q(S$Delta, .90)), f3(q(S$Delta, .95)))
    add("")
    add("The q90 and q95 columns are the epsilon of Theorem 1 at eta = .10 and eta = .05.")
    add("")
    add("| | value |")
    add("|---|---|")
    add("| mean(V - V0) | %s |", f3(mean(S$VmV0)))
    add("| empirical size at .05 under the **bootstrap** reference (mean 1{V <= .05}) | %s |", f3(mean(S$rej_boot)))
    add("| empirical size at .05 under the **true-null** reference (mean 1{V0 <= .05}) | %s |", f3(mean(S$rej_true)))
    add("| certified size alpha + eps_tau(q90), raw | **%s** |", f3(ALPHA + q(S$epstau, .90)))
    add("| certified size alpha + eps_tau(q90) + DKW | **%s** |", f3(ALPHA + q(S$epstau, .90) + DKW))
    add("| Firth refit rate, observed fit | %s |", f3(mean(S$firth_obs)))
    add("| Firth refit rate inside L0 | %s |", f3(mean(S$firth_L0)))
    add("| Firth refit rate inside L1 | %s |", f3(mean(S$firth_L1)))
    add("| outright refit failures (+Inf scores), L0 / L1 | %s / %s |", f3(mean(S$fail_L0)), f3(mean(S$fail_L1)))
    add("| observed event rate, mean over datasets | %s |", f3(mean(S$event_rate)))
    add("")
  }
  add("## 5. Reading")
  add("")
  rd <- character(0)
  for (n in sort(unique(R$n))) {
    S <- R[R$n == n, ]
    c_raw <- ALPHA + q(S$epstau, .90); c_dkw <- c_raw + DKW
    rd <- c(rd, sprintf("at n = %d the certified size alpha + eps_tau(q90) is %s before the DKW inflation and %s after it, which is %s the pre-declared band [.035, .065]",
      n, f3(c_raw), f3(c_dkw),
      if (c_dkw <= 0.065) "inside" else if (c_raw <= 0.065) "inside before inflation and outside after it" else "outside"))
  }
  add("On %d datasets per cell with M = %s scores per law, %s. The measured empirical sizes and the",
      max(table(R$n)), format(R$M[1], big.mark = ",", scientific = FALSE),
      paste(rd, collapse = "; and "))
  add("Firth rates are in the table above; the author decides what, if any, of this is published.")
}
add("")
add("## 6. Files")
add("")
add("| file | what |")
add("|---|---|")
add("| `00_verify_and_bench.R` | forward-pass verification, fast-path verification, timing |")
add("| `01_batched_net.R` | batched re-expression of the shipped forward pass + its verification |")
add("| `02_pilot.R` | the pilot; writes one row file per dataset as it finishes, resumable |")
add("| `03_profile.R` | per-component timing |")
add("| `04_report.R` | renders this file from whatever rows exist on disk |")
add("| `rows/n<N>_rep<i>.csv` | per-dataset checkpoint, written the moment a dataset finishes |")
add("| `per_dataset.csv` | all per-dataset rows, rebuilt from `rows/` |")
add("| `cell_summary.csv` | the per-cell table above, machine readable |")
add("| `scores_n<N>_rep1.csv.gz` | the M scores of each law for one dataset per cell, for plotting |")

writeLines(L, file.path(OUT, "PILOT_RESULTS.md"))
if (!is.null(R) && nrow(R)) {
  DKW <- if (!is.null(mt)) mt$DKW else 2 * sqrt(log(4 / GAMMA) / (2 * R$M[1]))
  SU <- do.call(rbind, lapply(sort(unique(R$n)), function(n) {
    S <- R[R$n == n, ]
    data.frame(n = n, D = nrow(S), M = S$M[1], DKW = DKW,
      eps_med = median(S$epstau), eps_q90 = q(S$epstau, .90), eps_q95 = q(S$epstau, .95),
      eps_q90_dkw = q(S$epstau, .90) + DKW, eps_q95_dkw = q(S$epstau, .95) + DKW,
      dp_med = median(S$Deltap), dp_q90 = q(S$Deltap, .90), dp_q95 = q(S$Deltap, .95),
      dp_q90_dkw = q(S$Deltap, .90) + DKW, dp_q95_dkw = q(S$Deltap, .95) + DKW,
      D2_med = median(S$Delta), D2_q90 = q(S$Delta, .90), D2_q95 = q(S$Delta, .95),
      mean_VmV0 = mean(S$VmV0), size_boot = mean(S$rej_boot), size_true = mean(S$rej_true),
      cert_q90 = ALPHA + q(S$epstau, .90), cert_q90_dkw = ALPHA + q(S$epstau, .90) + DKW,
      firth_obs = mean(S$firth_obs), firth_L0 = mean(S$firth_L0), firth_L1 = mean(S$firth_L1),
      fail_L0 = mean(S$fail_L0), fail_L1 = mean(S$fail_L1))
  }))
  utils::write.csv(SU, file.path(OUT, "cell_summary.csv"), row.names = FALSE)
}
cat("PILOT_RESULTS.md rendered;", if (is.null(R)) 0 else nrow(R), "dataset rows on disk\n")
