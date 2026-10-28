# 40_analyze.R -- Tier 3b analysis. Reads rows_main/; writes tables/*.csv; prints what REPORT.md quotes
# (log: logs/analyze.log). Nothing outside theory/tier3b is written.
# Conventions (as theory/tier3/40_analyze.R):
#   size-adjusted power at .05 by the randomized threshold rule of theory/shipped_level/rank_shipped_liu.R
#   (matched(), applied per data set as phi), calibrated on the 2,000 correct-model data sets of the same
#   design and n; each test evaluated on the data sets where it returns a p-value (NA left out; share reported).
#   Paired difference DG - arm: mean over the alternatives where both have a p-value of phi_DG - phi_arm,
#   paired SE = sd / sqrt(m). "Ahead" = |difference| > 2 paired SE.
T3B <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3b"
TAB <- file.path(T3B, "tables")
A <- do.call(rbind, lapply(list.files(file.path(T3B, "rows_main"), "\\.csv$", full.names = TRUE), utils::read.csv, stringsAsFactors = FALSE))
stopifnot(!anyDuplicated(A[, c("design", "shape", "n", "ci", "rep")]), !anyDuplicated(A$dseed))
CC <- utils::read.csv(file.path(TAB, "C_CHOSEN.csv"))
cat(sprintf("rows_main: %d rows (nulls %d, alternatives %d)\n", nrow(A), sum(A$shape == "null"), sum(A$shape != "null")))
ARMS <- c(`DeepGOF-1` = "p_dg", `chi-square of map` = "p_ss", `score test (bootstrap)` = "p_score",
          `score test (asymptotic)` = "p_score_asym", GRP = "p_grp", HL = "p_hl", `Stukel (fallback)` = "p_stukel_fb")
phi <- function(p0, p1, a = .05) {                       # rank_shipped_liu.R matched(), per data set
  p0 <- p0[!is.na(p0)]; if (!length(p0)) return(rep(NA_real_, length(p1)))
  u <- sort(unique(p0)); F <- vapply(u, function(c) mean(p0 <= c), 0)
  k <- which(F >= a)[1]; c <- u[k]; below <- if (k > 1) F[k - 1] else 0; g <- (a - below) / mean(p0 == c)
  ifelse(is.na(p1), NA, ifelse(p1 < c, 1, ifelse(p1 == c, g, 0)))
}
pdiff <- function(f1, f2) { ok <- !is.na(f1) & !is.na(f2); d <- f1[ok] - f2[ok]; m <- sum(ok)
  c(d = mean(d), se = stats::sd(d) / sqrt(m), m = m) }

## ---------------- type I error per (design, n) ----------------
TE <- do.call(rbind, lapply(split(A[A$shape == "null", ], list(A$design[A$shape == "null"], A$n[A$shape == "null"]), drop = TRUE), function(z)
  do.call(rbind, lapply(names(ARMS), function(a) { p <- z[[ARMS[[a]]]]; m <- sum(!is.na(p)); r <- mean(p <= .05, na.rm = TRUE)
    data.frame(design = z$design[1], n = z$n[1], test = a, R = m, size = r, se = sqrt(r * (1 - r) / m), missing = mean(is.na(p)),
               outside_2se_of_05 = abs(r - .05) > 2 * sqrt(.05 * .95 / m)) }))))
TE <- TE[order(TE$design, TE$n, match(TE$test, names(ARMS))), ]
utils::write.csv(TE, file.path(TAB, "B_TYPE1.csv"), row.names = FALSE)
cat("\n== type I error at .05 (2,000 correct-model data sets per design x n) ==\n"); print(format(TE, digits = 3), row.names = FALSE)

## ---------------- size-adjusted power per (design, shape, n, C) ----------------
PW <- list(); PD <- list()
cells <- unique(A[A$shape != "null", c("design", "shape", "n", "ci", "C")])
cells <- cells[order(match(cells$design, c("D2", "D10")), cells$shape, cells$n, cells$ci), ]
for (r in seq_len(nrow(cells))) { cc <- cells[r, ]
  X <- A[A$design == cc$design & A$shape == cc$shape & A$n == cc$n & A$ci == cc$ci, ]
  Z <- A[A$design == cc$design & A$shape == "null" & A$n == cc$n, ]
  stopifnot(nrow(X) == 500, nrow(Z) == 2000)
  ph <- lapply(ARMS, function(v) phi(Z[[v]], X[[v]]))
  for (a in names(ARMS)) { f <- ph[[a]]; m <- sum(!is.na(f))
    PW[[length(PW) + 1]] <- data.frame(design = cc$design, shape = cc$shape, n = cc$n, ci = cc$ci, C = cc$C, test = a, R = m,
      power_adj = mean(f, na.rm = TRUE), se = stats::sd(f, na.rm = TRUE) / sqrt(m), power_raw = mean(X[[ARMS[[a]]]] <= .05, na.rm = TRUE),
      missing_alt = mean(is.na(X[[ARMS[[a]]]])), missing_null = mean(is.na(Z[[ARMS[[a]]]]))) }
  for (a in names(ARMS)[-1]) { z <- pdiff(ph[["DeepGOF-1"]], ph[[a]])
    PD[[length(PD) + 1]] <- data.frame(design = cc$design, shape = cc$shape, n = cc$n, ci = cc$ci, C = cc$C, arm = a,
      dg_minus_arm = z[["d"]], paired_se = z[["se"]], m = z[["m"]]) }
}
PW <- do.call(rbind, PW); PD <- do.call(rbind, PD); PD$z <- PD$dg_minus_arm / PD$paired_se
PD$verdict <- ifelse(PD$z > 2, "DeepGOF-1 ahead", ifelse(PD$z < -2, "arm ahead", "within 2 SE"))
utils::write.csv(PW, file.path(TAB, "B_POWER.csv"), row.names = FALSE)
utils::write.csv(PD, file.path(TAB, "B_PAIRED.csv"), row.names = FALSE)
WIDE <- stats::reshape(PW[, c("design", "shape", "n", "ci", "C", "test", "power_adj")], idvar = c("design", "shape", "n", "ci", "C"),
                       timevar = "test", direction = "wide")
names(WIDE) <- sub("^power_adj\\.", "", names(WIDE))
utils::write.csv(WIDE, file.path(TAB, "B_POWER_WIDE.csv"), row.names = FALSE)
cat("\n== size-adjusted power (500 alternatives per cell) ==\n"); print(format(WIDE, digits = 3), row.names = FALSE)

## ---------------- counts ----------------
cnt <- function(z, by) data.frame(by, settings = nrow(z), DG_ahead = sum(z$z > 2), arm_ahead = sum(z$z < -2), within = sum(abs(z$z) <= 2),
  mean_dg_minus_arm = mean(z$dg_minus_arm))
CNT <- rbind(
  do.call(rbind, lapply(split(PD, PD$arm), function(z) cnt(z, data.frame(arm = z$arm[1], design = "all", shape = "all", n = "all")))),
  do.call(rbind, lapply(split(PD, list(PD$arm, PD$design), drop = TRUE), function(z) cnt(z, data.frame(arm = z$arm[1], design = z$design[1], shape = "all", n = "all")))),
  do.call(rbind, lapply(split(PD, list(PD$arm, PD$shape), drop = TRUE), function(z) cnt(z, data.frame(arm = z$arm[1], design = "all", shape = z$shape[1], n = "all")))),
  do.call(rbind, lapply(split(PD, list(PD$arm, PD$design, PD$shape), drop = TRUE), function(z) cnt(z, data.frame(arm = z$arm[1], design = z$design[1], shape = z$shape[1], n = "all")))),
  do.call(rbind, lapply(split(PD, list(PD$arm, PD$n), drop = TRUE), function(z) cnt(z, data.frame(arm = z$arm[1], design = "all", shape = "all", n = as.character(z$n[1]))))))
utils::write.csv(CNT, file.path(TAB, "B_COUNTS.csv"), row.names = FALSE)
cat("\n== settings where DeepGOF-1 / the arm is ahead by > 2 paired SE (of 48) ==\n"); print(format(CNT[CNT$design == "all" & CNT$shape == "all" & CNT$n == "all", ], digits = 3), row.names = FALSE)
cat("\n-- by design x shape, score test (bootstrap) and chi-square of map --\n")
print(format(CNT[CNT$arm %in% c("score test (bootstrap)", "chi-square of map") & CNT$design != "all" & CNT$shape != "all", ], digits = 3), row.names = FALSE)
## mean power by test over the 48 settings
OV <- aggregate(power_adj ~ test, PW, mean); OV <- OV[order(-OV$power_adj), ]
OVS2 <- aggregate(power_adj ~ test + design + shape, PW, mean)
OVS2 <- stats::reshape(OVS2, idvar = c("design", "shape"), timevar = "test", direction = "wide"); names(OVS2) <- sub("^power_adj\\.", "", names(OVS2))
utils::write.csv(OV, file.path(TAB, "B_MEAN_POWER.csv"), row.names = FALSE)
utils::write.csv(OVS2, file.path(TAB, "B_MEAN_POWER_BY_SHAPE.csv"), row.names = FALSE)
cat("\n== mean size-adjusted power over the 48 settings ==\n"); print(format(OV, digits = 3), row.names = FALSE)
cat("\n== mean over n and C, by design x shape ==\n"); print(format(OVS2, digits = 3), row.names = FALSE)
## missingness and diagnostics
cat(sprintf("\nGRP without a p-value: %d of %d data sets (%.2f%%); nulls %.2f%%, alternatives %.2f%%\n", sum(is.na(A$p_grp)), nrow(A),
            100 * mean(is.na(A$p_grp)), 100 * mean(is.na(A$p_grp[A$shape == "null"])), 100 * mean(is.na(A$p_grp[A$shape != "null"]))))
GM <- aggregate(cbind(grp_na = is.na(p_grp), hl_na = is.na(p_hl), stukel_na = is.na(p_stukel_fb), stukel_fallback1 = stukel_fallback == 1,
                      firth = firth == 1, refit_fail_any = n_fail > 0, score_df3 = score_df == 3, axes_x1x2 = axes == "x1+x2") ~ design + n, A, mean)
utils::write.csv(GM, file.path(TAB, "B_DIAGNOSTICS.csv"), row.names = FALSE)
cat("\ndiagnostic shares by design x n:\n"); print(format(GM, digits = 3), row.names = FALSE)
cat(sprintf("score statistic NA in a bootstrap refit: %d; refit failures: %d\n", sum(A$sc_na_boot), sum(A$n_fail)))
cat("\nGRP messages:\n"); print(table(A$grp_msg[is.na(A$p_grp)], useNA = "ifany"))
cat("\ntiming (s per data set):\n"); print(aggregate(cbind(t_combo, t_grp, t_other) ~ design + n, A, mean))
