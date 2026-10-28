# 40_analyze.R -- Tier 3 analysis. Reads rows_AB/ (items A and B) and rows_C/ (item C); writes CSVs in tables/
# and prints the numbers REPORT.md quotes (log: logs/analyze.log). Nothing outside theory/tier3 is written.
#
# Conventions (those of the published tables):
#   Liu  size-adjusted power at .05 by a randomized threshold on the held-out nulls 1001-2000 of the same
#        setting and n (theory/shipped_level/rank_shipped_liu.R, matched()); a test is evaluated on the data
#        sets where it returns a p-value (NA left out); power pooled over the 4 departure strengths
#        (2,000 alternatives per setting x n), as RANK_SHIPPED_LIU.csv. size = raw rejection at .05 on those nulls.
#   niche randomized threshold on the 2,000 nulls at the same n (theory/niche_paper/analyze_niche_paper.py,
#        matched_phi); 1,000 alternatives per cell. Stored arms: NA -> p = 1 (none occur). New arms: GRP NA left
#        out (author's rule) and its share reported; the other new arms have no NA (checked below).
#   Paired difference DG - arm: mean over the alternatives where both have a p-value of phi_DG - phi_arm,
#        paired SE = sd / sqrt(m); phi = the per-data-set randomized rejection (0, gamma or 1).
T3  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3"
PAP <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
SM  <- Sys.getenv("T3_SMOKE") == "1"                 # smoke: run the code paths on rows_*_smoke
TAB <- file.path(T3, if (SM) "tables_smoke" else "tables"); dir.create(TAB, showWarnings = FALSE)
NETS <- c(paste0("Control_s", 1:3), paste0("NoInter_s", 1:3), paste0("NoQuad_s", 1:3), paste0("NoQuadStrict_s", 1:3))
rd <- function(d) do.call(rbind, lapply(list.files(file.path(T3, d), "\\.csv$", full.names = TRUE), utils::read.csv, stringsAsFactors = FALSE))
A <- rd(if (SM) "rows_AB_smoke" else "rows_AB"); A$C <- round(A$C, 2)
cat(sprintf("rows_AB: %d rows (liu %d, niche %d) | p_dg == stored p-value on every row: %s\n", nrow(A),
            sum(A$design == "liu"), sum(A$design == "niche"), all(A$p_dg == A$p_dg_stored)))
stopifnot(!anyDuplicated(A[, c("design", "S", "n", "kind", "ci", "rep")]))

## ---------- the two randomized-threshold rules ----------
phi_liu <- function(p0, p1, a = .05) {                 # rank_shipped_liu.R matched(), per data set
  p0 <- p0[!is.na(p0)]; if (!length(p0)) return(rep(NA_real_, length(p1)))
  u <- sort(unique(p0)); F <- vapply(u, function(c) mean(p0 <= c), 0)
  k <- which(F >= a)[1]; c <- u[k]; below <- if (k > 1) F[k - 1] else 0; g <- (a - below) / mean(p0 == c)
  ifelse(is.na(p1), NA, ifelse(p1 < c, 1, ifelse(p1 == c, g, 0)))
}
phi_niche <- function(pn, pa, a = .05) {               # analyze_niche_paper.py matched_phi(), per data set
  pn <- pn[!is.na(pn)]; t <- sort(unique(pn)); F <- findInterval(t, sort(pn)) / length(pn)
  out <- rep(NA_real_, length(pa)); ok_a <- !is.na(pa); x <- pa[ok_a]
  ok <- which(F <= a)
  v <- if (!length(ok)) { g <- if (F[1] > 0) a / F[1] else 0; g * (x <= t[1]) } else {
    i <- max(ok); t0 <- t[i]; F0 <- F[i]
    if (i + 1 > length(t)) as.numeric(x <= t0) else { t1 <- t[i + 1]; F1 <- F[i + 1]
      g <- if (F1 <= F0) 0 else (a - F0) / (F1 - F0); (x <= t0) + g * ((x > t0) & (x <= t1)) } }
  out[ok_a] <- v; out
}
pdiff <- function(f1, f2) { ok <- !is.na(f1) & !is.na(f2); d <- f1[ok] - f2[ok]; m <- sum(ok)
  c(d = mean(d), se = stats::sd(d) / sqrt(m), m = m) }

## =====================================================================================================
## 1. Liu Settings 2-4 with the distributed network
## =====================================================================================================
L <- A[A$design == "liu", ]
CD <- if (SM) "rows_C_smoke" else "rows_C"
C_rows <- if (dir.exists(file.path(T3, CD))) rd(CD) else NULL
if (!is.null(C_rows)) { C_rows$C <- round(C_rows$C, 2)
  ## Stukel with fallback on the SAME data sets (block 101000): attach by key, check events agree
  cb <- C_rows[C_rows$seed0 == 101000 & C_rows$S %in% 2:4, ]
  m <- match(paste(L$S, L$n, L$kind, L$C, L$rep), paste(cb$S, cb$n, cb$kind, cb$C, cb$rep))
  L$p_stukel_fb <- cb$p_stukel_fb[m]; L$ev_C <- cb$events[m]
  cat(sprintf("Liu A rows matched to C rows: %d of %d | events agree on all matched: %s\n", sum(!is.na(m)), nrow(L),
              all(L$events[!is.na(m)] == L$ev_C[!is.na(m)])))
} else L$p_stukel_fb <- NA
ARMS1 <- c(`DeepGOF-1 shipped` = "p_dg", GRP = "p_grp", `chi-square of map` = "p_ss", `score test (bootstrap)` = "p_score",
           `score test (asymptotic)` = "p_score_asym", `Stukel with fallback` = "p_stukel_fb")
res <- list(); phis <- list()
for (S in 2:4) for (n in c(50, 100, 500)) {
  b <- L[L$S == S & L$n == n, ]; nul <- b$kind == 0; alt <- b$kind == 2
  ph <- list()
  for (a in names(ARMS1)) { p <- b[[ARMS1[[a]]]]
    ph[[a]] <- phi_liu(p[nul], p[alt])
    res[[length(res) + 1]] <- data.frame(setting = S, n = n, test = a, power = mean(ph[[a]], na.rm = TRUE),
      size = mean(p[nul] <= .05, na.rm = TRUE), missing = mean(is.na(p)), missing_null = mean(is.na(p[nul])),
      missing_alt = mean(is.na(p[alt]))) }
  phis[[sprintf("%d_%d", S, n)]] <- ph
}
NEW1 <- do.call(rbind, res)
OLD1 <- utils::read.csv(file.path(PAP, "theory/shipped_level/RANK_SHIPPED_LIU.csv"), check.names = FALSE)
chk <- merge(OLD1[OLD1$test == "DeepGOF-1 shipped", ], NEW1[NEW1$test == "DeepGOF-1 shipped", ], by = c("setting", "n"))
cat(sprintf("DeepGOF-1 shipped power recomputed vs RANK_SHIPPED_LIU.csv: max |diff| %.2e (size %.2e)\n",
            max(abs(chk$power.x - chk$power.y)), max(abs(chk$size.x - chk$size.y))))
T1 <- rbind(OLD1[OLD1$test != "DeepGOF-1 shipped", c("setting", "n", "test", "power", "size", "missing")],
            NEW1[, c("setting", "n", "test", "power", "size", "missing")])
T1$new_arm <- T1$test %in% names(ARMS1)[-1]
T1 <- T1[order(T1$setting, T1$n, -T1$power), ]
utils::write.csv(T1, file.path(TAB, "T1_LIU_RANK.csv"), row.names = FALSE)
utils::write.csv(NEW1, file.path(TAB, "T1_LIU_NEW_ARMS_DETAIL.csv"), row.names = FALSE)
O1 <- aggregate(cbind(power, size) ~ test, T1, mean); O1 <- O1[order(-O1$power), ]
BYN <- stats::reshape(aggregate(power ~ test + n, T1, mean), idvar = "test", timevar = "n", direction = "wide")
BYN <- merge(BYN, O1, by = "test"); BYN <- BYN[order(-BYN$power), ]
names(BYN) <- c("test", "power_n50", "power_n100", "power_n500", "power_all9", "size_all9")
utils::write.csv(BYN, file.path(TAB, "T1_LIU_OVERALL_BY_N.csv"), row.names = FALSE)
cat("\n== 1. Liu S2-4, distributed network: average size-adjusted power by n (9 blocks) ==\n"); print(format(BYN, digits = 3), row.names = FALSE)
## paired differences DG - new arm
pd <- list()
for (k in names(phis)) { ph <- phis[[k]]; S <- as.integer(substr(k, 1, 1)); n <- as.integer(sub("^._", "", k))
  for (a in names(ARMS1)[-1]) { z <- pdiff(ph[["DeepGOF-1 shipped"]], ph[[a]])
    pd[[length(pd) + 1]] <- data.frame(setting = S, n = n, arm = a, d = z[["d"]], se = z[["se"]], m = z[["m"]]) } }
PD1 <- do.call(rbind, pd); PD1$z <- PD1$d / PD1$se
agg <- function(g) do.call(rbind, lapply(split(PD1, g), function(z) data.frame(n = if (length(unique(z$n)) == 1) z$n[1] else NA,
  arm = z$arm[1], cells = nrow(z), d = mean(z$d), se = sqrt(sum(z$se^2)) / nrow(z))))
PDN <- rbind(agg(list(PD1$arm, PD1$n)), agg(list(PD1$arm))); PDN$z <- PDN$d / PDN$se
PDN$n <- ifelse(is.na(PDN$n), "all", PDN$n)
utils::write.csv(PD1, file.path(TAB, "T1_LIU_PAIRED_CELLS.csv"), row.names = FALSE)
utils::write.csv(PDN, file.path(TAB, "T1_LIU_PAIRED_BY_N.csv"), row.names = FALSE)
cat("\npaired DG - arm (averaged over settings; SE combines the cell SEs):\n"); print(format(PDN, digits = 3), row.names = FALSE)
cat("\nGRP: share of data sets without a p-value, by setting and n:\n")
print(format(stats::reshape(NEW1[NEW1$test == "GRP", c("setting", "n", "missing")], idvar = "setting", timevar = "n", direction = "wide"), digits = 3), row.names = FALSE)

## =====================================================================================================
## 2. Ten-covariate design
## =====================================================================================================
N <- A[A$design == "niche", ]
ARMS2 <- c(DG = "p_dg", GRP = "p_grp", `chi-square of map` = "p_ss", `score test (bootstrap)` = "p_score",
           `score test (asymptotic)` = "p_score_asym")
cat("\n== 2. ten-covariate design ==\nNA counts:", paste(names(ARMS2), sapply(ARMS2, function(v) sum(is.na(N[[v]]))), collapse = " | "), "\n")
SP <- utils::read.csv(file.path(PAP, "theory/niche_paper/niche_power.csv"))
SZ <- utils::read.csv(file.path(PAP, "theory/niche_paper/niche_size.csv"))
cells <- unique(N[N$kind == 1, c("scen", "n", "C")]); cells$fam <- sub("_.*", "", cells$scen)
cells <- cells[order(cells$fam, cells$n, cells$C), ]
pw <- list(); pd2 <- list(); phiN <- list()
for (r in seq_len(nrow(cells))) { cc <- cells[r, ]
  X <- N[N$scen == cc$scen, ]; Z <- N[N$kind == 0 & N$n == cc$n, ]; if (!SM) stopifnot(nrow(X) == 1000, nrow(Z) == 2000)
  ph <- lapply(ARMS2, function(v) phi_niche(Z[[v]], X[[v]])); names(ph) <- names(ARMS2)
  phiN[[cc$scen]] <- ph
  for (a in names(ARMS2)) pw[[length(pw) + 1]] <- data.frame(fam = cc$fam, C = cc$C, n = cc$n, arm = a,
      R = sum(!is.na(ph[[a]])), matched = mean(ph[[a]], na.rm = TRUE), se = stats::sd(ph[[a]], na.rm = TRUE) / sqrt(sum(!is.na(ph[[a]]))),
      raw = mean(X[[ARMS2[[a]]]] <= .05, na.rm = TRUE), missing_alt = mean(is.na(X[[ARMS2[[a]]]])), missing_null = mean(is.na(Z[[ARMS2[[a]]]])))
  for (a in names(ARMS2)[-1]) { z <- pdiff(ph$DG, ph[[a]])
    pd2[[length(pd2) + 1]] <- data.frame(fam = cc$fam, C = cc$C, n = cc$n, arm = a, d = z[["d"]], se = z[["se"]], m = z[["m"]]) }
}
PW2 <- do.call(rbind, pw); PD2 <- do.call(rbind, pd2); PD2$z <- PD2$d / PD2$se
chk2 <- merge(PW2[PW2$arm == "DG", ], SP[SP$arm == "DG", ], by = c("fam", "C", "n"))
cat(sprintf("DG matched power recomputed vs niche_power.csv: max |diff| %.2e over %d cells\n", max(abs(chk2$matched.x - chk2$matched.y)), nrow(chk2)))
stored <- SP[SP$arm != "DG", c("fam", "C", "n", "arm", "R", "matched", "se", "raw")]
stored$missing_alt <- 0; stored$missing_null <- 0; stored$source <- "stored (niche_paper)"
PW2$source <- ifelse(PW2$arm == "DG", "recomputed (= stored)", "Tier 3")
T2 <- rbind(PW2, stored); T2 <- T2[order(T2$fam, T2$n, T2$C, -T2$matched), ]
utils::write.csv(T2, file.path(TAB, "T2_NICHE_POWER.csv"), row.names = FALSE)
utils::write.csv(PD2, file.path(TAB, "T2_NICHE_PAIRED.csv"), row.names = FALSE)
CNT <- do.call(rbind, lapply(split(PD2, PD2$arm), function(z) data.frame(arm = z$arm[1], cells = nrow(z),
  DG_above_2SE = sum(z$z > 2), arm_above_2SE = sum(z$z < -2), within_2SE = sum(abs(z$z) <= 2))))
utils::write.csv(CNT, file.path(TAB, "T2_NICHE_COUNTS.csv"), row.names = FALSE)
## size of every arm at n = 200, 500 (raw, 2,000 nulls)
S2 <- do.call(rbind, lapply(c(200, 500), function(n) { Z <- N[N$kind == 0 & N$n == n, ]
  do.call(rbind, lapply(names(ARMS2), function(a) { p <- Z[[ARMS2[[a]]]]; r <- mean(p <= .05, na.rm = TRUE)
    data.frame(n = n, arm = a, R = sum(!is.na(p)), size = r, se = sqrt(r * (1 - r) / sum(!is.na(p))), missing = mean(is.na(p))) })) }))
utils::write.csv(S2, file.path(TAB, "T2_NICHE_SIZE.csv"), row.names = FALSE)
cat("\nmatched power (new arms; stored arms in T2_NICHE_POWER.csv):\n")
print(format(stats::reshape(PW2[, c("fam", "C", "n", "arm", "matched")], idvar = c("fam", "n", "C"), timevar = "arm", direction = "wide"), digits = 3), row.names = FALSE)
cat("\ncells (of 12) where DG exceeds the arm by > 2 paired SE, and vice versa:\n"); print(CNT, row.names = FALSE)
cat("\nsizes at .05:\n"); print(format(S2, digits = 3), row.names = FALSE)

## =====================================================================================================
## 3. Held-out departure families
## =====================================================================================================
net_cols <- c(distributed = "p_dg", setNames(paste0("p_", NETS), NETS))
cellsL <- function(Ss) do.call(rbind, lapply(Ss, function(S) do.call(rbind, lapply(c(50, 100, 500), function(n)
  data.frame(design = "liu", cell = sprintf("S%d n=%d", S, n), S = S, n = n, fam = NA)))))
cellsN <- function(f) { z <- cells[cells$fam == f, ]; data.frame(design = "niche", cell = z$scen, S = 0, n = z$n, fam = f) }
FAM <- list(interaction = list(cells = rbind(cellsL(3:4), cellsN("inter")), nox = "NoInter"),
            quadratic = list(cells = rbind(cellsL(2), cellsN("quad")), nox = c("NoQuad", "NoQuadStrict")))
H <- list(); HP <- list(); HS <- list()
for (fm in names(FAM)) for (r in seq_len(nrow(FAM[[fm]]$cells))) { cc <- FAM[[fm]]$cells[r, ]
  if (cc$design == "liu") { b <- L[L$S == cc$S & L$n == cc$n, ]; nul <- b[b$kind == 0, ]; alt <- b[b$kind == 2, ]
    ph <- lapply(net_cols, function(v) phi_liu(nul[[v]], alt[[v]]))
  } else { alt <- N[N$scen == cc$cell, ]; nul <- N[N$kind == 0 & N$n == cc$n, ]
    ph <- lapply(net_cols, function(v) phi_niche(nul[[v]], alt[[v]])) }
  names(ph) <- names(net_cols)
  for (v in names(net_cols)) H[[length(H) + 1]] <- data.frame(family = fm, design = cc$design, cell = cc$cell, n = cc$n,
     net = v, variant = sub("_s[123]$", "", v), power = mean(ph[[v]]), size = mean(nul[[net_cols[[v]]]] <= .05))
  vmean <- function(vr) Reduce(`+`, ph[paste0(vr, "_s", 1:3)]) / 3
  for (vr in FAM[[fm]]$nox) { z <- pdiff(vmean("Control"), vmean(vr))
    HP[[length(HP) + 1]] <- data.frame(family = fm, design = cc$design, cell = cc$cell, n = cc$n, held_out = vr,
       control_minus_heldout = z[["d"]], se = z[["se"]], z = z[["d"]] / z[["se"]]) }
}
HC <- do.call(rbind, H); HPD <- do.call(rbind, HP)
utils::write.csv(HC, file.path(TAB, "T3_HELDOUT_CELLS.csv"), row.names = FALSE)
SUM <- do.call(rbind, lapply(split(HC, list(HC$family, HC$cell, HC$variant), drop = TRUE), function(z) data.frame(
  family = z$family[1], design = z$design[1], cell = z$cell[1], n = z$n[1], variant = z$variant[1],
  power_mean = mean(z$power), power_min = min(z$power), power_max = max(z$power),
  size_mean = mean(z$size), size_min = min(z$size), size_max = max(z$size))))
SUM <- SUM[(SUM$family == "interaction" & SUM$variant %in% c("distributed", "Control", "NoInter")) |
           (SUM$family == "quadratic" & SUM$variant %in% c("distributed", "Control", "NoQuad", "NoQuadStrict")), ]
SUM <- SUM[order(SUM$family, SUM$design, SUM$cell, SUM$variant), ]
utils::write.csv(SUM, file.path(TAB, "T3_HELDOUT_SUMMARY.csv"), row.names = FALSE)
utils::write.csv(HPD, file.path(TAB, "T3_HELDOUT_PAIRED.csv"), row.names = FALSE)
W <- stats::reshape(SUM[, c("family", "cell", "variant", "power_mean")], idvar = c("family", "cell"), timevar = "variant", direction = "wide")
cat("\n== 3. held-out families: size-adjusted power, mean over 3 seeds (range in T3_HELDOUT_SUMMARY.csv) ==\n")
print(format(W, digits = 3), row.names = FALSE)
cat("\nControl - held-out (seed-averaged, paired SE):\n"); print(format(HPD, digits = 3), row.names = FALSE)
OV <- do.call(rbind, lapply(split(SUM, list(SUM$family, SUM$design, SUM$variant), drop = TRUE), function(z) data.frame(
  family = z$family[1], design = z$design[1], variant = z$variant[1], cells = nrow(z), power = mean(z$power_mean),
  size_mean = mean(z$size_mean), size_min = min(z$size_min), size_max = max(z$size_max))))
OVP <- do.call(rbind, lapply(split(HPD, list(HPD$family, HPD$design, HPD$held_out), drop = TRUE), function(z) data.frame(
  family = z$family[1], design = z$design[1], held_out = z$held_out[1], cells = nrow(z), control_minus_heldout = mean(z$control_minus_heldout),
  se = sqrt(sum(z$se^2)) / nrow(z))))
OVP$z <- OVP$control_minus_heldout / OVP$se
utils::write.csv(OV, file.path(TAB, "T3_HELDOUT_OVERALL.csv"), row.names = FALSE)
utils::write.csv(OVP, file.path(TAB, "T3_HELDOUT_OVERALL_PAIRED.csv"), row.names = FALSE)
cat("\noverall (mean over the family's cells):\n"); print(format(OV, digits = 3), row.names = FALSE); print(format(OVP, digits = 3), row.names = FALSE)
## type I error of every network on every null set
TE <- rbind(do.call(rbind, lapply(2:4, function(S) do.call(rbind, lapply(c(50, 100, 500), function(n) { z <- L[L$S == S & L$n == n & L$kind == 0, ]
          data.frame(design = "liu", nulls = sprintf("S%d n=%d", S, n), net = names(net_cols), R = nrow(z),
                     size = vapply(net_cols, function(v) mean(z[[v]] <= .05), 0)) })))),
            do.call(rbind, lapply(c(200, 500), function(n) { z <- N[N$kind == 0 & N$n == n, ]
          data.frame(design = "niche", nulls = sprintf("n=%d", n), net = names(net_cols), R = nrow(z),
                     size = vapply(net_cols, function(v) mean(z[[v]] <= .05), 0)) })))
TE$se <- sqrt(TE$size * (1 - TE$size) / TE$R)
utils::write.csv(TE, file.path(TAB, "T3_NETWORK_SIZE.csv"), row.names = FALSE)
cat(sprintf("\nnetwork type I error at .05: range %.3f-%.3f over %d (network x null set); outside .05 +/- 2 SE: %d\n",
            min(TE$size), max(TE$size), nrow(TE), sum(abs(TE$size - .05) > 2 * sqrt(.05 * .95 / TE$R))))

## =====================================================================================================
## 4. Stukel fallback and Osius-Rojek (item C)
## =====================================================================================================
if (!is.null(C_rows)) {
  Cc <- C_rows; ok <- Cc$status == "ok"
  cat(sprintf("\n== 4. item C: %d rows | status: %s\n", nrow(Cc), paste(names(table(Cc$status)), table(Cc$status), collapse = ", ")))
  both <- ok & Cc$fallback == 0
  cat(sprintf("Stukel, both terms defined (%d data sets): max |SstBoth - stored| %.2e | stored NA there %d | max |SstBoth - LogisticDx| %.2e\n",
              sum(both), max(abs(Cc$p_sst_both[both] - Cc$p_stukel_stored[both]), na.rm = TRUE), sum(is.na(Cc$p_stukel_stored[both])),
              max(abs(Cc$p_sst_both[both] - Cc$p_sst_ldx[both]), na.rm = TRUE)))
  cat(sprintf("stored Stukel is NA on every data set with fallback > 0: %s (%d data sets); fallback = 2 (no term): %d\n",
              all(is.na(Cc$p_stukel_stored[ok & Cc$fallback > 0])), sum(ok & Cc$fallback > 0), sum(ok & Cc$fallback == 2, na.rm = TRUE)))
  FB <- do.call(rbind, lapply(split(Cc[ok, ], list(Cc$S[ok], Cc$n[ok])), function(z) data.frame(S = z$S[1], n = z$n[1], R = nrow(z),
    fallback_share = mean(z$fallback == 1), none_share = mean(z$fallback == 2), stored_na_share = mean(is.na(z$p_stukel_stored)))))
  FBS <- do.call(rbind, lapply(split(Cc[ok, ], Cc$S[ok]), function(z) data.frame(S = z$S[1], n = "all", R = nrow(z),
    fallback_share = mean(z$fallback == 1), none_share = mean(z$fallback == 2), stored_na_share = mean(is.na(z$p_stukel_stored)))))
  FB <- rbind(FB[order(FB$S, FB$n), ], FBS)
  ## size and size-adjusted power, per block: nulls 1001-2000 calibrate, all 1000 alternatives per C evaluated
  SP4 <- list()
  for (bk in unique(Cc$block)) { z <- Cc[Cc$block == bk & ok, ]; nul <- z$kind == 0 & z$rep > 1000
    for (v in c("p_stukel_fb", "p_stukel_stored")) { pw <- c()
      for (Cv in c(.25, .5, .75, 1)) pw <- c(pw, mean(phi_liu(z[[v]][nul], z[[v]][z$kind == 2 & z$C == Cv]), na.rm = TRUE))
      SP4[[length(SP4) + 1]] <- data.frame(block = bk, S = z$S[1], n = z$n[1], seed0 = z$seed0[1], version = v,
        size_all_nulls = mean(z[[v]][z$kind == 0] <= .05, na.rm = TRUE), size_heldout = mean(z[[v]][nul] <= .05, na.rm = TRUE),
        power_C0.25 = pw[1], power_C0.5 = pw[2], power_C0.75 = pw[3], power_C1 = pw[4], power_mean = mean(pw),
        missing = mean(is.na(z[[v]]))) } }
  SP4 <- do.call(rbind, SP4)
  SPA <- aggregate(cbind(size_all_nulls, power_mean, missing) ~ S + n + version, SP4, mean)
  SPB <- aggregate(cbind(size_all_nulls, power_mean, missing) ~ S + version, SP4, mean); SPB$n <- "all"
  SPA <- rbind(SPA, SPB[, names(SPA)])
  SPA <- stats::reshape(SPA, idvar = c("S", "n"), timevar = "version", direction = "wide"); SPA <- SPA[order(SPA$S, SPA$n), ]
  FB <- merge(FB, SPA, by = c("S", "n"), all.x = TRUE)
  utils::write.csv(FB, file.path(TAB, "T4_STUKEL_FALLBACK.csv"), row.names = FALSE)
  utils::write.csv(SP4, file.path(TAB, "T4_STUKEL_BY_BLOCK.csv"), row.names = FALSE)
  cat("\nStukel fallback share, size and size-adjusted power (mean over 2 blocks and 4 C; 'stored' = NA left out):\n")
  print(format(FB, digits = 3), row.names = FALSE)
  e1 <- abs(Cc$p_or_hand - Cc$p_or_stored); e2 <- abs(Cc$p_or_ldx - Cc$p_or_stored); e3 <- abs(Cc$p_or_hand - Cc$p_or_ldx)
  ORT <- data.frame(pair = c("hand vs stored (ebrahim.gof)", "LogisticDx vs stored", "hand vs LogisticDx"), R = c(sum(ok & !is.na(e1)), sum(ok & !is.na(e2)), sum(ok & !is.na(e3))),
    max_abs = c(max(e1[ok], na.rm = TRUE), max(e2[ok], na.rm = TRUE), max(e3[ok], na.rm = TRUE)),
    n_above_1e8 = c(sum(e1[ok] > 1e-8, na.rm = TRUE), sum(e2[ok] > 1e-8, na.rm = TRUE), sum(e3[ok] > 1e-8, na.rm = TRUE)),
    max_abs_when_min_pq_ge_1e4 = c(max(e1[ok & Cc$min_pq >= 1e-4], na.rm = TRUE), max(e2[ok & Cc$min_pq >= 1e-4], na.rm = TRUE), max(e3[ok & Cc$min_pq >= 1e-4], na.rm = TRUE)),
    p_or_NA = c(sum(is.na(Cc$p_or_hand[ok])), sum(is.na(Cc$p_or_ldx[ok])), NA))
  utils::write.csv(ORT, file.path(TAB, "T4_OSIUS_ROJEK.csv"), row.names = FALSE)
  big <- ok & !is.na(e1) & e1 > 1e-8
  if (any(big)) utils::write.csv(Cc[big, c("block", "kind", "C", "rep", "firth", "min_pq", "p_or_hand", "p_or_ldx", "p_or_stored")],
                                 file.path(TAB, "T4_OSIUS_ROJEK_DISCREPANT.csv"), row.names = FALSE)
  cat("\nOsius-Rojek agreement:\n"); print(format(ORT, digits = 3), row.names = FALSE)
  if (any(big)) cat(sprintf("discrepant rows (> 1e-8): %d; their largest min(p, 1-p) is %.2e; Firth fits among them: %d\n",
                            sum(big), max(Cc$min_pq[big]), sum(Cc$firth[big] == 1)))
}
cat("\ntiming (seconds per data set inside a worker, 22 workers): \n")
print(aggregate(cbind(t_combo, t_nets, t_grp) ~ design + n, A, mean))
