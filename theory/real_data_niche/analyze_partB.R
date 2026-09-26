# =============================================================================
# analyze_partB.R -- the pre-declared analysis of Part B (PREDECLARATION.md (d)).
# NA -> p = 1. Size: raw rejection at .05 on the 1,000 nulls of each (data, design); calibrated if
# size <= .05 + 3 sqrt(.05 .95 / R). Matched level: matched_phi() copied VERBATIM from
# theory/niche/analyze_niche.R, each arm on its own nulls of the same (data, design).
# Paired difference DG minus arm: mean(phi_DG - phi_arm), SE = sd(.)/sqrt(R).
# Writes summary_size.csv, summary_power.csv, summary_diff.csv, result_tables.md.
# =============================================================================
RD <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/real_data_niche"
ROWD <- file.path(RD, Sys.getenv("PB_ROWS", "rows"))
CMP  <- c("HL", "Stukel_marginal", "Stukel_joint", "StuteZhu", "GiViTI", "mHL", "OsiusRojek", "proj")
ARMS <- c("DG", CMP)
ALPHA <- 0.05
fs <- list.files(ROWD, pattern = "\\.csv$", full.names = TRUE)
D <- do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
D <- D[order(D$id, D$rep), ]
if (anyDuplicated(D[, c("id", "rep")])) stop("duplicate (id, rep) rows")
na_rate <- aggregate(D[, ARMS], list(scen = D$scen), function(v) mean(!is.finite(v)))
for (a in ARMS) D[[a]][!is.finite(D[[a]])] <- 1

## ---- VERBATIM from theory/niche/analyze_niche.R ----
matched_phi <- function(pn, pa, alpha = ALPHA) {           # returns per-replicate phi for pa
  t <- sort(unique(pn)); F <- findInterval(t, sort(pn)) / length(pn)   # = searchsorted(side = right)/n
  ok <- which(F <= alpha)
  if (!length(ok)) { g <- if (F[1] > 0) alpha / F[1] else 0; return(list(phi = g * (pa <= t[1]), t0 = t[1], g = g)) }
  i <- max(ok); t0 <- t[i]; F0 <- F[i]
  if (i + 1 > length(t)) return(list(phi = as.numeric(pa <= t0), t0 = t0, g = 0))
  t1 <- t[i + 1]; F1 <- F[i + 1]; g <- if (F1 <= F0) 0 else (alpha - F0) / (F1 - F0)
  list(phi = as.numeric(pa <= t0) + g * as.numeric(pa > t0 & pa <= t1), t0 = t0, g = g)
}

DS <- unique(D[, c("data", "design")])
SZ <- do.call(rbind, lapply(seq_len(nrow(DS)), function(i) {
  Z <- D[D$data == DS$data[i] & D$design == DS$design[i] & D$fam == "null", ]
  do.call(rbind, lapply(ARMS, function(a) { r <- mean(Z[[a]] <= ALPHA)
    data.frame(data = DS$data[i], design = DS$design[i], arm = a, R = nrow(Z), size = r, se = sqrt(r * (1 - r) / nrow(Z)),
               calibrated = r <= ALPHA + 3 * sqrt(ALPHA * (1 - ALPHA) / nrow(Z))) }))
}))

cells <- unique(D[D$fam != "null", c("data", "design", "fam", "level", "C")])
cells <- cells[order(cells$design != "niche", cells$data, cells$fam, cells$level), ]
PW <- list(); DF <- list()
for (k in seq_len(nrow(cells))) {
  cc <- cells[k, ]
  A <- D[D$data == cc$data & D$design == cc$design & D$fam == cc$fam & D$level == cc$level, ]
  N0 <- D[D$data == cc$data & D$design == cc$design & D$fam == "null", ]
  PHI <- sapply(ARMS, function(a) matched_phi(N0[[a]], A[[a]])$phi)
  RAW <- sapply(ARMS, function(a) as.numeric(A[[a]] <= ALPHA))
  R <- nrow(A)
  for (a in ARMS) PW[[length(PW) + 1]] <- data.frame(cc, arm = a, R = R,
    matched = mean(PHI[, a]), matched_se = stats::sd(PHI[, a]) / sqrt(R), raw = mean(RAW[, a]), raw_se = stats::sd(RAW[, a]) / sqrt(R))
  for (a in CMP) {
    dm <- PHI[, "DG"] - PHI[, a]; dr <- RAW[, "DG"] - RAW[, a]
    DF[[length(DF) + 1]] <- data.frame(cc, arm = a, R = R,
      d_matched = mean(dm), se_matched = stats::sd(dm) / sqrt(R), z_matched = mean(dm) / (stats::sd(dm) / sqrt(R)),
      d_raw = mean(dr), se_raw = stats::sd(dr) / sqrt(R), z_raw = mean(dr) / (stats::sd(dr) / sqrt(R)))
  }
}
PW <- do.call(rbind, PW); DF <- do.call(rbind, DF)
DF$z_matched[!is.finite(DF$z_matched)] <- 0; DF$z_raw[!is.finite(DF$z_raw)] <- 0

## ---- the pre-declared reading ----
verdict <- function(design, zcol, datas = c("myopia", "Pima")) {
  X <- DF[DF$design == design & DF$data %in% datas, ]; ck <- unique(X[, c("data", "fam", "level", "C")])
  lead <- behind <- logical(nrow(ck)); nb <- integer(nrow(ck))
  for (k in seq_len(nrow(ck))) {
    Y <- X[X$data == ck$data[k] & X$fam == ck$fam[k] & X$level == ck$level[k], ]
    cal <- SZ$arm[SZ$data == ck$data[k] & SZ$design == design & SZ$calibrated]
    Pc <- Y[Y$arm %in% cal, ]
    lead[k] <- nrow(Pc) > 0 && all(Pc[[zcol]] > 2)
    behind[k] <- any(Y[[zcol]] < -2); nb[k] <- sum(Y[[zcol]] < -2)
  }
  dg_cal <- all(SZ$calibrated[SZ$arm == "DG" & SZ$design == design & SZ$data %in% datas])
  need <- ceiling(nrow(ck) / 2)
  v <- if (dg_cal && sum(lead) >= need && !any(behind)) "NICHE CONFIRMED on real designs" else if (sum(behind) >= need) "NOT CONFIRMED -- NEGATIVE" else "NOT CONFIRMED"
  list(verdict = v, dg_calibrated = dg_cal, n_cells = nrow(ck), n_lead = sum(lead), n_behind = sum(behind),
       cells = cbind(ck, lead = lead, behind = behind, n_arms_above = nb))
}
V <- list(primary_niche_matched = verdict("niche", "z_matched"),
          niche_myopia_matched = verdict("niche", "z_matched", "myopia"),
          niche_Pima_matched = verdict("niche", "z_matched", "Pima"),
          sensitivity_niche_raw = verdict("niche", "z_raw"),
          dense_control_matched = verdict("dense", "z_matched"))

utils::write.csv(SZ, file.path(RD, "summary_size.csv"), row.names = FALSE)
utils::write.csv(PW, file.path(RD, "summary_power.csv"), row.names = FALSE)
utils::write.csv(DF, file.path(RD, "summary_diff.csv"), row.names = FALSE)

f3 <- function(x) sprintf("%.3f", x)
L <- c(sprintf("Rows: %d; scenarios: %d", nrow(D), length(unique(D$scen))), "",
       "### Size at .05 (raw, binomial SE); * = not calibrated (> .05 + 3 SE)", "",
       paste0("| data | design | ", paste(ARMS, collapse = " | "), " |"), paste0("|---|---|", paste(rep("---", length(ARMS)), collapse = "|"), "|"))
for (i in seq_len(nrow(DS))) { z <- SZ[SZ$data == DS$data[i] & SZ$design == DS$design[i], ]; z <- z[match(ARMS, z$arm), ]
  L <- c(L, paste0("| ", DS$data[i], " | ", DS$design[i], " | ", paste(sprintf("%s (%s)%s", f3(z$size), f3(z$se), ifelse(z$calibrated, "", "*")), collapse = " | "), " |")) }
for (what in c("matched", "raw")) {
  L <- c(L, "", sprintf("### Power, %s (MC SE)", if (what == "matched") "MATCHED level (randomized to exactly .05 on own nulls)" else "RAW (p <= .05)"), "",
         paste0("| data | design | family | C | ", paste(ARMS, collapse = " | "), " |"), paste0("|---|---|---|---|", paste(rep("---", length(ARMS)), collapse = "|"), "|"))
  for (k in seq_len(nrow(cells))) { cc <- cells[k, ]; z <- PW[PW$data == cc$data & PW$design == cc$design & PW$fam == cc$fam & PW$level == cc$level, ]; z <- z[match(ARMS, z$arm), ]
    L <- c(L, paste0("| ", cc$data, " | ", cc$design, " | ", cc$fam, " | ", cc$C, " | ", paste(sprintf("%s (%s)", f3(z[[what]]), f3(z[[paste0(what, "_se")]])), collapse = " | "), " |")) }
}
for (what in c("matched", "raw")) {
  L <- c(L, "", sprintf("### DG minus arm, %s: difference (paired SE) [z]", what), "",
         paste0("| data | design | family | C | ", paste(CMP, collapse = " | "), " |"), paste0("|---|---|---|---|", paste(rep("---", length(CMP)), collapse = "|"), "|"))
  for (k in seq_len(nrow(cells))) { cc <- cells[k, ]; z <- DF[DF$data == cc$data & DF$design == cc$design & DF$fam == cc$fam & DF$level == cc$level, ]; z <- z[match(CMP, z$arm), ]
    L <- c(L, paste0("| ", cc$data, " | ", cc$design, " | ", cc$fam, " | ", cc$C, " | ",
      paste(sprintf("%+.3f (%s) [%+.1f]", z[[paste0("d_", what)]], f3(z[[paste0("se_", what)]]), z[[paste0("z_", what)]]), collapse = " | "), " |")) }
}
dg <- aggregate(cbind(sel_ab, firth, events, secs) ~ scen, data = D, FUN = mean)
L <- c(L, "", "### Diagnostics per scenario (means): sel_ab = DeepGOF-1 selected the pair (a, b) carrying the departure", "",
       "| scen | sel_ab | firth | events | secs/rep |", "|---|---|---|---|---|",
       sprintf("| %s | %.3f | %.4f | %.1f | %.2f |", dg$scen, dg$sel_ab, dg$firth, dg$events, dg$secs),
       "", "### NA rate per arm per scenario (NA -> p = 1)", "",
       paste0("| scen | ", paste(ARMS, collapse = " | "), " |"), paste0("|---|", paste(rep("---", length(ARMS)), collapse = "|"), "|"),
       sprintf("| %s | %s |", na_rate$scen, apply(na_rate[, ARMS], 1, function(r) paste(sprintf("%.3f", r), collapse = " | "))))
for (nm in names(V)) { v <- V[[nm]]
  L <- c(L, "", sprintf("### Verdict [%s]: **%s**", nm, v$verdict),
         sprintf("DG calibrated on the nulls: %s; cells: %d; DG above every calibrated arm by > 2 SE: %d; some arm above DG by > 2 SE: %d",
                 v$dg_calibrated, v$n_cells, v$n_lead, v$n_behind), "",
         "| data | family | level | C | lead | behind | # arms > 2 SE above DG |", "|---|---|---|---|---|---|---|",
         sprintf("| %s | %s | %s | %s | %s | %s | %s |", v$cells$data, v$cells$fam, v$cells$level, v$cells$C, v$cells$lead, v$cells$behind, v$cells$n_arms_above)) }
writeLines(L, file.path(RD, "result_tables.md"))
cat(paste(L, collapse = "\n"), "\n")
