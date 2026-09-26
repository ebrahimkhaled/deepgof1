# =============================================================================
# analyze_addon.R -- join the add-on rows to the stored rows and tabulate raw rejection at .05.
#
# 1. Pairing, on EVERY row: the add-on HL_check must equal the stored HL (|diff| <= 1e-12).
#    Any mismatch is printed and the script stops -- no rate is computed from unpaired rows.
# 2. Raw rejection rate at .05 with its binomial SE, per cell, for mHL, GiViTI, Stukel (both
#    forms) and, on the niche, BAGofT where its rows exist. NA p-values count as non-rejections
#    and the NA share is reported.
#    Liu: the join is an INNER join on (kind, rep, C), so only datasets the stored block kept are
#    used; the two corpora (101000, 102000) are pooled per (setting, n, C), as in the paper.
#    Niche: joined on (scen, rep).
# Writes: liu_addon_rates.csv, niche_addon_rates.csv (+ console).  Run: Rscript analyze_addon.R
# =============================================================================
HERE <- "/path/to/project/paper_deepgof/theory/addon_arms"
RES  <- "/path/to/project/paper_deepgof/release/deepgof1/results/liu"
NRO  <- "/path/to/project/paper_deepgof/theory/niche/rows"
A    <- 0.05
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-12 * pmax(1, abs(b)))
rate <- function(p) { r <- mean(!is.na(p) & p <= A); c(rate = r, se = sqrt(r * (1 - r) / length(p)),
                                                        na = mean(is.na(p)), N = length(p)) }

## ---- Liu grid ----
fs <- list.files(file.path(HERE, "liu_rows"), pattern = "^S[0-9]_n[0-9]+_s[0-9]+\\.csv$", full.names = TRUE)
if (length(fs)) {
  L <- do.call(rbind, lapply(fs, function(f) {
    b  <- sub("\\.csv$", "", basename(f)); ad <- utils::read.csv(f)
    st <- utils::read.csv(file.path(RES, b, "blk_pvalues.csv"), check.names = FALSE)[, c("kind", "rep", "C", "HL", "Stukel")]
    ad$Ckey <- ifelse(is.na(ad$C), -1, ad$C); st$Ckey <- ifelse(is.na(st$C), -1, st$C)
    m <- merge(st, ad[, setdiff(names(ad), "C")], by = c("kind", "rep", "Ckey"))
    nst <- sum(st$kind != 1); njn <- sum(m$kind != 1)
    if (njn != nst) cat(sprintf("  %s: %d stored rows (kinds 0/2) but %d joined\n", b, nst, njn))
    bad <- !same(m$HL_check, m$HL)
    if (any(bad)) { print(utils::head(m[bad, ])); stop("HL pairing FAILED in ", b) }
    m$block <- b; m
  }))
  cat(sprintf("Liu: %d blocks, %d joined rows, HL identical on every row\n", length(fs), nrow(L)))
  L <- L[L$kind != 1, ]; L$C[L$kind == 0] <- 0
  arms <- c("mHL", "GiViTI", "Stukel_joint", "Stukel_marginal")
  out <- list()
  for (g in split(L, list(L$setting, L$n, L$C), drop = TRUE)) for (a in arms)
    out[[length(out) + 1L]] <- data.frame(setting = g$setting[1], n = g$n[1], C = g$C[1], arm = a, t(rate(g[[a]])))
  R <- do.call(rbind, out); R <- R[order(R$arm, R$setting, R$n, R$C), ]
  utils::write.csv(R, file.path(HERE, "liu_addon_rates.csv"), row.names = FALSE)
  print(R[R$C == 0, ], row.names = FALSE, digits = 3)
}

## ---- niche ----
fs <- list.files(file.path(HERE, "niche_rows"), pattern = "\\.csv$", full.names = TRUE)
if (length(fs)) {
  ad <- do.call(rbind, lapply(fs, utils::read.csv))
  st <- do.call(rbind, lapply(list.files(NRO, pattern = "\\.csv$", full.names = TRUE), utils::read.csv))
  m  <- merge(st[, c("scen", "rep", "HL", "Stukel")], ad, by = c("scen", "rep"))
  bad <- !same(m$HL_check, m$HL)
  if (any(bad)) { print(utils::head(m[bad, ])); stop("HL pairing FAILED on the niche") }
  cat(sprintf("\nNiche: %d joined rows of %d stored, HL identical on every row\n", nrow(m), nrow(st)))
  bf <- list.files(file.path(HERE, "niche_rows_bagoft"), pattern = "\\.csv$", full.names = TRUE)
  if (length(bf)) {
    bg <- do.call(rbind, lapply(bf, utils::read.csv))
    m  <- merge(m, data.frame(scen = bg$scen, rep = bg$rep, BAGofT = bg$p), by = c("scen", "rep"), all.x = TRUE)
  }
  arms <- intersect(c("mHL", "GiViTI", "Stukel_joint", "Stukel_marginal", "BAGofT"), names(m))
  out <- list()
  for (g in split(m, m$scen)) for (a in arms) {
    p <- g[[a]]; if (a == "BAGofT") p <- p[!is.na(g$BAGofT)]      # BAGofT: only the reps it has run
    if (length(p)) out[[length(out) + 1L]] <- data.frame(scen = g$scen[1], arm = a, t(rate(p)))
  }
  R <- do.call(rbind, out)
  utils::write.csv(R, file.path(HERE, "niche_addon_rates.csv"), row.names = FALSE)
  print(R, row.names = FALSE, digits = 3)
}
