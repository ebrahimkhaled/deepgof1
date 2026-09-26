# =============================================================================
# analyze_proj.R -- joins the proj_addon.R rows to the stored proj p-values and reports
#   (1) pairing: proj_old (recomputed) must equal the stored proj on EVERY joined row;
#   (2) rejection at .05 per cell, with binomial SE, of three arms on the same data and bootstrap
#       draws: the stored (deployed) proj, proj_liu_raw (Liu et al.'s statistic; PRIMARY) and
#       proj_liu_std (same weight, standardised covariates; secondary), plus paired discordant counts.
# Liu grid: stored proj = release/deepgof1/results/liu/<block>/augment_pvalues.csv, key
#   (kind, rep, C). niche: stored proj = theory/niche/rows/*.csv, key (scen, rep).
# Reads chunk CSVs (finished) only; .part files are ignored.
# Run:  Rscript analyze_proj.R [liu|niche|both] [rows folder suffix, default none]
# Output: proj_pairing_<bench>.csv (every joined row), proj_rates_<bench>.csv, console summary.
# =============================================================================
HERE <- "/path/to/project/paper_deepgof/theory/addon_arms"
RES  <- "/path/to/project/paper_deepgof/release/deepgof1/results/liu"
NRO  <- "/path/to/project/paper_deepgof/theory/niche/rows"
args <- commandArgs(TRUE); WHICH <- if (length(args)) args[1] else "both"
SUF  <- if (length(args) > 1) args[2] else ""
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-12)
rate <- function(p) { p <- p[!is.na(p)]; r <- mean(p <= 0.05); c(rate = r, se = sqrt(r * (1 - r) / length(p)), m = length(p)) }
## one cell: the three arms (stored deployed proj; raw-X Liu = primary; standardised Liu = secondary)
## and the paired discordant counts of the primary arm against each of the other two
three <- function(s) {
  a <- rate(s$proj); r <- rate(s$proj_liu_raw); z <- rate(s$proj_liu_std)
  d <- function(x, y) sum(x <= .05 & y > .05, na.rm = TRUE)
  data.frame(m = a[["m"]], proj_stored = a[["rate"]], se_stored = a[["se"]],
             proj_liu_raw = r[["rate"]], se_raw = r[["se"]], proj_liu_std = z[["rate"]], se_std = z[["se"]],
             raw_only_vs_stored = d(s$proj_liu_raw, s$proj), stored_only_vs_raw = d(s$proj, s$proj_liu_raw),
             raw_only_vs_std = d(s$proj_liu_raw, s$proj_liu_std), std_only_vs_raw = d(s$proj_liu_std, s$proj_liu_raw))
}

if (WHICH %in% c("liu", "both")) {
  D <- file.path(HERE, paste0("proj_rows_liu", SUF))
  fs <- list.files(D, pattern = "^chunk_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(fs)) {
    R <- do.call(rbind, lapply(fs, utils::read.csv))
    J <- list()
    for (b in unique(paste0("S", R$setting, "_n", R$n, "_s", R$seed0))) {
      st <- utils::read.csv(file.path(RES, b, "augment_pvalues.csv"))
      r  <- R[paste0("S", R$setting, "_n", R$n, "_s", R$seed0) == b, ]
      r$Ck <- ifelse(r$kind == 2, r$C, NA); st$Ck <- ifelse(st$kind == 2, st$C, NA)
      m <- merge(r, st[, c("kind", "rep", "Ck", "proj")], by = c("kind", "rep", "Ck"), all.x = TRUE)
      m$block <- b; J[[b]] <- m
    }
    J <- do.call(rbind, J)
    ok <- J$status == "ok" & !is.na(J$proj); pk <- ok & !is.na(J$proj_old)
    cat(sprintf("LIU pairing: %d rows | %d joined with a stored proj | of those with proj_old recomputed, proj_old == stored on %d / %d (max |diff| %.2e) | dropped here %d, not in stored %d\n",
                nrow(J), sum(ok), sum(same(J$proj_old[pk], J$proj[pk])), sum(pk),
                max(abs(J$proj_old[pk] - J$proj[pk]), 0), sum(J$status != "ok"), sum(J$status == "ok" & is.na(J$proj))))
    utils::write.csv(J, file.path(HERE, paste0("proj_pairing_liu", SUF, ".csv")), row.names = FALSE)
    cells <- unique(J[ok, c("setting", "n", "kind", "C")])
    out <- do.call(rbind, lapply(seq_len(nrow(cells)), function(k) {
      cc <- cells[k, ]; s <- J[ok & J$setting == cc$setting & J$n == cc$n & J$kind == cc$kind &
                                (is.na(cc$C) | (!is.na(J$C) & J$C == cc$C)), ]
      cbind(data.frame(setting = cc$setting, n = cc$n, kind = cc$kind, C = cc$C), three(s))
    }))
    out <- out[order(out$setting, out$n, out$kind, out$C), ]
    utils::write.csv(out, file.path(HERE, paste0("proj_rates_liu", SUF, ".csv")), row.names = FALSE)
    print(out, row.names = FALSE, digits = 3)
  }
}

if (WHICH %in% c("niche", "both")) {
  D <- file.path(HERE, paste0("proj_rows_niche", SUF))
  fs <- list.files(D, pattern = "\\.csv$", full.names = TRUE)
  if (length(fs)) {
    R  <- do.call(rbind, lapply(fs, utils::read.csv))
    st <- do.call(rbind, lapply(list.files(NRO, pattern = "\\.csv$", full.names = TRUE), function(f)
      utils::read.csv(f)[, c("scen", "rep", "proj", "DG")]))
    J  <- merge(R, st, by = c("scen", "rep"), all.x = TRUE)
    ok <- J$status == "ok" & !is.na(J$proj); pk <- ok & !is.na(J$proj_old)
    cat(sprintf("NICHE pairing: %d rows | %d joined with a stored proj | of those with proj_old recomputed, proj_old == stored on %d / %d (max |diff| %.2e) | stored DG NA on joined rows %d\n",
                nrow(J), sum(ok), sum(same(J$proj_old[pk], J$proj[pk])), sum(pk),
                max(abs(J$proj_old[pk] - J$proj[pk]), 0), sum(ok & is.na(J$DG))))
    utils::write.csv(J, file.path(HERE, paste0("proj_pairing_niche", SUF, ".csv")), row.names = FALSE)
    out <- do.call(rbind, lapply(split(J[ok, ], J$scen[ok]), function(s) {
      cbind(data.frame(scen = s$scen[1]), three(s))
    }))
    utils::write.csv(out, file.path(HERE, paste0("proj_rates_niche", SUF, ".csv")), row.names = FALSE)
    print(out, row.names = FALSE, digits = 3)
  }
}
