# 11_choose_C.R -- choose three departure sizes C per (design, shape) from the pilot (n = 316, 100 data sets per
# grid point): DeepGOF-1 raw power (p <= .05) targets .2, .5, .8. The pilot power curve is made monotone in C by
# isotonic regression (pooled adjacent violators, weights = data sets), then C is interpolated linearly in log C.
# When the isotonic curve never reaches a target (power saturates: an indicator departure is bounded however
# large C is), the largest grid C is used and the row is flagged. C is rounded to 3 significant digits.
# Output: tables/PILOT_POWER.csv (every grid point, every arm's raw power) and tables/C_CHOSEN.csv.
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3b/b_common.R")
R <- do.call(rbind, lapply(list.files(file.path(T3B, "rows_pilot"), "\\.csv$", full.names = TRUE), utils::read.csv, stringsAsFactors = FALSE))
stopifnot(!anyDuplicated(R[, c("design", "shape", "ci", "rep")]))
arms <- c("p_dg", "p_ss", "p_score", "p_score_asym", "p_hl", "p_stukel_fb")
PP <- do.call(rbind, lapply(split(R, list(R$design, R$shape, R$ci), drop = TRUE), function(z) {
  o <- data.frame(design = z$design[1], shape = z$shape[1], ci = z$ci[1], C = z$C[1], R = nrow(z))
  for (a in arms) o[[sub("^p_", "raw_", a)]] <- mean(z[[a]] <= .05, na.rm = TRUE)
  o }))
PP <- PP[order(PP$design, PP$shape, PP$C), ]
utils::write.csv(PP, file.path(TABB, "PILOT_POWER.csv"), row.names = FALSE)
TG <- c(.2, .5, .8)
CC <- do.call(rbind, lapply(split(PP, list(PP$design, PP$shape), drop = TRUE), function(z) {
  ## corner: the C > 0 grid saturates (power <= .08 at every C); the C < 0 grid is used (curve in |C|)
  if (z$shape[1] == "corner" && any(z$C < 0)) z <- z[z$C < 0, ] else z <- z[z$C > 0, ]
  sgn <- sign(z$C[1]); z$C <- abs(z$C)
  z <- z[order(z$C), ]; iso <- stats::isoreg(log(z$C), z$raw_dg)$yf
  do.call(rbind, lapply(seq_along(TG), function(k) { t <- TG[k]
    flag <- ""
    if (t <= min(iso)) { Cv <- min(z$C); flag <- "target below the smallest grid C" }
    else if (t > max(iso)) { Cv <- max(z$C); flag <- sprintf("saturates: isotonic max %.2f < target", max(iso)) }
    else { i <- which(iso >= t)[1]
      if (iso[i] == t || i == 1) Cv <- z$C[i] else {
        lo <- i - 1L; while (lo > 1 && iso[lo] == iso[i]) lo <- lo - 1L
        w <- (t - iso[lo]) / (iso[i] - iso[lo]); Cv <- exp(log(z$C[lo]) + w * (log(z$C[i]) - log(z$C[lo]))) } }
    data.frame(design = z$design[1], shape = z$shape[1], ci = k, target_dg_power = t, C = sgn * signif(Cv, 3), flag = flag) }))
}))
CC <- CC[order(CC$design, CC$shape, CC$ci), ]
print(PP, digits = 3, row.names = FALSE); cat("\n"); print(CC, row.names = FALSE)
utils::write.csv(CC, file.path(TABB, "C_CHOSEN.csv"), row.names = FALSE)
