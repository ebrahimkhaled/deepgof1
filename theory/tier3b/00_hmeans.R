# 00_hmeans.R -- population mean of each departure shape h(x1, x2), x1, x2 iid N(0,1), by Monte Carlo
# (10 million draws, seed 20261024), with the exact value as a check. Also the share of h's variance the
# fitted linear terms cannot absorb (unweighted projection on 1, x1, x2), used only to space the pilot grid.
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3b/b_common.R")
dir.create(TABB, showWarnings = FALSE)
seed_at(20261024L); M <- 1e7
x1 <- stats::rnorm(M); x2 <- stats::rnorm(M)
exact <- c(threshold = 1 - stats::pnorm(1),
           bump = sqrt(0.5 / 2.5) * exp(-0.5^2 / 2.5) * sqrt(0.5 / 2.5) * exp(-1^2 / 2.5),
           oscillation = 0, corner = (1 - stats::pnorm(0.8))^2)
out <- do.call(rbind, lapply(SHAPES, function(s) {
  h <- h_raw(s, x1, x2); m <- mean(h)
  idx <- 1:1e6; hc <- h[idx] - m; r <- stats::lm.fit(cbind(1, x1[idx], x2[idx]), hc)$residuals
  data.frame(shape = s, mc_mean = m, mc_se = stats::sd(h) / sqrt(M), exact_mean = exact[[s]],
             sd_h = stats::sd(h), sd_resid_linear = stats::sd(r))
}))
print(out, digits = 6)
utils::write.csv(out, file.path(TABB, "H_MEANS.csv"), row.names = FALSE)
