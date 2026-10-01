# support_external_split.R -- external validation on the SUPPORT data. The public release has no hospital or date
# identifier, so the split is random: half of the 8,873 complete cases develop a model, which is then frozen and
# checked on the other half. Repeated over 20 random splits. Models: M0 (linear in every covariate, the model of
# the application) and M2 (natural splines in the seven physiological variables). On the validation half:
#   DeepGOF-1 for frozen predictions, deepgof1.external(), all three readings from one call (B = 199)
#   the GiViTI calibration test in its external form (givitiR, devel = "external")
#   Hosmer-Lemeshow, ten groups of predicted risk, against chi-square(10)
#   Stukel's two terms added to the frozen linear predictor, likelihood ratio against chi-square(2)
#   the calibration intercept (calibration-in-the-large, slope fixed at 1) and the calibration slope
# Uses the package source with deepgof1.external() (load_all of the 2.9.0 tree).
# Output screen_SUPPORT/external_split.csv; the maps of split 1 in screen_SUPPORT/external_split1_maps.rds
suppressMessages({ devtools::load_all("C:/Users/ebrah/.cursor-tutor/projects/ebrahim.gof-proj-bagoft", quiet = TRUE)
                   library(splines) })
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/medical_example/data_prep.R")
out <- file.path(HERE, "screen_SUPPORT")
p2 <- prep(2L); D <- p2$data
F0 <- p2$formula
F2 <- stats::formula(readRDS(file.path(out, "M2.rds")))
X <- data.frame(age = D$age, meanbp = D$meanbp, hrt = D$hrt, resp = D$resp, temp = D$temp, crea = D$crea,
                sod = D$sod, wblc = D$wblc, num.co = D$num.co, male = as.numeric(D$sex == "male"))
NSPLIT <- as.integer(Sys.getenv("SPLITS", "20"))

one_split <- function(s) {
  set.seed(20261001L + s)
  dev <- sort(sample(nrow(D), floor(nrow(D) / 2))); val <- setdiff(seq_len(nrow(D)), dev)
  rows <- list(); maps <- list()
  for (m in c("M0", "M2")) {
    fit <- suppressWarnings(stats::glm(if (m == "M0") F0 else F2, data = D[dev, ], family = stats::binomial()))
    eta <- as.numeric(stats::predict(fit, newdata = D[val, ], type = "link")); p <- stats::plogis(eta)
    y <- D$y[val]; Xv <- X[val, ]
    set.seed(30261001L + s)
    t0 <- Sys.time()
    dg <- deepgof1.external(y, p, Xv, B = 199L, reading = "combined")
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    set.seed(40261001L + s)
    dga <- deepgof1.external(y, p, Xv, B = 199L, reading = "allpairs")
    giv <- tryCatch(givitiR::givitiCalibrationTest(o = y, e = p, devel = "external")$p.value, error = function(e) NA_real_)
    g <- pmin(ceiling(rank(p, ties.method = "first") / (length(p) / 10)), 10)
    O <- tapply(y, g, sum); E <- tapply(p, g, sum); ng <- tabulate(g, 10)
    hl <- stats::pchisq(sum((O - E)^2 / (E * (1 - E / ng))), 10, lower.tail = FALSE)
    dev0 <- -2 * sum(y * log(p) + (1 - y) * log(1 - p))
    za <- 0.5 * eta^2 * (p >= .5); zb <- -0.5 * eta^2 * (p < .5)
    gs <- suppressWarnings(stats::glm(y ~ 0 + za + zb, offset = eta, family = stats::binomial()))
    stk <- stats::pchisq(dev0 - gs$deviance, 2, lower.tail = FALSE)
    citl <- stats::coef(stats::glm(y ~ 1, offset = eta, family = stats::binomial()))[[1]]
    slope <- stats::coef(stats::glm(y ~ eta, family = stats::binomial()))[[2]]
    rows[[m]] <- data.frame(split = s, model = m, n_val = length(val), deaths = sum(y),
                            DG_axes = dg$components[["axes"]], DG_allpairs = dg$components[["allpairs"]],
                            DG_combined = dg$p.value, axes = paste(dg$axes, collapse = " & "),
                            top_pair = paste(dga$axes, collapse = " & "), GiViTI = giv, HL = hl, Stukel = stk,
                            citl = citl, slope = slope, dg_seconds = secs)
    if (s == 1L) maps[[m]] <- list(axes = dg, allpairs = dga)
  }
  if (s == 1L) saveRDS(maps, file.path(out, "external_split1_maps.rds"))
  do.call(rbind, rows)
}

cl <- parallel::makeCluster(as.integer(Sys.getenv("SP_WORKERS", "4")))
parallel::clusterExport(cl, c("D", "X", "F0", "F2", "out", "one_split"))
parallel::clusterEvalQ(cl, suppressMessages({
  devtools::load_all("C:/Users/ebrah/.cursor-tutor/projects/ebrahim.gof-proj-bagoft", quiet = TRUE); library(splines) }))
R <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(NSPLIT), one_split, chunk.size = 1L))
parallel::stopCluster(cl)
utils::write.csv(R, file.path(out, "external_split.csv"), row.names = FALSE)
print(R[R$split == 1, ], row.names = FALSE)
for (m in c("M0", "M2")) {
  z <- R[R$model == m, ]
  cat(sprintf("\n%s over %d splits: rejection at .05 -- DG combined %.2f, DG axes %.2f, DG all-pairs %.2f, GiViTI %.2f, HL %.2f, Stukel %.2f\n",
              m, nrow(z), mean(z$DG_combined <= .05), mean(z$DG_axes <= .05), mean(z$DG_allpairs <= .05),
              mean(z$GiViTI <= .05, na.rm = TRUE), mean(z$HL <= .05), mean(z$Stukel <= .05)))
  cat(sprintf("   median p: DG combined %.3f, GiViTI %.3f, HL %.3f, Stukel %.3f; calibration intercept %.3f [%.3f, %.3f], slope %.3f [%.3f, %.3f]\n",
              stats::median(z$DG_combined), stats::median(z$GiViTI, na.rm = TRUE), stats::median(z$HL), stats::median(z$Stukel),
              stats::median(z$citl), min(z$citl), max(z$citl), stats::median(z$slope), min(z$slope), max(z$slope)))
  print(sort(table(z$top_pair), decreasing = TRUE)[1:3])
}
