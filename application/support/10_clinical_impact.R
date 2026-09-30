# stage1_impact_SUPPORT.R -- PLAN.md Sec. 4 step 4: what the map-guided repairs change clinically.
# Ladder: M0 (linear) -> M1 ns(meanbp) -> M2 ns(all physiology) -> M3 (+ ns(age)).
# (i) likelihood-ratio tests between rungs, an adjudicator independent of every GOF test;
# (ii) apparent and optimism-corrected (200 bootstrap refits) AUC, Brier score and calibration slope;
# (iii) predicted in-hospital death risk for patient profiles, with 95% intervals (delta method, logit scale).
suppressMessages({ library(splines) })
source("data_prep.R")
p <- prep(2L); D <- p$data
out <- file.path(HERE, "results", "screen_SUPPORT")
L <- list(M0 = readRDS(file.path(out, "M0.rds")), M1 = readRDS(file.path(out, "M1_b.rds"))$fit,
          M2 = readRDS(file.path(out, "M2.rds")), M3 = readRDS(file.path(out, "M3.rds")))

cat("== (i) likelihood-ratio tests\n")
for (i in 2:4) {
  a <- stats::anova(L[[i - 1]], L[[i]], test = "LRT")
  cat(sprintf("%s vs %s: LR = %.1f on %d df, p = %.2g\n", names(L)[i - 1], names(L)[i],
              a$Deviance[2], a$Df[2], a$`Pr(>Chi)`[2]))
}

auc <- function(y, s) { r <- rank(s); n1 <- sum(y); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }
perf <- function(fit, dat) {
  lp <- stats::predict(fit, newdata = dat); ph <- stats::plogis(lp)
  c(AUC = auc(dat$y, ph), Brier = mean((dat$y - ph)^2),
    slope = unname(stats::coef(stats::glm(dat$y ~ lp, family = stats::binomial()))[2]))
}
cat("\n== (ii) performance: apparent, optimism, corrected (200 bootstrap refits)\n")
set.seed(20261940L)
BS <- replicate(200, sample(nrow(D), replace = TRUE), simplify = FALSE)
perf_tab <- data.frame()
for (m in names(L)) {
  app <- perf(L[[m]], D)
  opt <- rowMeans(sapply(BS, function(i) {
    fb <- suppressWarnings(stats::glm(stats::formula(L[[m]]), data = D[i, ], family = stats::binomial()))
    perf(fb, D[i, ]) - perf(fb, D)
  }))
  perf_tab <- rbind(perf_tab, data.frame(model = m, AUC = app[1], AUC_corr = app[1] - opt[1],
                                         Brier = app[2], Brier_corr = app[2] - opt[2],
                                         slope_corr = app[3] - opt[3]))
}
print(format(perf_tab, digits = 4), row.names = FALSE)
utils::write.csv(perf_tab, file.path(out, "stage1_performance.csv"), row.names = FALSE)

cat("\n== (iii) predicted in-hospital death risk (95% CI) for patient profiles\n")
base <- as.list(sapply(D[, c("age", "hrt", "resp", "temp", "crea", "sod", "wblc", "num.co")], stats::median))
mk <- function(...) { x <- modifyList(base, list(...)); x$sex <- factor("male", levels = levels(D$sex)); as.data.frame(x) }
prof <- rbind(
  cbind(profile = "typical patient (all medians, MAP 77)", mk(meanbp = 77)),
  cbind(profile = "hypotension, MAP 45", mk(meanbp = 45)),
  cbind(profile = "hypotension, MAP 45, 3 comorbidities", mk(meanbp = 45, num.co = 3)),
  cbind(profile = "hypertension, MAP 130", mk(meanbp = 130)),
  cbind(profile = "hypertension, MAP 130, 3 comorbidities", mk(meanbp = 130, num.co = 3)),
  cbind(profile = "tachypnoea, resp 40, MAP 77", mk(meanbp = 77, resp = 40)))
pr <- function(fit) {
  z <- stats::predict(fit, newdata = prof, se.fit = TRUE)
  sprintf("%.2f (%.2f-%.2f)", stats::plogis(z$fit), stats::plogis(z$fit - 1.96 * z$se.fit),
          stats::plogis(z$fit + 1.96 * z$se.fit))
}
tab <- data.frame(profile = prof$profile, M0 = pr(L$M0), M1 = pr(L$M1), M3 = pr(L$M3))
print(tab, row.names = FALSE)
utils::write.csv(tab, file.path(out, "stage1_profiles.csv"), row.names = FALSE)
cat("\nmedians used:", paste(names(base), round(unlist(base), 2), collapse = ", "), "\n")
