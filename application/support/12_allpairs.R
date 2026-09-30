# support_allpairs.R -- the SUPPORT models read by ebrahim.gof 2.9.0's deepgof1(), both readings:
#   v1        (the default) one map over the two columns with the largest |b| sd
#   allpairs  the largest score over every pair of covariates, the same maximum in every replicate
# 2.9.0 refits the bootstrap on the design matrix. Up to 2.8.0 the refit re-evaluated the formula on the
# model frame, which fails for ns(meanbp, 3): every replicate scored +Inf and the repairs with a spline got
# p = 1 (stage1_repair_SUPPORT.log). This script replaces those numbers.
# Models: M0 (linear), the four repairs M1_a..M1_d (M1_b = M1 of the paper), M2, M3. Same seed as the
# screen (DG_SEED = 20261929 + 2). M0 is also read by all-pairs over 20 further tie-break seeds.
# Needs ebrahim.gof 2.9.0 or later.
suppressMessages({ library(ebrahim.gof); library(splines) })
stopifnot(utils::packageVersion("ebrahim.gof") >= "2.9.0")
source("data_prep.R")
out <- file.path(HERE, "results", "screen_SUPPORT")
DG_SEED <- 20261929L + 2L
fits <- list(M0 = readRDS(file.path(out, "M0.rds")))
for (r in c("a", "b", "c", "d")) fits[[paste0("M1_", r)]] <- readRDS(file.path(out, paste0("M1_", r, ".rds")))$fit
fits$M2 <- readRDS(file.path(out, "M2.rds")); fits$M3 <- readRDS(file.path(out, "M3.rds"))

res <- list(); top <- list()
for (m in names(fits)) for (rd in c("v1", "allpairs")) {
  t0 <- Sys.time(); set.seed(DG_SEED)
  dg <- deepgof1(fits[[m]], B = 199, reading = rd)
  stopifnot(all(is.finite(dg$boot)))                # the refit bug would show here as +Inf
  res[[length(res) + 1]] <- data.frame(model = m, reading = rd, statistic = dg$statistic, p = dg$p.value,
    boot_max = max(dg$boot), axes = paste(dg$axes, collapse = " & "), AIC = stats::AIC(fits[[m]]),
    seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  if (rd == "allpairs") top[[m]] <- cbind(model = m, utils::head(dg$pairs[order(-dg$pairs$score), ], 5))
  saveRDS(dg, file.path(out, sprintf("dg290_%s_%s.rds", m, rd)))
  cat(sprintf("%-5s %-8s p = %.3f  S = %6.2f  boot max %6.2f  axes %s\n", m, rd, dg$p.value, dg$statistic,
              max(dg$boot), paste(dg$axes, collapse = " & ")))
  flush.console()
}
R <- do.call(rbind, res)
utils::write.csv(R, file.path(out, "dg290_models.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, top), file.path(out, "dg290_allpairs_top.csv"), row.names = FALSE)
print(do.call(rbind, top), row.names = FALSE)

cat("\n== M0, all-pairs, over 20 further tie-break seeds\n")
pv <- t(sapply(1:20, function(s) {
  set.seed(20269000L + s); dg <- deepgof1(fits$M0, B = 199, reading = "allpairs")
  c(seed = s, p = dg$p.value, S = dg$statistic, meanbp = as.numeric("meanbp" %in% dg$axes))
}))
print(summary(pv[, "p"])); cat("score range", round(range(pv[, "S"]), 2),
                                "| max pair includes meanbp in", sum(pv[, "meanbp"]), "of 20\n")
utils::write.csv(pv, file.path(out, "dg290_allpairs_tiebreak_M0.csv"), row.names = FALSE)
