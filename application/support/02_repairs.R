# stage1_repair_SUPPORT.R -- PLAN.md Sec. 4 step 2: the fixed repair list, in order, DeepGOF-1 only.
# The map of M0 flags meanbp (rows) most: its first and last sixths are positive in every column.
#   (a) M0 + meanbp:num.co     (b) M0 with ns(meanbp, 3)     (c) ns(meanbp, 3) + ns(num.co, 3)
#   (d) (c) + meanbp:num.co
# The chosen repair is the first whose DeepGOF-1 p-value exceeds .05. Same seed as the screen.
suppressMessages({ library(ebrahim.gof); library(splines) })
source("data_prep.R")
k <- 2L
p <- prep(k)
DG_SEED <- 20261929L + k
rest <- "age + hrt + resp + temp + crea + sod + wblc + sex"
forms <- list(
  a = paste("y ~ meanbp + num.co + meanbp:num.co +", rest),
  b = paste("y ~ ns(meanbp, 3) + num.co +", rest),
  c = paste("y ~ ns(meanbp, 3) + ns(num.co, 3) +", rest),
  d = paste("y ~ ns(meanbp, 3) + ns(num.co, 3) + meanbp:num.co +", rest))
out <- file.path(HERE, "results", "screen_SUPPORT")
res <- data.frame()
for (r in names(forms)) {
  fit <- stats::glm(stats::as.formula(forms[[r]]), data = p$data, family = stats::binomial())
  set.seed(DG_SEED)
  dg <- deepgof1(fit, B = 199)
  res <- rbind(res, data.frame(repair = r, formula = forms[[r]], deepgof1_p = dg$p.value,
                               statistic = dg$statistic, axes = paste(dg$axes, collapse = " & "),
                               AIC = stats::AIC(fit), deviance = stats::deviance(fit)))
  saveRDS(list(fit = fit, dg = dg), file.path(out, paste0("M1_", r, ".rds")))
  cat(r, ": p =", dg$p.value, " S =", round(dg$statistic, 2), " axes =", paste(dg$axes, collapse = " & "),
      " AIC =", round(stats::AIC(fit), 1), "\n")
}
utils::write.csv(res, file.path(out, "stage1_repairs.csv"), row.names = FALSE)
cat("M0 AIC =", round(stats::AIC(readRDS(file.path(out, "M0.rds"))), 1), "\n")
cat("chosen repair:", if (any(res$deepgof1_p > .05)) res$repair[which(res$deepgof1_p > .05)[1]] else "none", "\n")
