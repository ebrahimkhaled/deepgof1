# 01_M0.R -- fit M0 and run the shipped DeepGOF-1 on it (axis rule, B = 199, fixed seed).
suppressMessages(library(ebrahim.gof))
source("data_prep.R")
p <- prep(2L)
fit <- stats::glm(p$formula, data = p$data, family = stats::binomial())
dir.create("results/screen_SUPPORT", showWarnings = FALSE, recursive = TRUE)
saveRDS(fit, "results/screen_SUPPORT/M0.rds")
set.seed(20261929L + 2L)
dg <- deepgof1(fit, B = 199)
cat("M0: n =", nrow(p$data), " deaths =", sum(p$data$y), " DeepGOF-1 p =", dg$p.value,
    " score =", round(dg$statistic, 2), " axes =", paste(dg$axes, collapse = " & "), "\n")
