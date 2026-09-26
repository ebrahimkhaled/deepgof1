# test_discriminates.R -- would the new unit test have caught the old (row-order) statistic?
source("/path/to/project/paper_deepgof/theory/tie_fix/common.R")
tied_sorted_data <- function(seed, n = 200) {
  set.seed(seed)
  x1 <- sample(1:4, n, replace = TRUE); x2 <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + 0.4 * x1 + 0.8 * x2))
  d <- data.frame(y, x1, x2); d[order(d$y), ]
}
d <- tied_sorted_data(23)
fit <- glm(y ~ x1 + x2, data = d, family = binomial())
for (v in c("old", "new")) {
  f <- if (v == "old") dg_old else dg_new
  p <- vapply(1:8, function(s) { set.seed(100 + s); f(fit, B = 19)$p.value }, 0)
  ds <- d[sample.int(nrow(d)), ]
  fs <- glm(y ~ x1 + x2, data = ds, family = binomial())
  s1 <- vapply(1:8, function(s) { set.seed(200 + s); f(fit, B = 1)$statistic }, 0)
  s2 <- vapply(1:8, function(s) { set.seed(200 + s); f(fs, B = 1)$statistic }, 0)
  cat(v, ": p =", round(p, 3), " mean", round(mean(p), 3), "| S sorted", round(s1, 2),
      "| S shuffled", round(s2, 2), "| Wilcoxon p", signif(suppressWarnings(wilcox.test(s1, s2))$p.value, 3), "\n")
}
