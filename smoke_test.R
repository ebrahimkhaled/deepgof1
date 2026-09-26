## smoke_test.R -- end-to-end check of the archived deepgof1(): one misspecified model
## (omitted quadratic, should reject), one correct model (should not, usually), and one model
## with a binary covariate (tied axis: ties are broken at random, once per call).
## Run from the archive root:  Rscript smoke_test.R
source("R/deepgof1.R")
set.seed(42); n <- 150
x1 <- runif(n,-3,3); x2 <- rnorm(n)
y  <- rbinom(n, 1, plogis(0.3 + 0.8*x1 - 0.5*x2 + 0.9*x1^2 - 0.9*mean(x1^2)))
d  <- data.frame(x1, x2, y)
fit <- glm(y ~ x1 + x2, data = d, family = binomial())
t0 <- Sys.time(); r <- deepgof1(fit); el <- as.numeric(Sys.time()-t0, units="secs")
print(r); cat(sprintf("runtime: %.1f s
", el))
y0 <- rbinom(n, 1, plogis(0.3 + 0.8*x1 - 0.5*x2))
d0 <- data.frame(x1, x2, y = y0)
fit0 <- glm(y ~ x1 + x2, data = d0, family = binomial())
r0 <- deepgof1(fit0)
cat(sprintf("correct model p-value: %.3f
", r0$p.value))
xb <- rbinom(n, 1, .5)
yb <- rbinom(n, 1, plogis(-0.5 + 0.9*x1 + 1.1*xb))
db <- data.frame(x1, xb, y = yb)
fitb <- glm(y ~ x1 + xb, data = db, family = binomial())
rb <- deepgof1(fitb)
cat(sprintf("binary-covariate model (random tie-break) p-value: %.3f
", rb$p.value))
