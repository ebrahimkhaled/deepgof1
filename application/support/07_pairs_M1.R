# stage1_fixedaxes_SUPPORT.R -- the repaired model read on the SAME axes the M0 map flagged.
# After repair (b) the deepgof1() axis rule picks two spline columns of meanbp, so its map no longer
# looks at meanbp x num.co. Here the axes are fixed in advance to the covariates meanbp and num.co
# (taken from the data, not from the model matrix) and everything else is the shipped test: the same
# pretrained scorer, K = 6, B = 199 refits of the model under its own fit, the same tie-break rule.
# Fixing the axes before looking at M1 keeps the test valid. Run for M0 too, as the check that the
# fixed-axis reading reproduces the shipped result there.
suppressMessages({ library(ebrahim.gof); library(splines) })
source("data_prep.R")
k <- 2L
p <- prep(k)
DG_SEED <- 20261929L + k
out <- file.path(HERE, "results", "screen_SUPPORT")
M <- ebrahim.gof:::deepgof1_weights
score <- function(m) ebrahim.gof:::.dg_score((m - M$mu) / M$sd, M)
AX <- NULL

## the residual map of .dg_map with the two axes given as covariate columns of the data
map_fixed <- function(fit, dat, tb, K = 6L) {
  ph <- as.numeric(stats::fitted(fit)); r <- dat$y - ph; n <- length(r)
  rk <- function(v) { o <- order(v, tb); rr <- integer(n); rr[o] <- seq_len(n); rr }
  bn <- function(v) pmin(K, 1L + floor(K * (rk(v) - 1) / n))
  idx <- factor((bn(dat[[AX[1]]]) - 1L) * K + bn(dat[[AX[2]]]), levels = 1:(K * K))
  s <- tapply(r, idx, sum); vv <- tapply(ph * (1 - ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  as.numeric(s) / sqrt(pmax(as.numeric(vv), 1e-8))
}

fixed_test <- function(fit, B = 199L) {
  dat <- p$data
  set.seed(DG_SEED)
  tb <- sample.int(nrow(dat))                     # both axes have ties, as in deepgof1()
  m0 <- map_fixed(fit, dat, tb); S <- score(m0)
  ph <- as.numeric(stats::fitted(fit)); frm <- stats::formula(fit)
  Sb <- numeric(B)
  for (b in seq_len(B)) {
    dat$y <- stats::rbinom(length(ph), 1L, ph)
    fb <- tryCatch(suppressWarnings(stats::glm(frm, data = dat, family = stats::binomial())),
                   error = function(e) NULL)
    Sb[b] <- if (is.null(fb)) Inf else score(map_fixed(fb, dat, tb))
  }
  list(map = m0, S = S, p = (1 + sum(Sb >= S)) / (B + 1), boot = Sb)
}

## Where is the misfit that remains in M1 = ns(meanbp, 3)? M1 read on every pair of the eight continuous
## covariates, axes fixed in advance (a map per pair; this is a search, reported as such, not a test).
M1 <- readRDS(file.path(out, "M1_b.rds"))$fit
cv <- c("age", "meanbp", "hrt", "resp", "temp", "crea", "sod", "wblc")
res <- data.frame()
for (pr in utils::combn(cv, 2, simplify = FALSE)) {
  AX <- pr
  r <- fixed_test(M1, B = 99L)
  res <- rbind(res, data.frame(ax1 = pr[1], ax2 = pr[2], S = r$S, p = r$p, maxcell = max(abs(r$map))))
}
res <- res[order(-res$S), ]
utils::write.csv(res, file.path(out, "ladder_M1_pairs.csv"), row.names = FALSE)
print(format(res, digits = 3), row.names = FALSE)
