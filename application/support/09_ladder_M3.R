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

## Ladder rung 2 (rule fixed before fitting): natural splines, 3 df, on every continuous physiological
## covariate (meanbp, hrt, resp, temp, crea, sod, wblc); age, num.co and sex as in M0.
F2 <- y ~ ns(age, 3) + ns(meanbp, 3) + ns(hrt, 3) + ns(resp, 3) + ns(temp, 3) + ns(crea, 3) + ns(sod, 3) +
  ns(wblc, 3) + num.co + sex
M2 <- stats::glm(F2, data = p$data, family = stats::binomial())
saveRDS(M2, file.path(out, "M3.rds"))
cat("AIC  M0", round(stats::AIC(readRDS(file.path(out, "M0.rds"))), 1), " M1",
    round(stats::AIC(readRDS(file.path(out, "M1_b.rds"))$fit), 1), " M3", round(stats::AIC(M2), 1), "\n")
cv <- c("age", "meanbp", "hrt", "resp", "temp", "crea", "sod", "wblc")
res <- data.frame()
for (pr in utils::combn(cv, 2, simplify = FALSE)) {
  AX <- pr
  r <- fixed_test(M2, B = 99L)
  res <- rbind(res, data.frame(ax1 = pr[1], ax2 = pr[2], S = r$S, p = r$p, maxcell = max(abs(r$map))))
}
res <- res[order(-res$S), ]
utils::write.csv(res, file.path(out, "ladder_M3_pairs.csv"), row.names = FALSE)
cat("pairs with p <= .05:", sum(res$p <= .05), "of", nrow(res), "\n")
print(format(head(res, 10), digits = 3), row.names = FALSE)
pv <- function(z) { v <- z$p_value; if (length(v) == 1 && is.finite(v)) as.numeric(v) else NA_real_ }
ctx <- ebrahim.gof:::.gof_context(M2, G = 10); set.seed(20261933L)
print(round(c(HL = pv(ebrahim.gof:::gof_hl(ctx)), mHL = pv(ebrahim.gof:::gof_ftest(ctx)),
  PH = pv(ebrahim.gof:::gof_ph_test(ctx)), HLew = pv(ebrahim.gof:::gof_hlw(ctx)),
  OR = pv(ebrahim.gof:::gof_osius(ctx)), Stukel = pv(ebrahim.gof:::gof_stukel(ctx, list(form = "marginal"))),
  GiViTI = givitiR::givitiCalibrationTest(o = as.numeric(ctx$y), e = as.numeric(ctx$ph), devel = "internal")$p.value,
  StuteZhu = pv(ebrahim.gof:::gof_stutezhu(ctx))), 4))
