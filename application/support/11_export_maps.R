# export_maps_allpairs_SUPPORT.R -- the application figure over two CONTINUOUS covariates (referee M3).
# The all-pairs reading of M0 reaches its largest score on mean arterial pressure x respiratory rate. This
# script exports that map, the map of M1 on the same pair (the pair chosen from M0, before M1 was fitted),
# the clinical range of each rank sixth, and a simultaneous 95% threshold for the largest |cell| of the M0
# map on this pair (999 parametric-bootstrap maps, refitted on the design matrix, same tie-break).
suppressMessages({ library(ebrahim.gof); library(splines) })
stopifnot(utils::packageVersion("ebrahim.gof") >= "2.9.0")
source("data_prep.R")
p <- prep(2L); d <- p$data; out <- file.path(HERE, "results", "screen_SUPPORT")
DG_SEED <- 20261929L + 2L; K <- 6L; AX <- c("meanbp", "resp")
M0 <- readRDS(file.path(out, "M0.rds")); M1 <- readRDS(file.path(out, "M1_b.rds"))$fit
M2 <- readRDS(file.path(out, "M2.rds"))

ap0 <- readRDS(file.path(out, "dg290_M0_allpairs.rds"))
stopifnot(identical(unname(ap0$axes), AX))
set.seed(DG_SEED); f0 <- deepgof1(M0, B = 199, reading = "allpairs", covariates = AX)
set.seed(DG_SEED); f1 <- deepgof1(M1, B = 199, reading = "allpairs", covariates = AX)
set.seed(DG_SEED); f2 <- deepgof1(M2, B = 199, reading = "allpairs", covariates = AX)
stopifnot(isTRUE(all.equal(f0$map, ap0$map)))           # the same map as the all-pairs maximum
cat(sprintf("M0 on %s x %s: S = %.2f, p = %.3f\nM1 on the same pair: S = %.2f, p = %.3f\n",
            AX[1], AX[2], f0$statistic, f0$p.value, f1$statistic, f1$p.value))

## the tie-break deepgof1() drew (the first draw after the seed), for the bins and the threshold
set.seed(DG_SEED); tb <- sample.int(nrow(d)); n <- nrow(d)
bin <- function(v) { o <- order(v, tb); rr <- integer(n); rr[o] <- seq_len(n); pmin(K, 1L + floor(K * (rr - 1) / n)) }
rng <- function(v) sapply(split(v, bin(v)), function(z) paste(range(z), collapse = "|"))

## simultaneous threshold: 95% quantile of max|cell| over bootstrap maps of M0 on this pair
X  <- ebrahim.gof:::.dg_covariates(M0)[, AX]
cl <- ebrahim.gof:::.dg_cells(X, K, tb)
mm <- stats::model.matrix(M0); ph <- as.numeric(stats::fitted(M0))
set.seed(20261951L)
mx <- vapply(1:999, function(b) {
  ys <- stats::rbinom(n, 1L, ph)
  fb <- suppressWarnings(stats::glm.fit(mm, ys, family = stats::binomial()))
  max(abs(ebrahim.gof:::.dg_cellmap(cl$idx[[1]], ys - fb$fitted.values, fb$fitted.values, K)))
}, 0)
thr <- unname(stats::quantile(mx, .95))
obs <- as.vector(t(f0$map))
cat(sprintf("simultaneous 95%% threshold for max|cell|: %.2f; M0 cells beyond it: %d (largest %.2f)\n",
            thr, sum(abs(obs) > thr), max(abs(obs))))
cat("M1 cells beyond +-2:", sum(abs(f1$map) > 2), "  M0 cells beyond +-2:", sum(abs(f0$map) > 2), "\n")
print(round(f0$map, 2)); print(round(f1$map, 2))

maps <- data.frame(model = rep(c("M0", "M1", "M2"), each = 36), cell = rep(1:36, 3),
                   meanbp_bin = rep(rep(1:6, each = 6), 3), resp_bin = rep(rep(1:6, 6), 3),
                   value = c(as.vector(t(f0$map)), as.vector(t(f1$map)), as.vector(t(f2$map))))
utils::write.csv(maps, file.path(out, "maps_allpairs.csv"), row.names = FALSE)
utils::write.csv(data.frame(bin = 1:6, meanbp = rng(d$meanbp), resp = rng(d$resp)),
                 file.path(out, "maps_allpairs_bins.csv"), row.names = FALSE)
cat(sprintf("M2 on the same pair: S = %.2f, p = %.3f; cells beyond +-2: %d\n",
            f2$statistic, f2$p.value, sum(abs(f2$map) > 2)))
print(round(f2$map, 2))
utils::write.csv(data.frame(S0 = f0$statistic, p0 = f0$p.value, S1 = f1$statistic, p1 = f1$p.value,
                            S2 = f2$statistic, p2 = f2$p.value,
                            threshold = thr, beyond = sum(abs(obs) > thr), largest = max(abs(obs))),
                 file.path(out, "maps_allpairs_stats.csv"), row.names = FALSE)
print(utils::read.csv(file.path(out, "maps_allpairs_bins.csv")))
