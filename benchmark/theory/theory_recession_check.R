## Run with the working directory theory/_support/rdeploy of this archive (the script reads
## deepgof_forward.R, weights.csv and scaler.csv from there).
## theory_recession_check.R -- the computable half of the consistency theorem, checked on the
## SHIPPED weights before a single line of proof is typed.
##
## CLAIM UNDER TEST. The scored statistic is S(z) with z = (m - mu)/sd and the CNN S built from
## affine maps, ReLU, max-pool and global max/mean -- a PIECEWISE AFFINE function. Then for any
## direction v the recession slope
##      Sinf(v) = lim_{t->inf} S(((t*v) - mu)/sd) / t
## exists, is finite, and S(t*v) is EXACTLY affine in t once t is large enough to fix the active
## linear piece. Under a fixed alternative the deployed map is m_n ~ sqrt(n)*mu_dir + O_p(1)
## (White-projection + CLT), the bootstrap scores stay O_p(1), so
##      Sinf(direction) > 0  <=>  power -> 1.
## If second differences of S along t*v do NOT vanish, the theorem's engine is wrong and no proof
## should be attempted.
source("deepgof_forward.R")
W  <- load_weights("weights.csv")
ms <- as.matrix(utils::read.csv("scaler.csv", header = FALSE))
mu <- ms[1, ]; sdv <- ms[2, ]
K  <- 6L
sc <- function(m) deepgof_score((m - mu) / sdv, W, K)

## map-vector index = (i1-1)*K + i2  (i1 = x1-bin -> rows after byrow=TRUE)
vec_rows <- function(val) as.numeric(sapply(1:K, function(i1) rep(val[i1], K)))
vec_cols <- function(val) rep(val, K)

stripe <- c(1.2, -0.4, -1.1, -1.1, -0.4, 1.2)     # omitted-quadratic curvature stripe
step   <- c(-1, -1, -1, 1, 1, 1)                  # threshold/step departure
dirs <- list(
  stripe_x1   = vec_rows(stripe),
  stripe_x2   = vec_cols(stripe),
  step_x1     = vec_rows(step),
  saddle      = as.numeric(outer(step, step)),    # omitted-interaction saddle
  diffuse_up  = rep(1, K*K),                      # uniform calibration drift (mean-pool bait)
  probit_tail = vec_rows(c(-1, .3, .7, .7, .3, -1) * .25),  # weak symmetric tail pattern
  random      = { set.seed(11); rnorm(K*K) }
)
TS <- c(5, 10, 20, 40, 80, 160, 320)
cat(sprintf("%-12s %10s %10s %12s %s\n", "direction", "S(160v)", "slope", "2nd-diff", "verdict"))
for (nm in names(dirs)) {
  v <- dirs[[nm]] / sqrt(sum(dirs[[nm]]^2))       # unit direction
  s  <- sapply(TS, function(t) sc(t * v))
  sl <- diff(s) / diff(TS)
  d2 <- abs(diff(sl))
  slope <- sl[length(sl)]
  lin   <- max(tail(d2, 2))
  cat(sprintf("%-12s %10.3f %10.4f %12.2e %s\n", nm, s[6], slope, lin,
      if (lin < 1e-8) sprintf("EXACTLY affine; Sinf %s 0 -> %s", 
                              ifelse(slope > 0, ">", "<="),
                              ifelse(slope > 0, "CONSISTENT", "blind/anti")) 
      else "NOT affine -- engine wrong"))
  ## and the mirror direction: blindness is direction-specific, record it
  s2 <- sapply(TS, function(t) sc(-t * v)); sl2 <- diff(s2)/diff(TS)
  cat(sprintf("%-12s %10s %10.4f  (mirror -v)\n", "", "", sl2[length(sl2)]))
}
