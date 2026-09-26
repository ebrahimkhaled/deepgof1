## Run with the working directory theory/_support/rdeploy of this archive (the script reads
## deepgof_forward.R, weights.csv and scaler.csv from there).
## theory_sphere_check.R -- Condition C of the consistency theorem, checked on the shipped weights.
##
## CLOSED FORM. For a ReLU/max-pool/linear network S, the recession slope along v is the ZERO-BIAS
## network evaluated at v:  Sinf(v) = S_nobias(v)  (a.e.; structural induction over layers -- each
## relu(a t + b)/t -> relu(a), which is the biasless unit). The deployed statistic is
## F(m) = S((m - mu)/sd), so its recession along a map direction u is  S_nobias(u / sd).
##
## CONDITION C (what the omnibus consistency theorem needs):
##      inf over unit map directions u of S_nobias(u/sd)  >  0.
## Then ANY fixed alternative with a nonzero limiting map direction mu_dir is detected: the
## observed score grows like sqrt(n) * Sinf(mu_dir) while the B=199 bootstrap scores stay O_p(1).
source("deepgof_forward.R")
W  <- load_weights("weights.csv")
ms <- as.matrix(utils::read.csv("scaler.csv", header = FALSE))
mu <- ms[1, ]; sdv <- ms[2, ]
K  <- 6L
W0 <- W
for (nm in names(W0)) if (grepl("bias", nm)) W0[[nm]] <- W0[[nm]] * 0
Sinf <- function(u) deepgof_score(u / sdv, W0, K)          # closed-form recession slope
Sfull<- function(m) deepgof_score((m - mu) / sdv, W, K)

## 1. closed form vs measured large-t slope, on yesterday's directions
vec_rows <- function(val) as.numeric(sapply(1:K, function(i) rep(val[i], K)))
stripe <- c(1.2,-0.4,-1.1,-1.1,-0.4,1.2); step <- c(-1,-1,-1,1,1,1)
dirs <- list(stripe_x1=vec_rows(stripe), step_x1=vec_rows(step),
             saddle=as.numeric(outer(step,step)), diffuse=rep(1,K*K),
             random={set.seed(11); rnorm(K*K)})
cat("closed-form check (Sinf vs slope measured between t=2000 and t=4000):\n")
for (nm in names(dirs)) {
  v  <- dirs[[nm]]/sqrt(sum(dirs[[nm]]^2))
  sl <- (Sfull(4000*v) - Sfull(2000*v)) / 2000
  cat(sprintf("  %-10s Sinf = %8.4f   slope = %8.4f   |diff| = %.2e\n",
              nm, Sinf(v), sl, abs(Sinf(v)-sl)))
}

## 2. Condition C over the sphere: dense random sampling + worst-case refinement
set.seed(2026)
N <- 200000L
w <- matrix(rnorm(N*K*K), N)
w <- w / sqrt(rowSums(w^2))
vals <- numeric(N)
for (i in seq_len(N)) vals[i] <- Sinf(w[i, ])
cat(sprintf("\nCondition C over %d random unit directions:\n", N))
cat(sprintf("  min = %.4f | 0.1%% quantile = %.4f | median = %.4f | max = %.4f\n",
            min(vals), quantile(vals,.001), median(vals), max(vals)))
cat(sprintf("  directions with Sinf <= 0:  %d\n", sum(vals <= 0)))
## refine around the worst direction: local random search on the sphere
v <- w[which.min(vals), ]; best <- min(vals)
for (r in c(.3,.1,.03,.01)) for (it in 1:2000) {
  cand <- v + r*rnorm(K*K); cand <- cand/sqrt(sum(cand^2))
  s <- Sinf(cand); if (s < best) { best <- s; v <- cand }
}
cat(sprintf("  after local refinement around the worst direction: min = %.4f\n", best))
cat(sprintf("\nVERDICT: Condition C %s on the shipped DeepGOF-1 weights.\n",
            if (best > 0) "HOLDS  (inf > 0 -- omnibus consistency)" else "FAILS -- theorem must be conditional"))
