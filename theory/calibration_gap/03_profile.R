## 03_profile.R -- where the per-score time goes, single process, at the real chunk size.
suppressMessages(library(brglm2))
RD  <- "/path/to/project/paper_deepgof/code/rdeploy"
OUT <- "/path/to/project/paper_deepgof/theory/calibration_gap"
source(file.path(RD, "deepgof_forward.R")); source(file.path(OUT, "01_batched_net.R"))
W  <- load_weights(file.path(RD, "weights.csv"))
ms <- as.matrix(utils::read.csv(file.path(RD, "scaler.csv"), header = FALSE))
MU <- ms[1, ]; SDV <- ms[2, ]; BFAM <- binomial()

gen <- function(n) {
  Z <- matrix(rnorm(n*2), n, 2); Z[,2] <- runif(n,-3,3)
  d <- as.data.frame(Z); names(d) <- c("x1","x2")
  MM <- model.matrix(~., d); bet <- c(0, runif(2,-.8,.8))
  lin <- as.numeric(MM %*% bet); lin <- lin/max(sd(lin),1e-8)*1.2
  b0 <- uniroot(function(b) mean(plogis(b+lin))-0.5, c(-25,25))$root
  list(mm = MM, p0 = plogis(b0+lin))
}
mkA <- function(mm, K = 6L) {
  n <- nrow(mm)
  bn <- function(v) pmin(K, 1L + floor(K*(rank(v, ties.method="first")-1)/n))
  cell <- (bn(mm[,"x1"])-1L)*K + bn(mm[,"x2"])
  A <- matrix(0, 36, n); A[cbind(cell, seq_len(n))] <- 1; A
}
for (n in c(50L, 200L)) {
  g <- gen(n); A <- mkA(g$mm); R <- 20000L
  ys <- lapply(seq_len(R), function(i) rbinom(n, 1L, g$p0))
  t <- Sys.time(); for (i in seq_len(R)) rbinom(n, 1L, g$p0)
  e_rb <- as.numeric(difftime(Sys.time(), t, units="secs"))/R
  t <- Sys.time(); fits <- lapply(ys, function(y) glm.fit(g$mm, y, family = BFAM))
  e_fit <- as.numeric(difftime(Sys.time(), t, units="secs"))/R
  t <- Sys.time()
  for (i in seq_len(R)) { f <- tryCatch(suppressWarnings(glm.fit(g$mm, ys[[i]], family=BFAM)), error=function(e) NULL) }
  e_fitw <- as.numeric(difftime(Sys.time(), t, units="secs"))/R
  phs <- lapply(fits, function(f) as.numeric(f$fitted.values))
  t <- Sys.time()
  mp <- matrix(0, R, 36)
  for (i in seq_len(R)) { ph <- phs[[i]]
    mp[i,] <- as.numeric(A %*% (ys[[i]]-ph))/sqrt(pmax(as.numeric(A %*% (ph*(1-ph))),1e-8)) }
  e_map <- as.numeric(difftime(Sys.time(), t, units="secs"))/R
  for (ch in c(2500L, 5000L)) {
    t <- Sys.time()
    invisible(deepgof_score_batch((mp - rep(MU, each=R))/rep(SDV, each=R), W, chunk=ch))
    cat(sprintf("n=%3d net chunk=%d : %.4f ms/map\n", n, ch,
                1000*as.numeric(difftime(Sys.time(), t, units="secs"))/R))
  }
  cat(sprintf("n=%3d rbinom %.4f | glm.fit %.4f | glm.fit+tryCatch %.4f | map %.4f  (ms)\n\n",
              n, 1000*e_rb, 1000*e_fit, 1000*e_fitw, 1000*e_map))
}
