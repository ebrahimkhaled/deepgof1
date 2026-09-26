## 01_batched_net.R -- a BATCHED re-expression of the shipped forward pass, plus its verification.
##
## The shipped deepgof_score() costs ~30 ms per map in base R, which at M = 1e5 x 2 laws x 40
## datasets x 2 cells would be ~130 core-hours. The batched version below does the identical
## arithmetic on N maps at once: the activation of a whole chunk is held as a (C) x (H*W*N)
## matrix, each conv layer is nine O x C matrix multiplies plus nine shifted adds (the zero
## padding is just the omission of the out-of-range positions), and both pooling layers are
## vectorised with pmax over spatial strides instead of a loop over maps.
##
## It is verified against deepgof_score() map by map; if the agreement is not at machine
## precision the pilot does not run.

## geometry of one conv layer: for each of the 9 kernel offsets, which spatial positions of the
## input feed which spatial positions of the output. Built once per (H, W).
.conv_geom <- function(H, W) {
  g <- vector("list", 9L); k <- 0L
  for (di in -1:1) for (dj in -1:1) {
    k <- k + 1L
    hs <- seq_len(H); ws <- seq_len(W)
    hv <- hs[hs + di >= 1L & hs + di <= H]
    wv <- ws[ws + dj >= 1L & ws + dj <= W]
    dst <- as.vector(outer(hv, wv, function(a, b) (b - 1L) * H + a))
    src <- as.vector(outer(hv + di, wv + dj, function(a, b) (b - 1L) * H + a))
    g[[k]] <- list(di = di, dj = dj, dst = dst, src = src)
  }
  g
}

## x: C x (H*W*n) ; Wt: O x C x 3 x 3 ; b: O  -->  O x (H*W*n)
.conv3b <- function(x, Wt, b, H, W, n, geom) {
  O <- dim(Wt)[1]; HW <- H * W
  out <- matrix(b, O, HW * n)                       # bias, recycled down columns
  off <- rep((0:(n - 1L)) * HW, each = length(geom[[1]]$dst))
  for (k in seq_len(9L)) {
    gk <- geom[[k]]
    ## PyTorch kernel index: weight[, , di+2, dj+2] multiplies the input at (h+di, w+dj)
    Wk <- Wt[, , gk$di + 2L, gk$dj + 2L, drop = FALSE][, , 1, 1, drop = FALSE]
    dim(Wk) <- c(O, dim(Wt)[2])
    ln <- length(gk$dst)
    sc <- rep(gk$src, n) + rep((0:(n - 1L)) * HW, each = ln)
    dc <- rep(gk$dst, n) + rep((0:(n - 1L)) * HW, each = ln)
    out[, dc] <- out[, dc] + Wk %*% x[, sc, drop = FALSE]
  }
  out
}

deepgof_score_batch <- function(M36, W, K = 6L, chunk = 4000L) {
  N <- nrow(M36)
  perm <- as.vector(t(matrix(seq_len(K * K), K, K)))   # row-major (PyTorch) -> R column-major
  g6 <- .conv_geom(K, K); g3 <- .conv_geom(K %/% 2L, K %/% 2L)
  out <- numeric(N)
  for (st in seq(1L, N, by = chunk)) {
    id <- st:min(st + chunk - 1L, N)
    nn <- length(id)
    x <- matrix(t(M36[id, , drop = FALSE])[perm, , drop = FALSE], 1L, K * K * nn)
    x <- pmax(.conv3b(x, W[["conv.0.weight"]], W[["conv.0.bias"]], K, K, nn, g6), 0)
    x <- pmax(.conv3b(x, W[["conv.2.weight"]], W[["conv.2.bias"]], K, K, nn, g6), 0)
    ## 2x2 stride-2 max pool on the K x K grid -> (K/2) x (K/2)
    H2 <- K %/% 2L; C <- nrow(x)
    ## column of (h,w,q) is (q-1)*K*K + (w-1)*K + h ; output (i,j,q) takes h in {2i-1,2i}, w in {2j-1,2j}
    base <- as.vector(outer(seq_len(H2), seq_len(H2), function(i, j) (2*j - 2L) * K + (2*i - 1L)))
    base <- rep(base, nn) + rep((0:(nn - 1L)) * K * K, each = H2 * H2)
    o <- pmax(x[, base, drop = FALSE], x[, base + 1L, drop = FALSE],
              x[, base + K, drop = FALSE], x[, base + K + 1L, drop = FALSE])
    x <- pmax(.conv3b(o, W[["conv.5.weight"]], W[["conv.5.bias"]], H2, H2, nn, g3), 0)
    sp <- H2 * H2; C <- nrow(x)
    gmax <- x[, seq(1L, sp * nn, by = sp), drop = FALSE]
    gsum <- gmax
    for (q in 2:sp) {
      sl <- x[, seq(q, sp * nn, by = sp), drop = FALSE]
      gmax <- pmax(gmax, sl); gsum <- gsum + sl
    }
    e <- rbind(gmax, gsum / sp)                        # GlobalMax then GlobalMean -> 64 x nn
    h <- pmax(W[["head.0.weight"]] %*% e + W[["head.0.bias"]], 0)
    out[id] <- as.numeric(W[["head.2.weight"]] %*% h + W[["head.2.bias"]])
  }
  out
}

if (sys.nframe() == 0L) {
  RD <- "/path/to/project/paper_deepgof/code/rdeploy"
  source(file.path(RD, "deepgof_forward.R"))
  W  <- load_weights(file.path(RD, "weights.csv"))
  Xt <- as.matrix(utils::read.csv(file.path(RD, "test_input.csv"), header = FALSE))
  Sp <- as.numeric(readLines(file.path(RD, "test_output.csv")))
  cat(sprintf("batched vs test_output.csv : max|diff| = %.3e\n",
              max(abs(deepgof_score_batch(Xt, W) - Sp))))
  set.seed(11)
  R <- 600L
  MM <- matrix(rnorm(R * 36, 0, 1.5), R, 36)
  MM[sample.int(R * 36, R * 6)] <- 0
  s1 <- apply(MM, 1, function(r) deepgof_score(r, W))
  s2 <- deepgof_score_batch(MM, W)
  cat(sprintf("batched vs deepgof_score() : %d maps, max|diff| = %.3e\n", R, max(abs(s1 - s2))))
  BIG <- matrix(rnorm(40000 * 36), 40000, 36)
  t0 <- Sys.time(); invisible(deepgof_score_batch(BIG, W))
  cat(sprintf("batched throughput         : %.4f ms / map\n",
              1000 * as.numeric(difftime(Sys.time(), t0, units = "secs")) / 40000))
  n <- 50L; mm <- cbind(1, rnorm(n), runif(n, -3, 3)); colnames(mm) <- c("(Intercept)","x1","x2")
  bfam <- binomial()
  ys <- lapply(1:3000, function(i) rbinom(n, 1, 0.5))
  t2 <- Sys.time(); for (i in 1:3000) glm.fit(mm, ys[[i]], family = bfam)
  cat(sprintf("glm.fit n=50 (family cached): %.3f ms\n",
              1000 * as.numeric(difftime(Sys.time(), t2, units = "secs")) / 3000))
  n <- 200L; mm <- cbind(1, rnorm(n), runif(n, -3, 3)); colnames(mm) <- c("(Intercept)","x1","x2")
  ys <- lapply(1:3000, function(i) rbinom(n, 1, 0.5))
  t3 <- Sys.time(); for (i in 1:3000) glm.fit(mm, ys[[i]], family = bfam)
  cat(sprintf("glm.fit n=200(family cached): %.3f ms\n",
              1000 * as.numeric(difftime(Sys.time(), t3, units = "secs")) / 3000))
}
