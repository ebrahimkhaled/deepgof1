## deepgof_forward.R -- the shipped test's forward pass, in BASE R. No PyTorch, no reticulate.
##
## WHY THIS FILE EXISTS. PyTorch is needed to TRAIN the network. It is not needed to RUN it. The
## shipped object is 18,273 frozen numbers and a forward pass that is a handful of small matrix
## multiplies -- so the test can live inside `ebrahim.gof` as pure R, with no Python dependency at
## the analyst's end. A goodness-of-fit test that required a PyTorch install would not be adopted by
## the people who run goodness-of-fit tests.
##
## The architecture (see code/map_train.py, K = 6):
##   1x6x6 -> Conv3x3(16, pad 1) -> ReLU -> Conv3x3(32, pad 1) -> ReLU -> MaxPool2
##         -> Conv3x3(32, pad 1) -> ReLU -> [GlobalMax(32) , GlobalMean(32)] = 64
##         -> Linear(64,64) -> ReLU -> Linear(64,1) -> scalar score
##
## Convolution is done by im2col: each 3x3 patch becomes a column, so a conv layer is ONE matrix
## multiply. That keeps base R fast enough for the 200 forward passes a deployed test needs
## (observed row + 199 bootstrap replicates).

## ---- load the exported weights -----------------------------------------------------------------
load_weights <- function(path) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  W <- list()
  for (i in seq_len(nrow(d))) {
    shp <- as.integer(strsplit(d$shape[i], "x")[[1]])
    v   <- as.numeric(strsplit(trimws(d$values[i]), " +")[[1]])
    ## PyTorch stores row-major; R arrays are column-major, so fill a reversed-dim array and
    ## permute back. Getting this wrong is the classic silent error in a port like this.
    W[[d$name[i]]] <- if (length(shp) == 1L) v else
      aperm(array(v, dim = rev(shp)), perm = length(shp):1)
  }
  W
}

## ---- im2col: (C,H,W) -> ((C*9) x (H*W)) with zero padding of 1 ---------------------------------
im2col3 <- function(x) {                      # x: C x H x W
  C <- dim(x)[1]; H <- dim(x)[2]; Wd <- dim(x)[3]
  xp <- array(0, c(C, H + 2L, Wd + 2L))
  xp[, 2:(H + 1L), 2:(Wd + 1L)] <- x
  out <- matrix(0, C * 9L, H * Wd)
  k <- 0L
  for (di in 0:2) for (dj in 0:2) {           # kernel row, kernel col (PyTorch order)
    k <- k + 1L
    patch <- xp[, (1L + di):(H + di), (1L + dj):(Wd + dj), drop = FALSE]
    out[((k - 1L) * C + 1L):(k * C), ] <- matrix(patch, C, H * Wd)
  }
  out
}

conv3 <- function(x, Wt, b) {                 # Wt: OUT x IN x 3 x 3 ; b: OUT
  C <- dim(x)[1]; H <- dim(x)[2]; Wd <- dim(x)[3]; O <- dim(Wt)[1]
  ## flatten the kernel to O x (C*9) in the SAME order im2col3 stacks its rows
  Wm <- matrix(0, O, C * 9L)
  k <- 0L
  for (di in 1:3) for (dj in 1:3) {
    k <- k + 1L
    Wm[, ((k - 1L) * C + 1L):(k * C)] <- Wt[, , di, dj, drop = FALSE][, , 1, 1]
  }
  z <- Wm %*% im2col3(x) + b                  # O x (H*W), b recycles down columns
  array(z, c(O, H, Wd))
}

relu    <- function(x) pmax(x, 0)
maxpool2 <- function(x) {                     # stride-2 2x2 max, C x H x W with H,W even
  C <- dim(x)[1]; H <- dim(x)[2]; Wd <- dim(x)[3]
  o <- array(0, c(C, H %/% 2L, Wd %/% 2L))
  for (i in seq_len(H %/% 2L)) for (j in seq_len(Wd %/% 2L))
    o[, i, j] <- apply(x[, (2*i-1):(2*i), (2*j-1):(2*j), drop = FALSE], 1, max)
  o
}

## ---- the score: one map (length 36) -> one scalar ----------------------------------------------
deepgof_score <- function(map_vec, W, K = 6L) {
  ## THE PORT TRAP: PyTorch's x.view(-1,1,K,K) fills ROW-major, R's array() fills COLUMN-major,
  ## so a naive array(v, c(1,K,K)) silently TRANSPOSES the map. byrow=TRUE restores PyTorch's
  ## ordering, which is also the order the R map builder emits (cell index = (i1-1)*K + i2).
  x <- array(matrix(as.numeric(map_vec), K, K, byrow = TRUE), c(1L, K, K))
  x <- relu(conv3(x, W[["conv.0.weight"]], W[["conv.0.bias"]]))
  x <- relu(conv3(x, W[["conv.2.weight"]], W[["conv.2.bias"]]))
  x <- maxpool2(x)
  x <- relu(conv3(x, W[["conv.5.weight"]], W[["conv.5.bias"]]))
  e <- c(apply(x, 1, max), apply(x, 1, mean))          # GlobalMax  then GlobalMean -> 64
  h <- relu(as.vector(W[["head.0.weight"]] %*% e + W[["head.0.bias"]]))
  as.numeric(W[["head.2.weight"]] %*% h + W[["head.2.bias"]])
}

deepgof_scores <- function(M, W, K = 6L)      # M: rows are maps
  apply(M, 1, function(r) deepgof_score(r, W, K))

## ---- verification against PyTorch --------------------------------------------------------------
if (sys.nframe() == 0L) {
  here <- "/path/to/project/paper_deepgof/code/rdeploy"
  W  <- load_weights(file.path(here, "weights.csv"))
  X  <- as.matrix(utils::read.csv(file.path(here, "test_input.csv"),  header = FALSE))
  Sp <- as.numeric(readLines(file.path(here, "test_output.csv")))
  Sr <- deepgof_scores(X, W)
  cat("PyTorch :", sprintf("%.6f", Sp), "\n")
  cat("base R  :", sprintf("%.6f", Sr), "\n")
  cat(sprintf("max |difference| = %.3e\n", max(abs(Sr - Sp))))
  t0 <- Sys.time(); for (i in 1:200) deepgof_score(X[1, ], W)
  cat(sprintf("200 forward passes (one deployed test): %.3f s\n",
              as.numeric(Sys.time() - t0, units = "secs")))
}
