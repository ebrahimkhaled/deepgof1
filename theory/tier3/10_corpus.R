# 10_corpus.R -- Tier 3 item B: one full-size training corpus for a held-out-family network.
# The deposited generator (release/deepgof1_v2_final/training/01_corpus.R) is read and its functions are
# evaluated VERBATIM, except the departure prior's KINDS/PR lines (and, for NoQuadStrict, the term draw):
#   Control       all nine kinds (the deposited prior, unchanged)
#   NoInter       without "inter" and "indexq"                    (author's decision)
#   NoQuad        without "quad"  and "indexq"                    (author's decision)
#   NoQuadStrict  NoQuad, and an "inter" term is redrawn when its two columns coincide. In the deposited
#                 generator xj = X[, j %% ncol(X) + 1] with j != i, which equals xi whenever
#                 j %% ncol + 1 == i: always on the 2-covariate design, 1 time in 3 on the 4-covariate
#                 design. Such an "inter" term is xi^2, a quadratic, so NoQuad still trains on squares.
#                 NoQuadStrict removes those too (added by the agent; see REPORT.md, deviations).
# Sizes as deposited: 100 families x 42 alternatives + 4,200 bootstrap-replicate nulls, n = 200, K = 6.
# Seeds: SEED0 = 80101 + 1000 * s (s = 1, 2, 3); the same SEED0 for every variant at seed s (paired).
# Usage: Rscript 10_corpus.R <variant> <s>      Output: nets/<variant>_s<s>/corpus.csv (+ corpus_info.csv)
a <- commandArgs(TRUE); VAR <- a[1]; S <- as.integer(a[2])
stopifnot(VAR %in% c("Control", "NoInter", "NoQuad", "NoQuadStrict"), S %in% 1:3)
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3"
SRC  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/release/deepgof1_v2_final/training/01_corpus.R"
OUTD <- file.path(HERE, "nets", sprintf("%s_s%d", VAR, S)); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
OUT  <- file.path(OUTD, "corpus.csv")
if (file.exists(OUT)) { cat("exists:", OUT, "\n"); quit(save = "no") }
EXCL <- switch(VAR, Control = character(0), NoInter = c("inter", "indexq"), NoQuad = c("quad", "indexq"),
               NoQuadStrict = c("quad", "indexq"))

src <- readLines(SRC)
ki  <- grep("^  KINDS <- c[(]", src); pi_ <- grep("^  PR    <- c[(]", src); stopifnot(length(ki) == 1, length(pi_) == 1)
src[ki]  <- sprintf('  KINDS <- setdiff(c("quad","cubic","inter","sine","bump","step","abs","indexq","indexs"), c(%s))',
                    paste0('"', EXCL, '"', collapse = ","))
if (!length(EXCL)) src[ki] <- '  KINDS <- c("quad","cubic","inter","sine","bump","step","abs","indexq","indexs")'
src[pi_] <- '  PR    <- c(quad=.15,cubic=.10,inter=.15,sine=.12,bump=.12,step=.10,abs=.06,indexq=.10,indexs=.10)[KINDS]'
if (VAR == "NoQuadStrict") {
  kl <- grep("^    kind <- sample[(]KINDS, 1, prob = PR[)]$", src); il <- grep("^    i <- sample[(]seq_len[(]vars[)], 1[)]; j <- ", src)
  stopifnot(length(kl) == 1, length(il) == 1, il == kl + 1)
  src[kl] <- "    repeat { kind <- sample(KINDS, 1, prob = PR)"
  src[il] <- paste0(src[il], '; if (!(kind == "inter" && (j %% vars + 1L) == i)) break }')
}
grab <- function(name) { a <- grep(paste0("^", name, " <- function"), src); b <- a - 1 + grep("^}", src[a:length(src)])[1]
  eval(parse(text = src[a:b]), envir = globalenv()) }
grab("map_of"); grab("make_family"); grab("gen_random")

N_FAM <- 100L; PER_FAM <- 42L; N_NULL <- 4200L; NOBS <- 200L; K <- 6L; SEED0 <- 80101L + 1000L * S
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
## the prior actually used: every term of the 100 families, with its effective form
TT <- do.call(rbind, lapply(seq_len(N_FAM), function(f) { fam <- make_family(SEED0 + f); v <- if (fam$four) 4L else 2L
  do.call(rbind, lapply(fam$terms, function(t) data.frame(f = f, vars = v, kind = t$kind, i = t$i, jeff = t$j %% v + 1L))) }))
TT$self <- TT$kind == "inter" & TT$i == TT$jeff
stopifnot(!any(TT$kind %in% EXCL))
if (VAR == "NoQuadStrict") stopifnot(!any(TT$self))

t0 <- proc.time()[["elapsed"]]
rows <- vector("list", N_FAM * PER_FAM + N_NULL); ri <- 0L
for (f in seq_len(N_FAM)) { fam <- make_family(SEED0 + f)
  for (i in seq_len(PER_FAM)) { RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(SEED0 + f * 10000L + i)
    g <- gen_random(fam, NOBS)
    if (is.null(g$fit) || !g$fit$converged) next
    m <- tryCatch(map_of(g$fit, K), error = function(e) NULL); if (is.null(m)) next
    ri <- ri + 1L; rows[[ri]] <- c(1, m) } }
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(SEED0 * 7L); b <- 0L
while (b < N_NULL) {                                   # verbatim null loop of 01_corpus.R
  x1 <- runif(NOBS, -3, 3); x2 <- rnorm(NOBS); eta0 <- 0.3 + 0.8 * x1 - 0.5 * x2
  d0 <- data.frame(x1, x2, y = rbinom(NOBS, 1, plogis(eta0)))
  f0 <- suppressWarnings(glm(y ~ x1 + x2, data = d0, family = binomial())); if (!f0$converged) next
  d0$y <- rbinom(NOBS, 1, fitted(f0))
  fb <- tryCatch(suppressWarnings(glm(y ~ x1 + x2, data = d0, family = binomial())), error = function(e) NULL)
  if (is.null(fb) || !fb$converged) next
  m <- tryCatch(map_of(fb, K), error = function(e) NULL); if (is.null(m)) next
  b <- b + 1L; ri <- ri + 1L; rows[[ri]] <- c(0, m) }
M <- do.call(rbind, rows[seq_len(ri)]); colnames(M) <- c("y", paste0("c", seq_len(K * K)))
utils::write.csv(M, paste0(OUT, ".part"), row.names = FALSE); file.rename(paste0(OUT, ".part"), OUT)
info <- data.frame(variant = VAR, s = S, seed0 = SEED0, rows = nrow(M), alt = sum(M[, 1] == 1), null = sum(M[, 1] == 0),
  terms = nrow(TT), inter_terms = sum(TT$kind == "inter"), inter_self_square = sum(TT$self),
  quad_terms = sum(TT$kind == "quad"), indexq_terms = sum(TT$kind == "indexq"),
  secs = proc.time()[["elapsed"]] - t0)
utils::write.csv(info, file.path(OUTD, "corpus_info.csv"), row.names = FALSE)
print(info)
