# run_ushape.R -- referee M4: how often does the axis rule pick a covariate whose effect is a pure U-shape,
# and what does that cost? Five independent N(0,1) covariates; x1 has no linear effect, only a U-shape:
#   eta = -0.5 + 0.8 x2 + 0.6 x3 + 0.4 x4 + 0.3 x5 + C (x1^2 - 1),  C in {0.3, 0.6};  C = 0 is the null.
# The fitted model is linear in all five. n in {200, 500}; 1,000 nulls per n, 500 alternatives per cell.
# Arms on every dataset: deepgof1() with the axis rule (v1, the default) and with the all-pairs reading,
# Hosmer-Lemeshow and Stukel (run.all.gof), and the projection test (projection.gof), all with B = 199.
# Recorded: whether each DeepGOF-1 reading's map contains x1. ebrahim.gof 2.9.0 from the private library.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/ushape"
LIB  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290c"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
OUT <- file.path(HERE, if (SMOKE) "rows_smoke" else "rows"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
B <- 199L

SC <- rbind(data.frame(C = 0, n = c(200L, 500L), reps = 1000L),
            expand.grid(C = c(0.3, 0.6), n = c(200L, 500L), reps = 500L))
SC$id <- seq_len(nrow(SC))
if (SMOKE) SC$reps <- 2L
jobs <- list()
for (k in seq_len(nrow(SC))) for (ch in split(seq_len(SC$reps[k]), ceiling(seq_len(SC$reps[k]) / 25))) {
  f <- file.path(OUT, sprintf("sc%02d_r%04d.csv", SC$id[k], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(sc = SC[k, ], reps = ch, f = f)
}

one_chunk <- function(j) {
  suppressMessages(library(ebrahim.gof, lib.loc = LIB))
  stopifnot(utils::packageVersion("ebrahim.gof", lib.loc = LIB) == "2.9.0")
  R <- do.call(rbind, lapply(j$reps, function(r) {
    s <- 5000000L + j$sc$id * 100000L + r
    set.seed(s)
    n <- j$sc$n; X <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta <- -0.5 + X %*% c(0, 0.8, 0.6, 0.4, 0.3) + j$sc$C * (X[, 1]^2 - 1)
    d <- data.frame(X, y = stats::rbinom(n, 1, stats::plogis(eta)))
    fit <- suppressWarnings(stats::glm(y ~ x1 + x2 + x3 + x4 + x5, data = d, family = stats::binomial()))
    set.seed(s + 1L); a <- deepgof1(fit, B = B)
    set.seed(s + 2L); b <- deepgof1(fit, B = B, reading = "allpairs")
    set.seed(s + 3L)
    rc <- tryCatch(suppressWarnings(suppressMessages(run.all.gof(fit, tests = c("HL", "Stukel"),
                     include_slow = FALSE, install = "no"))), error = function(e) NULL)
    pv <- if (is.null(rc)) c(HL = NA, Stukel = NA) else stats::setNames(rc$p_value, rc$Test)
    set.seed(s + 4L)
    pj <- tryCatch(projection.gof(fit, B = B)$p.value, error = function(e) NA_real_)
    data.frame(id = j$sc$id, C = j$sc$C, n = n, rep = r, events = sum(d$y),
               p_v1 = a$p.value, v1_x1 = "x1" %in% a$axes,
               p_all = b$p.value, all_x1 = "x1" %in% b$axes,
               p_HL = unname(pv["HL"]), p_Stukel = unname(pv["Stukel"]), p_proj = pj)
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

## a second runner can take the chunks from the other end (US_REV=1); a chunk both reach is computed
## twice with the same seeds and written atomically, so the result is the same file
if (Sys.getenv("US_REV") == "1") jobs <- rev(jobs)
cat(sprintf("ushape%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else 8L)
  parallel::clusterExport(cl, c("LIB", "B", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk))
  parallel::stopCluster(cl)
}
R <- do.call(rbind, lapply(list.files(OUT, "\\.csv$", full.names = TRUE), utils::read.csv))
cat("rows:", nrow(R), " NA:", colSums(is.na(R[, grep("^p_", names(R))])), "\n")
print(R)
