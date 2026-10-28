# 30_verify.R -- fingerprint check of rows_main: (1) regenerate 400 random data sets from their recorded seeds
# and assert events and xfp (sum of all covariates) are reproduced exactly; (2) the smoke rows that share a
# seed with the main run (nulls, reps 1-2) must agree on every p-value (the whole pipeline is reproducible).
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3b/b_common.R")
rd <- function(d) do.call(rbind, lapply(list.files(file.path(T3B, d), "\\.csv$", full.names = TRUE), utils::read.csv, stringsAsFactors = FALSE))
A <- rd("rows_main"); S <- rd("rows_smoke"); HM <- h_means()
seed_at(4242L); idx <- sample.int(nrow(A), min(400L, nrow(A)))
ok <- vapply(idx, function(i) { r <- A[i, ]
  stopifnot(r$dseed == data_seed(r$design, r$shape, r$n, r$ci, r$rep), r$bseed == boot_seed(r$dseed))
  seed_at(r$dseed); d <- gen_b(r$design, r$shape, r$n, r$C, HM)
  sum(d$y) == r$events && sprintf("%.17g", sum(as.matrix(d[, setdiff(names(d), "y")]))) == sprintf("%.17g", r$xfp) }, TRUE)
cat(sprintf("regenerated %d data sets from recorded seeds: events and xfp reproduced on %d\n", length(idx), sum(ok)))
## C recorded in the rows equals C_CHOSEN.csv
CC <- utils::read.csv(file.path(TABB, "C_CHOSEN.csv")); alt <- A[A$shape != "null", ]
m <- match(paste(alt$design, alt$shape, alt$ci), paste(CC$design, CC$shape, CC$ci))
cat(sprintf("alternatives whose C equals C_CHOSEN.csv: %d of %d\n", sum(abs(alt$C - CC$C[m]) < 1e-9), nrow(alt)))
k <- function(z) paste(z$design, z$shape, z$n, z$ci, z$rep)
S0 <- S[S$shape == "null", ]; mm <- match(k(S0), k(A))
cols <- c("events", "xfp", "p_dg", "p_ss", "p_score", "p_score_asym", "p_hl", "p_stukel_fb", "p_grp")
same <- sapply(cols, function(v) sum(S0[[v]] == A[[v]][mm] | (is.na(S0[[v]]) & is.na(A[[v]][mm])), na.rm = TRUE))
cat(sprintf("smoke null rows found in main: %d of %d; identical per column:\n", sum(!is.na(mm)), nrow(S0))); print(same)
