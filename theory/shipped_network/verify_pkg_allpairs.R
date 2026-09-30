# verify_pkg_allpairs.R -- does ebrahim.gof 2.9.0's deepgof1() (reading = "allpairs") reproduce the
# p-values of the shipped-network study in phd/ (all-pairs on Setting 2; the one map on Settings 1 and 4)?
# Same seeds, same generator. The study refitted separated samples with Firth's method and scored a
# failed refit 0; the package refits with glm() only, so rows whose fit or bootstrap needed Firth can
# differ. Those are reported, not hidden.
PKG <- "C:/Users/ebrah/.cursor-tutor/projects/ebrahim.gof-proj-bagoft"
PAP <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
suppressMessages(devtools::load_all(PKG, quiet = TRUE))
source(file.path(PAP, "phd", "run_phd_defs.R"))           # gen_liu, K
CS_ALL <- c(0, 0.25, 0.5, 0.75, 1.0)
NPER <- as.integer(Sys.getenv("VP_NPER", "12"))

check <- function(SET, NOBS, col, rowsdir, file) {
  R <- utils::read.csv(file.path(PAP, "phd", rowsdir, file))
  R <- R[R$firth == 0, ]
  pick <- rbind(utils::head(R[R$kind == 0, ], NPER), utils::head(R[R$kind == 2 & R$C == 0.5, ], NPER))
  do.call(rbind, lapply(seq_len(nrow(pick)), function(k) {
    w <- pick[k, ]; ci <- if (w$kind == 0) 0L else match(w$C, CS_ALL)
    set.seed(w$seed + SET * 500000L + w$kind * 50000L + ci * 5000L + w$rep)
    d <- gen_liu(NOBS, SET, if (w$kind == 0) 0 else w$C)
    xn <- setdiff(names(d), "y")
    fit <- glm(stats::as.formula(paste("y ~", paste(xn, collapse = " + "))), data = d, family = binomial())
    p <- deepgof1(fit)$p.value
    data.frame(setting = SET, n = NOBS, kind = w$kind, rep = w$rep, study = w[[col]], package = p)
  }))
}
out <- rbind(check(2L, 100L, "allpairs", "rows", "S2_n100_s101000.csv"),
             check(2L, 500L, "allpairs", "rows", "S2_n500_s101000.csv"),
             check(4L, 100L, "v1", "rows", "S4_n100_s101000.csv"),
             check(1L, 100L, "v1", "rows_s1", "S1_n100_s101000.csv"))
out$same <- abs(out$study - out$package) < 1e-12
print(out, row.names = FALSE)
cat(sprintf("\nidentical p-values: %d of %d\n", sum(out$same), nrow(out)))
print(aggregate(same ~ setting + n, out, function(z) sprintf("%d/%d", sum(z), length(z))))
utils::write.csv(out, file.path(PAP, "theory", "shipped_level", "VERIFY_PKG_ALLPAIRS.csv"), row.names = FALSE)
