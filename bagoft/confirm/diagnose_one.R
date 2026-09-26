# The one dataset of 240 where the RELEASED BAGofT package did not return a value while the
# reimplementation did (Setting 2, n = 50, C = 1, rep 35). The paper claims the two are the same
# computation, so a disagreement -- even one about whether an answer exists at all -- has to be
# named rather than averaged away.
PAP  <- "/path/to/project/paper_deepgof"
EXPL <- file.path(PAP, "explore_bagoft")
source(file.path(EXPL, "common.R"))
source(file.path(EXPL, "parRF_dropfix.R"))
source(file.path(EXPL, "bagoft_fast.R"))
suppressPackageStartupMessages(library(BAGofT))

dd  <- liu_dataset(101000L, 2L, 50L, 1.0, 35L)
dat <- as.data.frame(lapply(dd$data, as.vector))
frm <- FRM(dat); environment(frm) <- globalenv()
cat("n =", nrow(dat), " events =", sum(dat$y), " of", nrow(dat), "\n")
cat("covariate ranges:\n"); print(summary(dat[, setdiff(names(dat), "y")]))
g <- suppressWarnings(glm(frm, data = dat, family = binomial()))
cat("\nglm converged:", g$converged, "  max |beta| (non-intercept):",
    round(max(abs(coef(g)[-1])), 3), "\n")

set.seed(dd$seed + 500000000L)
cat("\n--- released BAGofT 1.0.0 at published defaults ---\n")
res <- try(invisible(capture.output(out <- suppressMessages(suppressWarnings(
  BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = frm, link = "logit"),
                 parFun = BAGofT::parRF(ntree = 60L), data = dat,
                 nsplits = 100L, nsim = 100L, ne = floor(5 * sqrt(nrow(dat)))))))), silent = TRUE)
if (inherits(res, "try-error")) {
  cat("ERROR:\n"); cat(as.character(res))
} else {
  cat("returned: p =", out$p.value, " p.value2 =", out$p.value2, " pmin =", out$pmin, "\n")
}

set.seed(dd$seed + 500000000L)
cat("\n--- bagoft_fast (the deposit's implementation) ---\n")
f <- try(bagoft_fast(dat, frm, nsplits = 100L, nsim = 100L, ne = floor(5 * sqrt(nrow(dat))),
                     ntree = 60L), silent = TRUE)
if (inherits(f, "try-error")) cat("ERROR:\n", as.character(f)) else
  cat("returned: p =", f$p.value, " p.value2 =", f$p.value2, " pmin =", f$pmin, "\n")
