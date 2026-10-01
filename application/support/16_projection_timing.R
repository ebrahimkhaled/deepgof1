# projection_SUPPORT.R -- referee M6: run the projection test of Liu et al. (2024) on the SUPPORT model M0 at
# n = 8,873 and record the wall-clock time and peak memory, instead of asserting that it is impractical.
# projection.gof(fit, B = 199), ebrahim.gof 2.9.0 (private library); M0 is linear in the ten raw covariates,
# which is the test as Liu et al. define it. Output: screen_SUPPORT/projection_M0.csv
LIB <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
suppressMessages(library(ebrahim.gof, lib.loc = LIB))
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/medical_example/data_prep.R")
out <- file.path(HERE, "screen_SUPPORT")
M0 <- readRDS(file.path(out, "M0.rds"))
cat("n =", nrow(stats::model.matrix(M0)), " start", format(Sys.time()), "\n"); flush.console()
gc(reset = TRUE)
set.seed(20261960L)
t0 <- Sys.time()
r <- projection.gof(M0, B = 199)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
peak <- sum(gc()[, ncol(gc())])            # "max used" in MB, summed over cons and vector cells
res <- data.frame(model = "M0", n = nrow(stats::model.matrix(M0)), B = 199, p = r$p.value,
                  statistic = unname(r$statistic), seconds = secs, peak_MB = peak)
print(res)
utils::write.csv(res, file.path(out, "projection_M0.csv"), row.names = FALSE)
