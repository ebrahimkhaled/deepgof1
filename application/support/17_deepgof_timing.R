# dg_timing_SUPPORT.R -- referee M6 companion: wall-clock time of deepgof1() on the SUPPORT model M0 at
# n = 8,873, B = 199, for each reading, on the same machine and under the same load as projection_SUPPORT.R,
# so the two times can be compared directly. Output: screen_SUPPORT/deepgof_timing_M0.csv
LIB <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
suppressMessages(library(ebrahim.gof, lib.loc = LIB))
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/medical_example/data_prep.R")
out <- file.path(HERE, "screen_SUPPORT")
M0 <- readRDS(file.path(out, "M0.rds"))
rows <- list()
for (rd in c("axes", "allpairs", "combined")) {
  set.seed(20261960L)
  t0 <- Sys.time()
  r <- deepgof1(M0, B = 199L, reading = rd)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rows[[rd]] <- data.frame(model = "M0", n = nrow(stats::model.matrix(M0)), B = 199, reading = rd,
                           p = r$p.value, seconds = secs)
  print(rows[[rd]]); flush.console()
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(out, "deepgof_timing_M0.csv"), row.names = FALSE)
