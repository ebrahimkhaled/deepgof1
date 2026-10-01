# support_ladder_290d.R -- the SUPPORT modelling ladder read by the final 2.9.0 deepgof1(): the axis rule over
# COVARIATES (each covariate scored by the SD of its terms' contribution to the linear predictor), the all-pairs
# reading and the combined reading, all from one deepgof1() call each. Replaces dg290_models.csv, which used the
# earlier column rule (its "v1" picked two spline columns of one covariate after the repair).
# Models M0, M1 (= M1_b), M2, M3; same seed as support_allpairs.R. Output screen_SUPPORT/dg290d_ladder.csv
LIB <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
suppressMessages({ library(ebrahim.gof, lib.loc = LIB); library(splines) })
stopifnot(utils::packageVersion("ebrahim.gof", lib.loc = LIB) == "2.9.0")
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/medical_example/data_prep.R")
out <- file.path(HERE, "screen_SUPPORT")
DG_SEED <- 20261929L + 2L
fits <- list(M0 = readRDS(file.path(out, "M0.rds")), M1 = readRDS(file.path(out, "M1_b.rds"))$fit,
             M2 = readRDS(file.path(out, "M2.rds")), M3 = readRDS(file.path(out, "M3.rds")))
res <- list()
for (m in names(fits)) {
  t0 <- Sys.time(); set.seed(DG_SEED)
  dg <- deepgof1(fits[[m]], B = 199, reading = "combined")
  stopifnot(all(is.finite(dg$boot)))
  cp <- dg$components
  res[[m]] <- data.frame(model = m, p_axes = cp[["axes"]], p_allpairs = cp[["allpairs"]], p_combined = dg$p.value,
                         axes = paste(dg$axes, collapse = " & "),
                         top_pair = paste(unlist(dg$pairs[which.max(dg$pairs$score), 1:2]), collapse = " & "),
                         seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  print(res[[m]], row.names = FALSE); flush.console()
}
utils::write.csv(do.call(rbind, res), file.path(out, "dg290d_ladder.csv"), row.names = FALSE)
