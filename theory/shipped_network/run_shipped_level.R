# run_shipped_level.R -- the level of the SHIPPED DeepGOF-1 (ebrahim.gof 2.8.0, deepgof1()) on the null
# datasets of the sixty-cell level study (Simulation I). The paper's Simulation I scored those datasets with
# networks retrained per block; this run answers the obvious referee question: what is the level of the
# network that is actually deployed, run exactly as an analyst runs it (plain glm fit, B = 199, the
# deployed axis rule and tie-break)?
#
# Datasets: regenerated from their seeds by the archived harness (benchmark/grid/s3_designs.R,
# gen_design_null, set.seed(SEED0 + i)) with SEED0 = 101000 (the first corpus) and i = 1001..2000, the
# held-out half. One network, so 20 cells (5 n x 4 designs) of 1,000 nulls each.
# Output: rows/<design>_n<n>_c<chunk>.csv, one file per chunk of 50 datasets, written when the chunk ends
# (resume-safe: finished chunks are skipped). Usage: Rscript run_shipped_level.R [smoke]
suppressMessages(library(ebrahim.gof))
stopifnot(utils::packageDescription("ebrahim.gof")$Version == "2.8.0", exists("projection.gof"))
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level"
GRID <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/release/deepgof1_v2_final/benchmark/grid"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
ROWS <- file.path(HERE, if (SMOKE) "rows_smoke" else "rows")
dir.create(ROWS, showWarnings = FALSE, recursive = TRUE)

SEED0 <- 101000L
CFGS <- list(D1 = list(p = 2L, law = "indep", rho = 0.5, cat = 0L, rate = 0.50),
             D2 = list(p = 3L, law = "indep", rho = 0.5, cat = 3L, rate = 0.50),
             D3 = list(p = 3L, law = "ar1",   rho = 0.7, cat = 0L, rate = 0.50),
             D4 = list(p = 2L, law = "indep", rho = 0.5, cat = 0L, rate = 0.15))
NS <- c(50L, 75L, 100L, 150L, 200L)
IDS <- if (SMOKE) 1001:1004 else 1001:2000
CH <- if (SMOKE) 2L else 50L

jobs <- list()
for (des in names(CFGS)) for (n in NS) {
  chunks <- split(IDS, ceiling(seq_along(IDS) / CH))
  for (c in seq_along(chunks))
    jobs[[length(jobs) + 1]] <- list(des = des, n = n, c = c, ids = chunks[[c]],
                                     file = file.path(ROWS, sprintf("%s_n%03d_c%02d.csv", des, n, c)))
}
todo <- Filter(function(j) !file.exists(j$file), jobs)

one_chunk <- function(j) {
  suppressMessages(library(ebrahim.gof))
  source(file.path(GRID, "s3_designs.R"))
  cfg <- CFGS[[j$des]]
  out <- lapply(j$ids, function(i) {
    set.seed(SEED0 + i)                                   # the archived null dataset, exactly
    g <- tryCatch(gen_design_null(j$n, cfg$p, cfg$law, cfg$rho, cfg$cat, cfg$rate), error = function(e) NULL)
    if (is.null(g) || length(unique(g$data$y)) < 2L)
      return(data.frame(design = j$des, n = j$n, id = i, events = NA, converged = NA, p = NA, axes = NA))
    fit <- g$fit
    set.seed(770000L + 1000L * match(j$n, NS) + 100000L * match(j$des, names(CFGS)) + i)
    dg <- tryCatch(deepgof1(fit, B = 199), error = function(e) NULL)
    data.frame(design = j$des, n = j$n, id = i, events = sum(g$data$y), converged = fit$converged,
               p = if (is.null(dg)) NA else dg$p.value,
               axes = if (is.null(dg)) NA else paste(dg$axes, collapse = "|"))
  })
  r <- do.call(rbind, out)
  tmp <- paste0(j$file, ".part"); utils::write.csv(r, tmp, row.names = FALSE); file.rename(tmp, j$file)
  j$file
}

cat(sprintf("shipped level%s | %d of %d chunks to run | %s\n", if (SMOKE) " SMOKE" else "",
            length(todo), length(jobs), format(Sys.time())))
if (length(todo)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else 20L)
  parallel::clusterExport(cl, c("GRID", "CFGS", "NS", "SEED0", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, todo, one_chunk))
  parallel::stopCluster(cl)
}
R <- do.call(rbind, lapply(Filter(file.exists, vapply(jobs, `[[`, "", "file")), utils::read.csv))
cat("rows:", nrow(R), " failed:", sum(is.na(R$p)), " |", format(Sys.time()), "\n")
if (SMOKE) print(R)
