# =============================================================================
# calibrate.R -- the pre-declared calibration of C (PREDECLARATION.md (d)).
# Pilot arm HL only. For each (data set, design, family), pid = 1..8 in the order
# myopia/Pima x niche/dense x inter/quad; grid C = 0.1..3.0 (gi = 1..30); 200 replicates each,
# seed 91000000 + 100000 pid + 1000 gi + r. HL raw rejection at .05 on the fitted main-effects model.
# C1 = smallest grid C whose pilot power and that of the next grid value are both >= .10;
# C3 = the same with >= .40 (3.0 if never reached); C2 = grid value nearest (C1 + C3)/2, lower on a tie;
# strict order enforced by +0.1. Writes pilot/pilot_<pid>.csv, calibration.csv, scenarios.csv.
# Env: CAL_W (workers, default 3).
# =============================================================================
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/real_data_niche/common.R")
NW <- as.integer(Sys.getenv("CAL_W", "3"))
dir.create(file.path(RD, "pilot"), showWarnings = FALSE)
GRID <- round(seq(0.1, 3.0, by = 0.1), 1)
PT <- expand.grid(fam = c("inter", "quad"), design = c("niche", "dense"), k = 1:2, stringsAsFactors = FALSE)[, c("k", "design", "fam")]
PT$pid <- seq_len(nrow(PT))

todo <- PT$pid[!file.exists(file.path(RD, "pilot", sprintf("pilot_%d.csv", PT$pid)))]
if (length(todo)) {
  cl <- parallel::makeCluster(min(NW, length(todo)))
  parallel::clusterExport(cl, c("PT", "GRID"))
  invisible(parallel::clusterEvalQ(cl, source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/real_data_niche/common.R")))
  parallel::parLapplyLB(cl, todo, function(pid) {
    pt <- PT[PT$pid == pid, ]; D <- design_setup(pt$k)
    out <- do.call(rbind, lapply(seq_along(GRID), function(gi) {
      rej <- vapply(1:200, function(r) {
        set.seed(91000000L + 100000L * pid + 1000L * gi + r)
        dat <- data.frame(y = gen_y(D, pt$design, pt$fam, GRID[gi]), D$X); f <- fit_pol(dat)
        if (is.null(f)) return(0)
        p <- tryCatch(pv(ebrahim.gof:::gof_hl(ebrahim.gof:::.gof_context(f, G = 10))), error = function(e) NA_real_)
        as.numeric(is.finite(p) && p <= 0.05)
      }, 0)
      data.frame(pid = pid, data = D$name, design = pt$design, fam = pt$fam, C = GRID[gi], R = 200, hl_power = mean(rej))
    }))
    f <- file.path(RD, "pilot", sprintf("pilot_%d.csv", pid)); utils::write.csv(out, paste0(f, ".part"), row.names = FALSE); file.rename(paste0(f, ".part"), f)
  })
  parallel::stopCluster(cl)
}

pick <- function(pw, thr) {                    # smallest grid C with power >= thr at C and the next grid value
  ok <- which(pw >= thr & c(pw[-1], Inf) >= thr)
  if (length(ok)) GRID[min(ok)] else NA_real_
}
CAL <- do.call(rbind, lapply(PT$pid, function(pid) {
  P <- utils::read.csv(file.path(RD, "pilot", sprintf("pilot_%d.csv", pid))); P <- P[order(P$C), ]
  c1 <- pick(P$hl_power, 0.10); c3 <- pick(P$hl_power, 0.40); c3_reached <- is.finite(c3)
  if (!is.finite(c3)) c3 <- 3.0
  if (!is.finite(c1)) c1 <- 3.0
  at <- function(c) { v <- P$hl_power[abs(P$C - c) < 1e-9]; if (length(v)) v else NA_real_ }
  mid <- (c1 + c3) / 2; dd <- abs(GRID - mid); c2 <- GRID[which(dd == min(dd))[1]]
  if (c2 <= c1) c2 <- round(c1 + 0.1, 1)
  if (c3 <= c2) c3 <- round(c2 + 0.1, 1)
  data.frame(pid = pid, data = P$data[1], design = P$design[1], fam = P$fam[1], C1 = c1, C2 = c2, C3 = c3, C3_target_reached = c3_reached,
             hl_at_C1 = at(c1), hl_at_C2 = at(c2), hl_at_C3 = at(c3))
}))
utils::write.csv(CAL, file.path(RD, "calibration.csv"), row.names = FALSE)
print(CAL, row.names = FALSE)

## scenario table: ids in the order data (myopia, Pima) x design (niche, dense) x (null, inter C1..C3, quad C1..C3)
SC <- list(); id <- 0L
for (k in 1:2) for (design in c("niche", "dense")) {
  id <- id + 1L; SC[[length(SC) + 1]] <- data.frame(id = id, k = k, data = DNAMES[k], design = design, fam = "null", level = 0, C = 0, reps = 1000L)
  for (fam in c("inter", "quad")) { cc <- CAL[CAL$data == DNAMES[k] & CAL$design == design & CAL$fam == fam, ]
    for (lv in 1:3) { id <- id + 1L
      SC[[length(SC) + 1]] <- data.frame(id = id, k = k, data = DNAMES[k], design = design, fam = fam, level = lv, C = cc[[paste0("C", lv)]], reps = 500L) } }
}
SC <- do.call(rbind, SC)
SC$scen <- ifelse(SC$fam == "null", sprintf("%s_%s_null", SC$data, SC$design), sprintf("%s_%s_%s_C%d", SC$data, SC$design, SC$fam, SC$level))
utils::write.csv(SC, file.path(RD, "scenarios.csv"), row.names = FALSE)
print(SC, row.names = FALSE)
