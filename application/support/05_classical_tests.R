# stage1_tests_SUPPORT.R -- the paper's comparison set on SUPPORT, for M0 (linear meanbp) and the repaired
# model M1 = ns(meanbp, 3) (stage1_repair_SUPPORT.R, repair (b)). Replaces the run.all.gof() call of the
# screen, which also runs tests outside the paper (le Cessie, the GAM tests) and did not finish in 70 min
# at n = 8,873. The projection test is not run: its O(n^3) weight is refused above n = 3,000 by the
# package (gof_proj max_n). DeepGOF-1 values are the ones already computed with the fixed seeds.
suppressMessages({ library(ebrahim.gof); library(splines) })
source("data_prep.R")
out <- file.path(HERE, "results", "screen_SUPPORT")
fits <- list(M0 = readRDS(file.path(out, "M0.rds")), M1 = readRDS(file.path(out, "M1_b.rds"))$fit)
pv <- function(z) { v <- z$p_value; if (length(v) == 1 && is.finite(v)) as.numeric(v) else NA_real_ }
res <- data.frame()
for (m in names(fits)) {
  fit <- fits[[m]]
  ctx <- ebrahim.gof:::.gof_context(fit, G = 10)
  set.seed(20261931L)
  t0 <- Sys.time()
  r <- c(HL = pv(ebrahim.gof:::gof_hl(ctx)),
         mHL = pv(ebrahim.gof:::gof_ftest(ctx)),
         `Pigeon-Heyse` = pv(ebrahim.gof:::gof_ph_test(ctx)),
         `HL equal-width` = pv(ebrahim.gof:::gof_hlw(ctx)),
         `Osius-Rojek` = pv(ebrahim.gof:::gof_osius(ctx)),
         Stukel = pv(ebrahim.gof:::gof_stukel(ctx, list(form = "marginal"))),
         GiViTI = givitiR::givitiCalibrationTest(o = as.numeric(ctx$y), e = as.numeric(ctx$ph),
                                                 devel = "internal")$p.value,
         `Stute-Zhu` = pv(ebrahim.gof:::gof_stutezhu(ctx)))
  cat(m, "classical tests done (", format(Sys.time() - t0), ")\n")
  res <- rbind(res, data.frame(model = m, test = names(r), p_value = as.numeric(r)))
  utils::write.csv(res, file.path(out, "stage1_tests.csv"), row.names = FALSE)
}
fx <- readRDS(file.path(out, "stage1_fixedaxes.rds"))
res <- rbind(res, data.frame(model = c("M0", "M1"), test = "DeepGOF-1 (axes meanbp x num.co)",
                             p_value = c(fx$M0$p, fx$M1$p)))
utils::write.csv(res, file.path(out, "stage1_tests.csv"), row.names = FALSE)
print(stats::reshape(res, idvar = "test", timevar = "model", direction = "wide"), row.names = FALSE)
