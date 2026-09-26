# check_separation.R -- inclusion gate "no separation" for the data sets that pass the metadata screen.
# Prints ONLY separation diagnostics (convergence, extreme fitted risks, largest |b_j| sd_j) and the event
# count. No coefficient table, no goodness-of-fit test. Rule (declared in PREDECLARATION.md (a)):
# separation (or quasi-separation) is declared if glm(maxit = 100) does not converge, or any fitted risk
# is below 1e-8 or above 1 - 1e-8, or the largest |b_j| sd_j exceeds 15.
suppressMessages({ library(aplore3); library(MASS); library(mlbench); library(carData) })
chk <- function(name, d, f) {
  g <- suppressWarnings(glm(f, data = d, family = binomial(), control = glm.control(maxit = 100)))
  mm <- model.matrix(g); b <- coef(g)[-1]; s <- apply(mm[, -1, drop = FALSE], 2, sd)
  ph <- fitted(g); ext <- sum(ph < 1e-8 | ph > 1 - 1e-8); mx <- max(abs(b) * s, na.rm = TRUE)
  sep <- !g$converged || ext > 0 || !is.finite(mx) || mx > 15 || any(is.na(b))
  cat(sprintf("%-12s n = %4d  events = %4d  p = %2d  converged = %s  extreme fitted = %3d  aliased = %d  max|b|sd = %8.2f  SEPARATION = %s\n",
              name, nrow(d), sum(model.response(model.frame(g))), ncol(mm) - 1, g$converged, ext, sum(is.na(b)), mx, sep))
}
m <- aplore3::myopia
dm <- data.frame(y = as.integer(m$myopic == "Yes"), age = m$age, gender = as.integer(m$gender == "Female"),
                 al = m$al, acd = m$acd, lt = m$lt, vcd = m$vcd, sporthr = m$sporthr, readhr = m$readhr,
                 comphr = m$comphr, studyhr = m$studyhr, tvhr = m$tvhr,
                 mommy = as.integer(m$mommy == "Yes"), dadmy = as.integer(m$dadmy == "Yes"))
chk("myopia", dm, y ~ .)
pm <- rbind(MASS::Pima.tr, MASS::Pima.te); dp <- data.frame(y = as.integer(pm$type == "Yes"), pm[, 1:7])
chk("Pima", dp, y ~ .)
mr <- carData::Mroz
dr <- data.frame(y = as.integer(mr$lfp == "yes"), k5 = mr$k5, k618 = mr$k618, age = mr$age,
                 wc = as.integer(mr$wc == "yes"), hc = as.integer(mr$hc == "yes"), lwg = mr$lwg, inc = mr$inc)
chk("Mroz", dr, y ~ .)
data(Ionosphere, package = "mlbench"); io <- Ionosphere
di <- data.frame(y = as.integer(io$Class == "good"), V1 = as.integer(io$V1 == "1"), io[, 3:34])
chk("Ionosphere", di, y ~ .)
data(Sonar, package = "mlbench"); so <- Sonar
ds <- data.frame(y = as.integer(so$Class == "M"), so[, 1:60])
chk("Sonar", ds, y ~ .)
cr <- MASS::crabs
dc <- data.frame(y = as.integer(cr$sp == "O"), sex = as.integer(cr$sex == "M"), cr[, c("FL", "RW", "CL", "CW", "BD")])
chk("crabs_sp", dc, y ~ .)
dc2 <- data.frame(y = as.integer(cr$sex == "M"), sp = as.integer(cr$sp == "O"), cr[, c("FL", "RW", "CL", "CW", "BD")])
chk("crabs_sex", dc2, y ~ .)
# collinearity of the myopia eye-length columns (design only, no outcome): AL against ACD + LT + VCD
cat(sprintf("myopia: cor(AL, ACD + LT + VCD) = %.4f; max |AL - (ACD + LT + VCD)| = %.3f\n",
            cor(dm$al, dm$acd + dm$lt + dm$vcd), max(abs(dm$al - dm$acd - dm$lt - dm$vcd))))
# is any file sorted by the outcome? (first index of the second outcome level)
for (nm in c("myopia", "Pima", "Mroz")) { y <- switch(nm, myopia = dm$y, Pima = dp$y, Mroz = dr$y)
  cat(sprintf("%-7s rows sorted by y: %s  (number of runs of y = %d)\n", nm, !is.unsorted(y) || !is.unsorted(rev(y)), length(rle(y)$lengths))) }
