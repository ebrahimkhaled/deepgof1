# data_prep.R -- the SUPPORT data and the linear in-hospital mortality model M0 (Section 9).
# prep(2) returns the complete cases, permuted once by the fixed seed, and the formula of M0.
# The index 2 is kept so that every seed matches the one used in the paper (PERM_SEED + k, DG_SEED + k).
HERE <- "."
DNAMES <- c(NA, "SUPPORT")
PERM_SEED <- 20260929L

prep <- function(k = 2L) {
  stopifnot(k == 2L)
  s <- utils::read.csv(file.path(HERE, "data", "support2.csv"), na.strings = c("NaN", "NA", ""))
  d <- data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt,
                  resp = s$resp, temp = s$temp, crea = s$crea, sod = s$sod, wblc = s$wblc,
                  num.co = s$num.co, sex = factor(s$sex))
  d <- stats::na.omit(d)
  set.seed(PERM_SEED + k)
  d <- d[sample(nrow(d)), ]
  rownames(d) <- NULL
  list(data = d, formula = stats::as.formula(paste("y ~", paste(setdiff(names(d), "y"), collapse = " + "))))
}
