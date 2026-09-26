# =============================================================================
# common.R -- shared pieces for the EXPLORATORY BAGofT small-n trial.
#
# *** EXPLORATORY, NOT CONFIRMATORY. ***
# These numbers may never be published as a head-to-head. If the result looks
# promising a separate PRE-DECLARED confirmatory run on FRESH seeds must follow.
#
# Nothing here is reimplemented: gen_liu() and the dataset seeding are copied
# VERBATIM from paper_deepgof/code/liu_grid.R so the datasets are bit-identical
# to the ones the existing DeepGOF-1 numbers were computed on, and the BAGofT
# call itself is the registry's verified wrapper bagoft_p() at PUBLISHED
# DEFAULTS (nsplits=100, nsim=100, parRF(ntree=60), ne=floor(5*sqrt(n))).
# =============================================================================

ROOT <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def"
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
EXPL <- file.path(PAP, "explore_bagoft")

source(file.path(ROOT, "jasa_study", "config.R"))
source(file.path(ROOT, "jasa_study", "R", "methods_registry.R"))   # <- bagoft_p()
suppressMessages(library(BAGofT))

## ---- Liu's four DGPs, VERBATIM from liu_grid.R lines 55-72 ------------------
gen_liu <- function(n, setting, C) {
  if (setting == 1L) {
    x1 <- runif(n,-3,3); d <- data.frame(x1)
    eta <- -1.963 + 0.982*x1 + C*x1^2
  } else if (setting == 2L) {
    x1 <- runif(n,-3,3); x2 <- rnorm(n,0,1.5); x3 <- rt(n,4)
    d <- data.frame(x1,x2,x3)
    eta <- -2 + 0.3*x1 + 0.3*x2 + 0.3*x3 + C*x1^2
  } else if (setting == 3L) {
    x1 <- runif(n,-3,3); x2 <- rbinom(n,1,.5); d <- data.frame(x1,x2)
    eta <- -1.972 + 0.135*x1 + 1.117*x2 + C*x1*x2
  } else {
    x1 <- runif(n,-3,3); x2 <- runif(n,-3,3); d <- data.frame(x1,x2)
    eta <- 0.3*x1 + 0.3*x2 + C*x1*x2                      # no intercept in the truth
  }
  d$y <- rbinom(n, 1, plogis(eta)); d
}
FRM <- function(d) stats::as.formula(paste("y ~", paste(setdiff(names(d),"y"), collapse=" + ")))

## ---- the EXACT dataset seeding of liu_grid.R::one() -------------------------
## one(): set.seed(SEED0 + SET*500000L + kind*50000L + ci*5000L + i)
##   kind = 0 -> the NULL block  (C = 0, ci = 0),  i = 1..NNULL (2000)
##   kind = 2 -> the C > 0 block, ci = index of C in CS, i = 1..NEVAL (1000)
LIU_CS <- c(0, 0.25, 0.5, 0.75, 1.0)

liu_seed <- function(SEED0, SET, C, i) {
  if (C == 0) {
    kind <- 0L; ci <- 0L
  } else {
    kind <- 2L; ci <- match(C, LIU_CS)
    if (is.na(ci)) stop("C not on the pre-declared Liu grid: ", C)
  }
  as.integer(SEED0) + as.integer(SET)*500000L + kind*50000L + as.integer(ci)*5000L + as.integer(i)
}

## returns the dataset liu_grid.R would have used, or NULL if liu_grid.R would
## have DROPPED it (degenerate response) -- same filter, same datasets.
liu_dataset <- function(SEED0, SET, n, C, i) {
  s <- liu_seed(SEED0, SET, C, i)
  set.seed(s)
  d <- gen_liu(n, as.integer(SET), C)
  if (length(unique(d$y)) < 2L) return(list(seed = s, data = NULL))
  list(seed = s, data = d)
}

## ---- BAGofT's own split arithmetic at these n (no simulation needed) --------
bag_geometry <- function(n) {
  ne   <- floor(5 * sqrt(n))        # validation set size  (published default)
  nt   <- n - ne                    # training set size
  nmin <- ceiling(sqrt(ne))         # parRF(): min cell count in the validation set
  Kmax <- floor(ne / nmin)          # parRF(): max number of groups searched
  data.frame(n = n, ne = ne, n_train = nt, nmin = nmin, Kmax = Kmax,
             frac_train = round(nt / n, 3))
}
