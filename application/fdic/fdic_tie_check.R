## fdic_tie_check.R -- the FDIC adequacy calls under both tie rules, against the deposited values.
## Same fits and seeds as fdic_deepgof.R (771000 + cell index), B = 199. Takes about a minute.
## The deposit (results/fdic/deepgof_observed.csv) holds the fixed rule (ties = "random"), which is
## what the paper reports; the row-order rule gives the values computed before the fix.
## Run from the archive root:  Rscript application/fdic/fdic_tie_check.R
## Output recorded in application/fdic/fdic_tie_check.log.
ARC <- normalizePath(Sys.getenv("DEEPGOF_ARCHIVE", "."))
source(file.path(ARC, "R", "deepgof1.R"))
d  <- read.csv(file.path(ARC, "application", "fdic", "data", "fdic_2008_banks.csv"), stringsAsFactors = FALSE)
st <- read.csv(file.path(ARC, "application", "fdic", "data", "fdic_state_counts.csv"), stringsAsFactors = FALSE)
CELLS <- c("GA", setdiff(st$state[st$failures >= 15], "GA"))
FRM <- failed ~ eq_ta + npa_ta + roa + bro_dep + cre_ta
ROWORDER <- c(GA = .320, FL = .185, IL = .805, CA = .710, WA = .350, MN = .905)  # before the fix
dep <- read.csv(file.path(ARC, "results", "fdic", "deepgof_observed.csv"))
ok <- TRUE
for (S in CELLS) {
  ds  <- d[d$state == S, ]
  fit <- suppressWarnings(glm(FRM, data = ds, family = binomial()))
  set.seed(771000L + match(S, CELLS)); rn <- deepgof1(fit, ties = "random")
  set.seed(771000L + match(S, CELLS)); ro <- deepgof1(fit, ties = "roworder")
  pd <- dep$p_deepgof[dep$state == S]; sd_ <- dep$S_obs[dep$state == S]
  same_n <- isTRUE(all.equal(rn$p.value, pd)) && isTRUE(all.equal(rn$statistic, sd_))
  same_o <- isTRUE(all.equal(ro$p.value, ROWORDER[[S]]))
  ok <- ok && same_n && same_o
  cat(sprintf("%s  random ties: S %.6f p %.3f (deposited S %.6f p %.3f, %s) | row-order ties: p %.3f (%s)\n",
              S, rn$statistic, rn$p.value, sd_, pd, if (same_n) "identical" else "DIFFERENT",
              ro$p.value, if (same_o) "as before the fix" else "DIFFERENT"))
}
cat(if (ok) "fixed rule reproduces every deposited FDIC adequacy call; row-order rule reproduces the pre-fix values\n" else "MISMATCH\n")
