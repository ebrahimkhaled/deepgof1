# check_liu_regen.R -- the regenerated benchmark datasets are the deposited ones (HL p-values match).
suppressMessages(library(ebrahim.gof))
src <- readLines("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/run_shipped_liu.R")
a <- grep("^gen_liu <- function", src); b <- a + grep("^}", src[a:length(src)])[1] - 1
eval(parse(text = src[a:b]))
dep <- read.csv("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/release/deepgof1_v2_final/results/liu/S4_n100_s101000/blk_pvalues.csv")
hl <- function(d) { f <- glm(y ~ x1 + x2, data = d, family = binomial())
  r <- suppressWarnings(run.all.gof(f, tests = "HL", include_slow = FALSE, install = "no")); r$p_value[r$Test == "HL"] }
for (i in 1001:1003) { set.seed(101000L + 4L * 500000L + i); d <- gen_liu(100L, 4L, 0)
  cat("null", i, signif(hl(d), 6), signif(dep$HL[dep$kind == 0 & dep$rep == i], 6), "\n") }
set.seed(101000L + 4L * 500000L + 2L * 50000L + 3L * 5000L + 7L); d <- gen_liu(100L, 4L, 0.5)
cat("alt C=.5 rep 7", signif(hl(d), 6), signif(dep$HL[dep$kind == 2 & dep$rep == 7 & abs(dep$C - 0.5) < 1e-9], 6), "\n")
