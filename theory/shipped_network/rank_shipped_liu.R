# rank_shipped_liu.R -- does the shipped network's power change DeepGOF-1's place among the tests? Every test on
# the SAME datasets (Settings 2-4, n in {50,100,500}, corpus 101000; held-out nulls 1001-2000; the first 500
# alternatives at each C), each at exact matched level .05 by a randomized threshold on those nulls and
# evaluated on the datasets where it returns a p-value (the paper's convention). BAGofT is not included (its
# benchmark rows cover a hundred alternatives per cell).
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory"
DEP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/release/deepgof1_v2_final/results/liu"
R <- do.call(rbind, lapply(list.files(file.path(HERE, "shipped_level", "rows_liu"), "\\.csv$", full.names = TRUE), utils::read.csv))
matched <- function(p0, p1, a = .05) {
  p0 <- p0[!is.na(p0)]; ok <- !is.na(p1); u <- sort(unique(p0)); F <- vapply(u, function(c) mean(p0 <= c), 0)
  k <- which(F >= a)[1]; c <- u[k]; below <- if (k > 1) F[k - 1] else 0; g <- (a - below) / mean(p0 == c)
  mean(ifelse(p1[ok] < c, 1, ifelse(p1[ok] == c, g, 0)))
}
key <- function(k, r, C) paste(k, r, ifelse(is.na(C), "", sprintf("%.2f", C)))
res <- list()
for (S in 2:4) for (n in c(50, 100, 500)) {
  b <- R[R$setting == S & R$n == n, ]; b$key <- key(b$kind, b$rep, b$C)
  blk <- sprintf("S%d_n%d_s101000", S, n)
  dp <- utils::read.csv(file.path(DEP, blk, "deep_pvalues.csv")); dp$key <- key(dp$kind, dp$rep, dp$C)
  cl <- utils::read.csv(file.path(DEP, blk, "blk_pvalues.csv"), check.names = FALSE); cl$key <- key(cl$kind, cl$rep, cl$C)
  au <- utils::read.csv(file.path(DEP, blk, "augment_pvalues.csv")); au$key <- key(au$kind, au$rep, au$C)
  ad <- utils::read.csv(file.path(HERE, "addon_arms", "liu_rows", paste0(blk, ".csv"))); ad$key <- key(ad$kind, ad$rep, ad$C)
  pr <- do.call(rbind, lapply(list.files(file.path(HERE, "addon_arms", "proj_rows_liu", blk), full.names = TRUE), utils::read.csv))
  pr$key <- key(pr$kind, pr$rep, pr$C)
  ## ebrahim.gof 2.9.0's default (all-pairs) reading: Setting 2 from run_allpairs_liu.R; with two covariates
  ## (Settings 3 and 4) it is bit-identical to the 2.8.0 rule, so those p-values are the 2.8.0 ones
  ap <- if (S == 2) {
    a <- do.call(rbind, lapply(list.files(file.path(HERE, "shipped_level", "rows_liu_ap"),
                                          sprintf("^S2_n%03d_.*\\.csv$", n), full.names = TRUE), utils::read.csv))
    a$key <- key(a$kind, a$rep, a$C); a$p_allpairs[match(b$key, a$key)]
  } else b$p_shipped
  P <- data.frame(key = b$key, kind = b$kind,
    `DeepGOF-1 all-pairs` = ap,
    `DeepGOF-1 shipped` = b$p_shipped, `DeepGOF-1 per-block` = dp$p_deep[match(b$key, dp$key)],
    projection = pr$proj_liu_raw[match(b$key, pr$key)], Stukel = cl$Stukel[match(b$key, cl$key)],
    GiViTI = ad$GiViTI[match(b$key, ad$key)], HL = cl$HL[match(b$key, cl$key)],
    mHL = ad$mHL[match(b$key, ad$key)], `Stute-Zhu` = au$SZ[match(b$key, au$key)],
    `Pigeon-Heyse` = cl$`Pigeon-Heyse`[match(b$key, cl$key)], `HL equal-width` = cl$`HL-equalwidth`[match(b$key, cl$key)],
    check.names = FALSE)
  for (tst in names(P)[-(1:2)]) res[[length(res) + 1]] <- data.frame(setting = S, n = n, test = tst,
    power = matched(P[[tst]][P$kind == 0], P[[tst]][P$kind == 2]), size = mean(P[[tst]][P$kind == 0] <= .05, na.rm = TRUE),
    missing = mean(is.na(P[[tst]])))
}
X <- do.call(rbind, res)
O <- aggregate(cbind(power, size) ~ test, X, mean)
O <- O[order(-O$power), ]
cat("Settings 2-4, same datasets, matched level .05, averaged over the 9 blocks:\n")
print(format(O, digits = 3), row.names = FALSE)
cat("\nby n:\n"); print(format(stats::reshape(aggregate(power ~ test + n, X, mean), idvar = "test", timevar = "n",
                                           direction = "wide")[order(-O$power[match(unique(X$test), O$test)]), ], digits = 3), row.names = FALSE)
utils::write.csv(X, file.path(HERE, "shipped_level", "RANK_SHIPPED_LIU.csv"), row.names = FALSE)
