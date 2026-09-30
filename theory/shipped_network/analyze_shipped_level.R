# analyze_shipped_level.R -- size of the SHIPPED DeepGOF-1 on the held-out nulls of the level study, by cell,
# against the pre-declared band [.035, .065] of Simulation I. Rejection at p <= .05 (the deployed rule).
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level"
R <- do.call(rbind, lapply(list.files(file.path(HERE, "rows"), "\\.csv$", full.names = TRUE), utils::read.csv))
cat("datasets:", nrow(R), " with a p-value:", sum(!is.na(R$p)), " non-converged fits:", sum(!R$converged, na.rm = TRUE), "\n")
S <- do.call(rbind, lapply(split(R, list(R$design, R$n), drop = TRUE), function(d) {
  ok <- !is.na(d$p); s <- mean(d$p[ok] <= .05)
  data.frame(design = d$design[1], n = d$n[1], nulls = sum(ok), size = s, se = sqrt(s * (1 - s) / sum(ok)),
             in_band = s >= .035 & s <= .065, axes_top = names(sort(table(d$axes), decreasing = TRUE))[1])
}))
S <- S[order(S$design, S$n), ]
print(format(S, digits = 3), row.names = FALSE)
cat(sprintf("\ncells in band: %d of %d | mean size %.4f | range %.4f-%.4f | pooled size %.4f (SE %.4f)\n",
            sum(S$in_band), nrow(S), mean(S$size), min(S$size), max(S$size),
            mean(R$p[!is.na(R$p)] <= .05), sqrt(.05 * .95 / sum(!is.na(R$p)))))
utils::write.csv(S, file.path(HERE, "SHIPPED_LEVEL.csv"), row.names = FALSE)
