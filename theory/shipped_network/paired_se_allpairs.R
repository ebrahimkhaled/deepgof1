# paired standard errors for the all-pairs option against Stukel, the default reading and the per-block networks,
# on the same datasets as rank_shipped_liu.R: per-dataset randomized rejection at each test's matched threshold
src <- readLines("rank_shipped_liu.R"); src <- src[seq_len(grep("^X <- do.call", src) - 1)]
i <- grep("for (tst in names(P)", src, fixed = TRUE)
src <- c(src[seq_len(i - 1)], "  PP[[length(PP) + 1]] <- cbind(S = S, n = n, P)", "}")
PP <- list(); eval(parse(text = src))
rej <- function(p0, p1, a = .05) { p0 <- p0[!is.na(p0)]; u <- sort(unique(p0)); F <- vapply(u, function(c) mean(p0 <= c), 0)
  k <- which(F >= a)[1]; c <- u[k]; below <- if (k > 1) F[k - 1] else 0; g <- (a - below) / mean(p0 == c)
  ifelse(p1 < c, 1, ifelse(p1 == c, g, 0)) }
D <- list()
for (B in PP) { al <- B$kind == 2
  r <- sapply(c("DeepGOF-1 all-pairs", "Stukel", "DeepGOF-1 shipped", "DeepGOF-1 per-block"),
              function(t) rej(B[[t]][B$kind == 0], B[[t]][al]))
  D[[length(D) + 1]] <- data.frame(S = B$S[1], n = B$n[1], r) }
D <- do.call(rbind, D); names(D)[3:6] <- c("ap", "stukel", "v1", "perblock")
for (o in c("stukel", "v1", "perblock")) { d <- D$ap - D[[o]]; ok <- !is.na(d)
  cat(sprintf("all-pairs minus %-9s %+.3f  (paired SE %.3f, %d datasets)\n", o, mean(d[ok]), sd(d[ok]) / sqrt(sum(ok)), sum(ok))) }
print(aggregate(cbind(ap, v1, perblock) ~ S, D, mean), digits = 3)
