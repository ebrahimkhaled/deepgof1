# glow_roworder.R -- does the fix remove the order dependence on aplore3::glow500?
#   part A: new statistic, file order + the same 30 random orders as posthoc_roworder.R, seed 20260925 + 1.
#   part B (in law): new statistic on the FILE order with 100 different seeds, against the new
#           statistic on 100 random orders with the fixed seed. The two sets should have one law.
# Resumable: one row per unit in glow_roworder.csv.
source("/path/to/project/paper_deepgof/theory/tie_fix/common.R")
suppressMessages(library(aplore3))
g <- aplore3::glow500
d <- data.frame(y = as.integer(g$fracture == "Yes"), age = g$age, height = g$height,
                priorfrac = as.integer(g$priorfrac == "Yes"), momfrac = as.integer(g$momfrac == "Yes"),
                armassist = as.integer(g$armassist == "Yes"), raterisk3 = as.integer(g$raterisk == "Greater"))
f1 <- y ~ age + height + priorfrac + momfrac + armassist + raterisk3
f2 <- y ~ age + height + priorfrac + momfrac + armassist + raterisk3 + age:priorfrac + momfrac:armassist
set.seed(20260925L + 500L)
perms <- replicate(30, sample(nrow(d)), simplify = FALSE)          # identical to posthoc_roworder.R
set.seed(20260925L + 900L)
NB <- as.integer(Sys.getenv("GLOW_NB", "100"))
perms2 <- replicate(300, sample(nrow(d)), simplify = FALSE)
OUT <- file.path(OUTD, "glow_roworder.csv"); dk <- done_keys(OUT)
run1 <- function(key, part, dd, mn, seed, order) {
  if (key %in% dk) return(invisible())
  fit <- glm(if (mn == "M1") f1 else f2, data = dd, family = binomial())
  set.seed(seed); r <- dg_new(fit, B = 199L)
  append_row(OUT, data.frame(key = key, part = part, model = mn, order = order, seed = seed,
                             axes = paste(r$axes, collapse = "|"), S = r$statistic, p = r$p.value))
}
for (mn in strsplit(Sys.getenv("GLOW_MODELS", "M1,M2"), ",")[[1]]) {
  for (k in 0:30) run1(paste("A", mn, k, sep = "_"), "A", if (k == 0) d else d[perms[[k]], ], mn,
                       20260925L + 1L, k)
  for (j in 1:NB) run1(paste("Bseed", mn, j, sep = "_"), "B_fileorder_seeds", d, mn, 30000L + j, 0L)
  for (j in 1:NB) run1(paste("Bperm", mn, j, sep = "_"), "B_randomorder_fixedseed", d[perms2[[j]], ],
                        mn, 20260925L + 1L, j)
  cat(mn, "done\n"); flush.console()
}
cat("GLOW ROWORDER COMPLETE\n")
