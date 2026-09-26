## VERIFY.R -- the cross-language gate. Run after ANY retraining/export, and once on install.
##
## Scores the shipped verify maps (inst/extdata/verify_input.csv) with the base-R forward pass and
## compares against the PyTorch scores written by 03_export_to_R.py. PASS = max |difference| < 1e-6.
## DeepGOF-1's shipped pair agrees to 4.6e-11. If this fails, the row/column-major port trap is the first suspect (see the note
## in 03_export_to_R.py) -- our first port agreed to 7.5e-3, which LOOKS fine and is not.
##
## Usage:  Rscript training/VERIFY.R        (from the release root)

## robust script-location: works under Rscript from ANY working directory, and under source()
args <- commandArgs(trailingOnly = FALSE)
fa   <- sub("^--file=", "", args[grepl("^--file=", args)])
here <- if (length(fa)) dirname(normalizePath(fa[1])) else
        tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "training")
root <- normalizePath(file.path(here, ".."))
source(file.path(root, "R", "deepgof1_forward.R"))

ed <- file.path(root, "inst", "extdata")
W  <- load_weights(file.path(ed, "deepgof1_weights.csv"))
X  <- as.matrix(utils::read.csv(file.path(ed, "verify_input.csv"),  header = FALSE))
Sp <- as.numeric(readLines(file.path(ed, "verify_output.csv")))
## verify_input.csv rows are ALREADY-STANDARDIZED maps (z-space) -- score them directly.
Sr <- apply(X, 1, function(r) deepgof_score(r, W))
dmax <- max(abs(Sr - Sp))
cat(sprintf("base R vs PyTorch on %d verify maps: max |difference| = %.3e  ->  %s\n",
            length(Sp), dmax, if (dmax < 1e-6) "PASS" else "FAIL"))

t0 <- Sys.time(); for (i in 1:200) deepgof_score(X[1, ], W)
cat(sprintf("200 forward passes (one deployed test): %.2f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))
if (dmax >= 1e-6) stop("VERIFY FAILED -- do not deploy these weights")
