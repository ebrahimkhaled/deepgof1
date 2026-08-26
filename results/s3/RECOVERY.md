# S3 recovery notes (if the run dies overnight)

The driver is IDEMPOTENT. Just re-run it:

    bash code/run_s3.sh

- A (design, n, seed) block with `summary.csv` is skipped entirely.
- A block with `blk_meta.csv` + `blk_maps.bin` but no summary is REUSED (not regenerated) --
  generation is the expensive part, scoring is cheap.
- A group is only scored once all four designs are present, because the pre-declared rule is
  trained on D1 and applied to D1-D4 together (speckit/21 section 5).
- `.bin` files are deleted after a group scores, so peak disk stays near 600 MB.

If a group has some blocks scored and others not, delete that group's `summary.csv` files and
re-run -- the block data will be reused and the group rescored consistently.

Verdict: `python code/s3_verdict.py` (to be written) or read `results/s3/*/summary.csv`.
Criteria are frozen in `speckit/21_S3_PREDECLARATION.md` section 6.
