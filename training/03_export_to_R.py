"""03_export_to_R.py -- export a trained CNN + scaler to the base-R deployment format.

Writes four files consumed by R/deepgof1_forward.R:
  weights.csv   name,shape,values  -- one row per parameter tensor, values space-separated,
                shape as "AxBxC" in PyTorch's own dimension order
  scaler.csv    two headerless rows: mu then sd (36 columns each)
  verify_input.csv / verify_output.csv -- 20 held-out maps and their PyTorch scores, the
                cross-language test bench that VERIFY.R must reproduce to < 1e-6

THE PORT TRAP (learned the hard way; do not remove this note). PyTorch stores tensors ROW-major
and x.view(-1,1,K,K) fills row-major; R's array() fills COLUMN-major. The R loader therefore
rebuilds each tensor with aperm(array(v, rev(shape)), n:1), and the R scorer rebuilds the map
with matrix(v, K, K, byrow=TRUE). Exporting in any other layout produces a port that is
plausibly wrong: our first attempt agreed to 7.5e-3 -- close enough to look right, wrong enough
to change p-values. The verify bench exists because of that.

Usage (after 02_train.py has produced a trained `net`, `mu`, `sd`):
    from importlib import import_module
    export_to_R(net, mu, sd, test_maps, out_dir="../inst/extdata", prefix="deepgof1_")
"""
import os
import numpy as np
import torch


def export_to_R(net, mu, sd, test_maps, out_dir=".", prefix="deepgof1_", n_test=20):
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    for name, t in net.state_dict().items():
        a = t.detach().cpu().numpy()
        shape = "x".join(str(d) for d in a.shape) if a.ndim else "1"
        vals = " ".join(repr(float(v)) for v in a.reshape(-1))   # row-major, full precision
        rows.append((name, shape, vals))
    with open(os.path.join(out_dir, prefix + "weights.csv"), "w", encoding="utf-8") as f:
        f.write("name,shape,values\n")
        for name, shape, vals in rows:
            f.write(f"{name},{shape},{vals}\n")

    with open(os.path.join(out_dir, prefix + "scaler.csv"), "w", encoding="utf-8") as f:
        f.write(",".join(repr(float(v)) for v in np.asarray(mu).reshape(-1)) + "\n")
        f.write(",".join(repr(float(v)) for v in np.asarray(sd).reshape(-1)) + "\n")

    # the verify bench is stored in Z-SPACE (standardized maps): VERIFY.R scores rows directly
    X = np.asarray(test_maps, dtype=np.float64)[:n_test]
    Z = (X - np.asarray(mu)) / np.asarray(sd)
    net.eval()
    with torch.no_grad():
        S = net(torch.from_numpy(Z.astype(np.float32))).cpu().numpy().reshape(-1)
    np.savetxt(os.path.join(out_dir, "verify_input.csv"), Z, delimiter=",")
    np.savetxt(os.path.join(out_dir, "verify_output.csv"), S, fmt="%.10f")
    print(f"exported {len(rows)} tensors + scaler + {len(S)} verify vectors -> {out_dir}")
    print("now run training/VERIFY.R -- it must print max |difference| < 1e-6")


if __name__ == "__main__":
    print(__doc__)
