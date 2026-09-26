"""clean_jsonl.py -- after a hard power-off, drop any line of the jsonl files that does not parse
(a truncated last line), and de-duplicate net points by t.  Run: python clean_jsonl.py"""
import os, glob, json
HERE = os.path.dirname(os.path.abspath(__file__))
for path in sorted(glob.glob(os.path.join(HERE, "*.jsonl"))):
    raw = open(path, encoding="utf-8", errors="replace").read().split("\n")
    good, bad, seen = [], 0, set()
    for line in raw:
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            bad += 1; continue
        key = r.get("t", None) if os.path.basename(path).startswith("net_") else None
        if key is not None:
            if key in seen: bad += 1; continue
            seen.add(key)
        good.append(json.dumps(r))
    with open(path, "w") as fh:
        fh.write("\n".join(good) + ("\n" if good else ""))
    print("%-28s kept %4d  dropped %d" % (os.path.basename(path), len(good), bad))
