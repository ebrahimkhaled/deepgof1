"""run_all.py -- driver: 01 null quantiles, then 02 net / 03 tail / 04 truth per direction in parallel
(at most MAXP processes while an R simulation holds the box: 4 while tasklist shows >= 3 Rscript
processes, 5 afterwards), then 05 certify.  Every stage is resumable; rerun after a crash.
Run: python -u run_all.py
"""
import os, sys, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from smooth_core import ORDER, rscript_count

PY = sys.executable
def run(script, *args):
    return subprocess.Popen([PY, "-u", os.path.join(HERE, script), *args], cwd=HERE,
                            stdout=open(os.path.join(HERE, "log_%s_%s.txt" % (script[:2], "_".join(args) or "all")), "a"),
                            stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS="2"))

def pool(jobs, log):
    running = []
    jobs = list(jobs); t0 = time.time()
    while jobs or running:
        maxp = int(os.environ.get("MAXP", "8"))
        running = [p for p in running if p.poll() is None]
        while jobs and len(running) < maxp:
            script, arg = jobs.pop(0)
            log("  start %s %s (%d running, Rscript=%d)" % (script, arg, len(running) + 1, rscript_count()))
            running.append(run(script, arg))
        time.sleep(10)
    log("  stage done in %.0f s" % (time.time() - t0))

if __name__ == "__main__":
    log = lambda s: print(time.strftime("%H:%M:%S ") + s, flush=True)
    log("stage 01: null quantiles")
    if not os.path.exists(os.path.join(HERE, "NULL_UCB.json")):
        p = run("01_null_quantiles.py"); p.wait()
    log("stage 02+03+04: net, tail and truth per direction (independent jobs, one pool)")
    from smooth_core import read_jsonl
    import importlib
    design = importlib.import_module("02_mc_net").design
    jobs = []
    for d in ORDER:
        if len(read_jsonl(os.path.join(HERE, "net_%s.jsonl" % d))) < len(design(d)):
            jobs.append(("02_mc_net.py", d))
    for d in ORDER:
        s = os.path.join(HERE, "tail_%s_summary.json" % d)
        if not (os.path.exists(s) and __import__("json").load(open(s))["closed"]):
            jobs.append(("03_tail.py", d))
    for d in ORDER:
        if len(read_jsonl(os.path.join(HERE, "truth_%s.jsonl" % d))) < 5:
            jobs.append(("04_truth.py", d))
    log("  jobs still needed: %s" % jobs)
    pool(jobs, log)
    log("stage 05: certify")
    p = subprocess.run([PY, "-u", os.path.join(HERE, "05_certify.py")], cwd=HERE, capture_output=True, text=True)
    print(p.stdout); print(p.stderr)
    log("all done")
