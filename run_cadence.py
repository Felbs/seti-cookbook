#!/usr/bin/env python3
"""run_cadence.py — run the panel across a whole ON/OFF cadence, in parallel.

Built for a big machine: the six fine-frequency scans of a BL cadence are ~12.8
GB each, and a box with 128 GB of RAM can hold several of them resident at once
instead of streaming each one twice. Scans are independent until the very last
step, so they are processed CONCURRENTLY and only the candidate lists — a few
kilobytes — come back to be compared.

    python run_cadence.py --target GJ699                 # auto-find the cadence
    python run_cadence.py --target GJ699 --workers 3     # cap memory
    python run_cadence.py --target GJ699 --recipe doppler_coherence

CPU WORKERS OR GPU? DO THE ARITHMETIC, NOT THE INTUITION
--------------------------------------------------------
"Six scans, sixty-four cores, fan them out" is the obvious move and it is the
WRONG one for the fine-frequency products. Each scan needs ~2.2x its 12.8 GB
resident, so 111 GB of free RAM supports about TWO concurrent workers, not
sixty-four — memory, not cores, is the binding constraint:

    CPU, 2 workers  : 6 scans / 2   = 3 rounds x T_cpu   = 3.0 T_cpu
    GPU, sequential : 6 x T_cpu/21                       = 0.29 T_cpu

so the single GPU beats the fan-out by roughly 10x, and --gpu is the DEFAULT
whenever a device is present. CPU workers are the right answer only for the
small products (0001/0002), where many scans fit in RAM at once and the GPU
would be starved feeding them.

The 21x figure is measured, not assumed (Voyager, 92.1 s CPU vs 4.4 s GPU,
identical result). Note also that the CPU and GPU take DIFFERENT code paths on
purpose: an index-gather beats np.roll on a GPU and loses badly to it on a CPU.

MEMORY GUARD: each worker holds roughly one scan, so peak RSS is about
workers x filesize. The default worker count is chosen from free RAM rather than
core count, because running out of memory on a 77 GB cadence is a much worse
failure than being one core short.

THE POINT OF A CADENCE: a signal from the target appears in the ON scans and NOT
in the OFF scans. Anything present in an OFF scan is local — this is exactly how
BLC-1 died (Sheikh et al. 2021, Nature Astronomy 5, 1153). Zero survivors is the
normal, honest result.
"""
import argparse
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
sys.path.insert(0, str(HERE))


def _scan_one(args):
    """Worker: run recipes on ONE scan. Returns small dicts, never arrays."""
    path, recipe, use_gpu = args
    sys.path.insert(0, str(HERE))
    import seti_io
    import recipe_api as R
    t0 = time.time()
    try:
        spec = seti_io.open_any(path)
    except Exception as e:
        return dict(path=path, error=f"open failed: {type(e).__name__}: {e}")
    recs, _ = R.discover()
    if recipe and recipe != "all":
        recs = [r for r in recs if r.info().get("name") == recipe]
    cands = []
    errors = []
    for r in recs:
        name = r.info().get("name")
        try:
            # Recipe.run returns (candidates, elapsed) and already applied
            # as_candidates — re-wrapping raises "argument after ** must be a
            # mapping, not list".
            params = {"gpu": True} if (use_gpu and name == "doppler_coherence") else None
            cs, _dt = r.run(spec, params)
            cands += cs
        except Exception as e:
            # Do NOT swallow these silently — a recipe that dies on real data
            # while the summary still prints "0 candidates" is how you conclude
            # "nothing there" from a pipeline that never ran.
            errors.append(f"{name}: {type(e).__name__}: {e}")
    cands = R.explain(cands, spec)
    return dict(
        path=path, name=Path(path).name,
        kind="OFF" if "_OFF_" in Path(path).name.upper() else "ON",
        f_lo=float(spec.f_lo), f_hi=float(spec.f_hi),
        res_hz=float(abs(spec.res_hz)), duration_s=float(spec.duration_s),
        secs=round(time.time() - t0, 1), errors=errors,
        candidates=[dict(freq_mhz=c.freq_mhz, score=c.score,
                         drift_hz_s=c.drift_hz_s, recipe=c.recipe,
                         kind=c.kind, label=c.label, verdict=c.verdict)
                    for c in cands])


def _auto_workers(files, cap=None):
    """Pick worker count from FREE RAM, not core count — see MEMORY GUARD."""
    try:
        import psutil
        free = psutil.virtual_memory().available
    except Exception:
        free = 8 << 30
    biggest = max((Path(f).stat().st_size for f in files), default=1 << 30)
    # ~2.2x the file size per worker: the array, plus room for the bandpass copy
    n = max(1, int(free * 0.75 // (biggest * 2.2)))
    n = min(n, len(files), os.cpu_count() or 4)
    return min(n, cap) if cap else n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="GJ699")
    ap.add_argument("--product", default="0000")
    ap.add_argument("--recipe", default="all")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--gpu", action="store_true", default=None,
                    help="force the single-worker GPU path (default: auto — GPU "
                         "whenever a device exists, because it wins ~10x on the "
                         "fine products; see the module docstring)")
    ap.add_argument("--cpu", dest="gpu", action="store_false",
                    help="force CPU fan-out even if a GPU is present")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    files = sorted(str(p) for p in DATA.glob(f"*{a.target}*.gpuspec.{a.product}.h5"))
    if not files:
        print(f"no {a.target} product-{a.product} files in {DATA}", file=sys.stderr)
        return 1
    if a.gpu is None:                     # auto-detect
        try:
            import cupy
            cupy.zeros(1)
            a.gpu = True
        except Exception:
            a.gpu = False
    if a.gpu:
        os.environ["SETI_COOKBOOK_GPU"] = "1"
    nw = 1 if a.gpu else (a.workers or _auto_workers(files))
    total_gb = sum(Path(f).stat().st_size for f in files) / 1e9
    print(f"{a.target}: {len(files)} scan(s), {total_gb:.1f} GB, "
          f"{nw} worker(s){' [GPU]' if a.gpu else ''}")
    for f in files:
        print(f"  {'OFF' if '_OFF_' in Path(f).name.upper() else 'ON ':<3} "
              f"{Path(f).stat().st_size/1e9:5.2f} GB  {Path(f).name[:58]}")

    t0 = time.time()
    jobs = [(f, a.recipe, a.gpu) for f in files]
    if nw == 1:
        results = [_scan_one(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=nw) as ex:
            results = list(ex.map(_scan_one, jobs))
    print(f"\nall scans done in {time.time()-t0:.0f} s")

    for r in results:
        if r.get("error"):
            print(f"  {r['path']}: {r['error']}")
            continue
        un = [c for c in r["candidates"] if "unexplained" in (c["verdict"] or "")]
        print(f"  {r['kind']:<3} {r['name'][:52]:<54} {r['secs']:6.1f}s  "
              f"{len(r['candidates']):>4} cand  {len(un):>3} unexplained")
        for e in r.get("errors", []):
            print(f"        RECIPE ERROR {e}")

    # the cadence test itself
    import recipe_api as R
    per = []
    for r in results:
        if r.get("error"):
            continue
        cs = [R.Candidate(freq_mhz=c["freq_mhz"], score=c["score"],
                          drift_hz_s=c["drift_hz_s"], kind=c["kind"],
                          label=c["label"]) for c in r["candidates"]]
        per.append((r["kind"], cs))
    try:
        verdicts = R.cadence_verify([c for _k, c in per],
                                    pattern=[k for k, _c in per])
    except Exception as e:
        verdicts = {"error": f"{type(e).__name__}: {e}"}
    print("\n=== CADENCE VERDICT (ON-only survives; anything in an OFF is local) ===")
    print(json.dumps(verdicts, indent=2, default=str)[:2500])

    out = Path(a.out) if a.out else HERE / f"cadence_{a.target}.json"
    out.write_text(json.dumps(dict(target=a.target, results=results,
                                   cadence=verdicts), indent=2, default=str),
                   encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
