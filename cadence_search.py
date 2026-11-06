#!/usr/bin/env python3
"""cadence_search.py — windowed ON/OFF search across a whole BL cadence.

WHY WINDOWED
------------
A fine-frequency BL scan is the entire 750 MHz L-band at 2.836 Hz: 264,503,296
channels x 16 integrations = 4.2 billion floats, ~17 GB as one array. seti_io
rightly REFUSES to open that in one piece (max_chans guard). Even with 128 GB of
RAM it is the wrong shape to work in: the search is local in frequency, so the
natural unit is ONE COARSE CHANNEL — 2.93 MHz, 1,033,215 channels, the same
dimensions as the Voyager tutorial file the GPU path does in ~4.4 s.

The payoff of that choice: all SIX scans of one window are only ~400 MB
together, so every window is compared ON-vs-OFF *in memory, at the same time*,
which is what the cadence test actually needs. Big RAM is what lets the whole
cadence be resident per window; the GPU is what makes each window fast.

    python cadence_search.py --target GJ699                    # water hole
    python cadence_search.py --target GJ699 --f-start 1126 --f-stop 1876  # all
    python cadence_search.py --target GJ699 --recipe narrowband_drift

THE TEST
--------
A signal from the target appears in the ON scans and NOT in the OFF scans.
Anything present in an OFF scan is local. This is how BLC-1 died (Sheikh et al.
2021, Nature Astronomy 5, 1153) and zero survivors is the normal, honest result.
A survivor here is NOT a detection — it is a thing that has not yet been
explained, which is a different and much weaker claim.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
sys.path.insert(0, str(HERE))

COARSE_MHZ = 2.9296875


def scans_for(target, product="0000"):
    fs = sorted(DATA.glob(f"*{target}*.gpuspec.{product}.h5"))
    return [(("OFF" if "_OFF_" in p.name.upper() else "ON"), p) for p in fs]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="GJ699")
    ap.add_argument("--product", default="0000")
    ap.add_argument("--recipe", default="all")
    # default band = the "water hole", 1420 (HI) to 1666 (OH): the classic SETI
    # window, quiet and the one every survey looks at first.
    ap.add_argument("--f-start", type=float, default=1400.0)
    ap.add_argument("--f-stop", type=float, default=1500.0)
    ap.add_argument("--window-mhz", type=float, default=COARSE_MHZ)
    ap.add_argument("--tol-hz", type=float, default=600.0,
                    help="ON/OFF frequency match tolerance")
    ap.add_argument("--gpu", action="store_true", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if a.gpu is None:
        try:
            import cupy
            cupy.zeros(1)
            a.gpu = True
        except Exception:
            a.gpu = False
    if a.gpu:
        os.environ["SETI_COOKBOOK_GPU"] = "1"

    import seti_io
    import recipe_api as R

    scans = scans_for(a.target, a.product)
    if len(scans) < 2:
        print(f"need >=2 scans, found {len(scans)}", file=sys.stderr)
        return 1
    recs, _ = R.discover()
    if a.recipe != "all":
        # comma-separated list. Profiled 2026-07-30 on one 2.93 MHz window:
        # narrowband_drift 56.1 s and dispersion_sweep 53.8 s dominate a 114 s
        # total, while doppler_coherence on the GPU does a comparable de-drift
        # search in 3.25 s. Being able to pick the set is the difference between
        # a 20-minute search and a 6.5-hour one.
        want = {x.strip() for x in a.recipe.split(",") if x.strip()}
        recs = [r for r in recs if r.info().get("name") in want]
        missing = want - {r.info().get("name") for r in recs}
        if missing:
            print(f"WARNING: unknown recipe(s) {sorted(missing)}", flush=True)

    edges = np.arange(a.f_start, a.f_stop, a.window_mhz)
    print(f"{a.target}: {len(scans)} scans, {len(edges)} window(s) of "
          f"{a.window_mhz:.4f} MHz over {a.f_start}-{a.f_stop} MHz"
          f"{'  [GPU]' if a.gpu else '  [CPU]'}")
    print(f"pattern: {' '.join(k for k, _ in scans)}")
    print(f"recipes: {', '.join(r.info()['name'] for r in recs)}\n")

    survivors, per_window, t0 = [], [], time.time()
    for wi, f0 in enumerate(edges):
        f1 = min(f0 + a.window_mhz, a.f_stop)
        on_hits, off_hits, errs = [], [], []
        for kind, path in scans:
            try:
                sp = seti_io.open_any(str(path), f_start=float(f0), f_stop=float(f1))
            except Exception as e:
                errs.append(f"{path.name[:30]}: open {type(e).__name__}: {e}")
                continue
            for r in recs:
                nm = r.info()["name"]
                try:
                    # Recipe.run returns (candidates, elapsed) and has ALREADY
                    # called as_candidates internally — wrapping it again is
                    # what produced "Candidate() argument after ** must be a
                    # mapping, not list". Unpack; do not re-wrap.
                    cs, _dt = r.run(
                        sp, {"gpu": True} if (a.gpu and nm == "doppler_coherence")
                        else None)
                    cs = R.explain(cs, sp)
                except Exception as e:
                    # never silent: a recipe that dies here would otherwise look
                    # exactly like a recipe that found nothing
                    errs.append(f"{nm}: {type(e).__name__}: {e}")
                    continue
                for c in cs:
                    c.provenance["_scan"] = path.name
                (on_hits if kind == "ON" else off_hits).extend(cs)

        # ── THE CADENCE CRITERION, PROPERLY ────────────────────────────────
        # A real signal from the target must be in EVERY ON scan and in NO OFF
        # scan. The first version of this only asked "in an ON, not in an OFF",
        # which a one-off blip in a single ON passes — and that is exactly what
        # intermittent local interference looks like. On the 1400-1500 MHz run
        # that laxity produced 51 "survivors", 41 of them piled into 1442-1452
        # MHz, which is an RFI cluster wearing a candidate's clothes.
        # turboSETI's find_event requires presence in all ON scans for the same
        # reason; require it here too.
        off_f = np.array([c.freq_mhz for c in off_hits]) if off_hits else np.array([])
        n_on_scans = sum(1 for k, _ in scans if k == "ON")
        keep = []
        for c in on_hits:
            if "unexplained" not in (c.verdict or ""):
                continue
            if off_f.size and np.min(np.abs(off_f - c.freq_mhz)) * 1e6 <= a.tol_hz:
                continue                      # seen in an OFF scan -> local
            # how many DISTINCT ON scans contain this frequency?
            seen = {oc.provenance.get("_scan") for oc in on_hits
                    if abs(oc.freq_mhz - c.freq_mhz) * 1e6 <= a.tol_hz}
            if len(seen) < n_on_scans:
                continue                      # not in every ON -> intermittent
            keep.append(c)
        survivors += keep
        per_window.append(dict(f0=float(f0), f1=float(f1), n_on=len(on_hits),
                               n_off=len(off_hits), n_survive=len(keep),
                               errors=errs[:4]))
        if wi % 5 == 0 or keep or errs:
            el = time.time() - t0
            eta = el / (wi + 1) * (len(edges) - wi - 1)
            print(f"  [{wi+1:>3}/{len(edges)}] {f0:8.3f}-{f1:8.3f} MHz  "
                  f"ON {len(on_hits):>4} / OFF {len(off_hits):>4}  "
                  f"survive {len(keep):>2}   eta {eta/60:5.1f} min", flush=True)
            for e in errs[:2]:
                print(f"        ERROR {e}", flush=True)

    el = time.time() - t0
    print(f"\n=== {a.target} {a.f_start}-{a.f_stop} MHz done in {el/60:.1f} min ===")
    print(f"ON-only survivors: {len(survivors)}")
    for c in sorted(survivors, key=lambda x: -x.score)[:25]:
        print(f"  {c.freq_mhz:12.6f} MHz  score {c.score:9.1f}  "
              f"drift {c.drift_hz_s:+7.3f}  {c.recipe:<18} {c.label[:40]}")
    if not survivors:
        print("  (none — the normal, honest result)")

    out = Path(a.out) if a.out else HERE / f"cadence_{a.target}_{int(a.f_start)}_{int(a.f_stop)}.json"
    out.write_text(json.dumps(dict(
        target=a.target, f_start=a.f_start, f_stop=a.f_stop,
        window_mhz=a.window_mhz, scans=[(k, p.name) for k, p in scans],
        elapsed_s=el, windows=per_window,
        survivors=[dict(freq_mhz=c.freq_mhz, score=c.score, drift_hz_s=c.drift_hz_s,
                        recipe=c.recipe, label=c.label, verdict=c.verdict)
                   for c in survivors]), indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
