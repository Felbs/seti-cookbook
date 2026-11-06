#!/usr/bin/env python3
"""fastsearch_recovery.py - fastsearch vs turboSETI on the SAME planted real data.

Takes one real coarse channel (a hi_starframe .npz slice), plants a grid of
Doppler-smeared setigen tones (inject.plant), and searches the same array with
turboSETI (max_drift 4, min_drift 0, snr 10 - the GJ699 cadence settings) and
with fastsearch. The two are compared at the SAME false-alarm level: on the
UNPLANTED data, fastsearch's threshold is set so that it reports no more hits
than turboSETI does (the data's own lines), and plants are kept away from
every null hit of either search.

Scoring is inject.match (one-to-one along the track) + inject.drift_ok.
Strength is TRUE radiometer SNR (inject.radiometer_noise).

The npz slice of 1,048,576 channels does not line up with this file's coarse
channels (1,033,216 = 2.9296875 MHz each, as blimpy computes), so only the
part inside ONE coarse channel is used, and its DC spike (at the real coarse
centre) is repaired the same way for both searches.

    python fastsearch_recovery.py data/hi_..._GJ699_0003.gpuspec.0000.npz --seed 1 --json r1.json
    python fastsearch_recovery.py selftest     # synthetic, fastsearch only, no turboSETI needed
"""
import argparse
import json
import shutil
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fastsearch as F  # noqa: E402
import inject  # noqa: E402

DRIFTS = [0.0, 0.1, -0.1, 0.5, -0.5, 1.0, -1.0, 2.0, -2.0, 3.5, -3.5]
SNRS = [10, 15, 20, 25, 30, 40, 50, 70, 100, 140, 200, 300, 400, 600, 800, 1200, 1600, 2400]
GBT_NFPC = 1033216


def coarse_segment(data, m, nfpc=GBT_NFPC):
    """Largest part of the slice inside one coarse channel; returns (seg, lo, dc_in_seg)."""
    ch0 = int(m.get("ch0", 0))
    nf = data.shape[1]
    off = ch0 % nfpc
    first = nfpc - off                       # channels of the first coarse channel in the slice
    if first >= nf - first:
        lo, hi, dc = 0, min(first, nf), nfpc // 2 - off
    else:
        lo, hi, dc = first, nf, nfpc // 2 - off + nfpc - first
    seg = np.ascontiguousarray(data[:, lo:hi])
    if 0 <= dc < seg.shape[1]:               # repair DC like turboSETI's blank_dc does on a real file
        seg[:, dc] = 0.5 * (seg[:, dc - 1] + seg[:, dc + 1])
    else:
        dc = None
    return seg, lo, dc


def ts_search(tmp, name, data, m, snr):
    import recovery
    fn = tmp / f"{name}.h5"
    recovery.write_h5(fn, data, m)
    return inject.read_dat(recovery.search(fn, 4.0, snr), m["fch1"], m["foff"])


def fs_search(data, m, thr, dc):
    hits, info = F.search_slab(data, m["foff"] * 1e6, m["tsamp"], 4.0, thr,
                               dc_idx=() if dc is None else (dc,))
    return hits, info


def cell_table(plants, why):
    t = {}
    for i, p in enumerate(plants):
        key = f"{p['drift']:+.2f}|{p['snr']:g}"
        c = t.setdefault(key, dict(found=0, n=0, statuses=[]))
        c["n"] += 1
        c["found"] += why[i]["status"] == "found"
        c["statuses"].append(why[i]["status"])
    return t


def snr90(tables, drift_abs, snrs=SNRS):
    """True SNR where >= 90 % of plants at |drift| are found (pooled over signs and
    all tables); log-interpolated between grid points; None if never reached."""
    frac = []
    for s in snrs:
        f = n = 0
        s = int(s) if float(s).is_integer() else s
        for t in tables:
            for sign in ((+1, -1) if drift_abs > 0 else (+1,)):
                c = t.get(f"{sign * drift_abs:+.2f}|{s:g}")
                if c:
                    f += c["found"]
                    n += c["n"]
        frac.append(f / n if n else np.nan)
    frac = np.array(frac)
    for i, fr in enumerate(frac):
        if fr >= 0.9:
            if i == 0:
                return snrs[0], frac
            a, b = frac[i - 1], fr
            x = np.log(snrs[i - 1]) + (0.9 - a) / (b - a) * (np.log(snrs[i]) - np.log(snrs[i - 1]))
            return float(np.exp(x)), frac
    return None, frac


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--per-cell", type=int, default=4)
    ap.add_argument("--ts-snr", type=float, default=10.0)
    ap.add_argument("--fs-snr", type=float, help="fastsearch threshold (default: matched to turboSETI's null count)")
    ap.add_argument("--json")
    ap.add_argument("--no-turboseti", action="store_true")
    ap.add_argument("--snrs", help="comma list of TRUE SNRs (default: the standard grid)")
    a = ap.parse_args()
    global SNRS
    if a.snrs:
        SNRS = [float(v) for v in a.snrs.split(",")]
    import recovery
    data, m = recovery.load(a.npz)
    seg, lo, dc = coarse_segment(data, m)
    del data
    m = dict(m)
    m["fch1"] = float(m["fch1"]) + float(m["foff"]) * lo
    m["foff"], m["tsamp"] = float(m["foff"]), float(m["tsamp"])
    nt, nc = seg.shape
    foff_hz, ts = m["foff"] * 1e6, m["tsamp"]
    print(f"{Path(a.npz).name}: coarse-channel segment {nc} ch from slice ch {lo}; DC at {dc}")
    tmp = Path(tempfile.mkdtemp(prefix="fsrec_"))
    out = dict(file=Path(a.npz).name, seed=a.seed, per_cell=a.per_cell, nchan=nc, dc=dc, snrs=SNRS)
    try:
        # 1. null probes
        t0 = time.time()
        fs0, info = fs_search(seg, m, 5.0, dc)
        t_fs = time.time() - t0
        fs_snrs = np.sort([h["snr"] for h in fs0])[::-1]
        if not a.no_turboseti:
            t0 = time.time()
            ts0 = ts_search(tmp, "null", seg, m, a.ts_snr)
            t_ts = time.time() - t0
            n_ts = len(ts0)
        else:
            ts0, n_ts, t_ts = [], 0, 0
        if a.fs_snr:
            thr = a.fs_snr
        else:   # smallest threshold (>= 6) at which fastsearch's null count <= turboSETI's
            thr = 6.0
            while (fs_snrs >= thr).sum() > n_ts:
                thr += 0.25
        fs_null = [h for h in fs0 if h["snr"] >= thr]
        print(f"null probe: turboSETI {n_ts} hits at its SNR {a.ts_snr:g} ({t_ts:.0f} s); fastsearch "
              f"{len(fs_null)} hits at true SNR {thr:g} ({t_fs:.0f} s, +-{info['drift_searched']:.3f} Hz/s); "
              f"fastsearch counts at 6/7/8/10: "
              + "/".join(str(int((fs_snrs >= v).sum())) for v in (6, 7, 8, 10)))
        out.update(ts_null=n_ts, fs_null=len(fs_null), fs_thr=thr, t_ts_null=t_ts, t_fs_null=t_fs,
                   fs_null_counts={str(v): int((fs_snrs >= v).sum()) for v in (6, 7, 8, 10)})
        T = (nt - 1) * ts
        radius = int(np.ceil(4.0 * T / abs(foff_hz)))
        avoid = [h["ch"] for h in ts0] + [h["ch"] for h in fs0 if h["snr"] >= 6]
        cells = [(d, s) for d in DRIFTS for s in SNRS]
        n = len(cells) * a.per_cell
        starts = recovery.layout(nc, n, avoid, guard=2 * radius + 50, seed=a.seed)
        rng = np.random.default_rng(a.seed)
        order = rng.permutation(n)
        plants = [dict(ch=float(starts[i]), drift=cells[order[i] % len(cells)][0],
                       snr=cells[order[i] % len(cells)][1]) for i in range(n)]
        d = seg.copy()
        inject.plant(d, foff_hz, ts, plants)
        res = {}
        # 2. fastsearch
        t0 = time.time()
        hf, _ = fs_search(d, m, thr, dc)
        out["t_fs"] = time.time() - t0
        hf = [h for h in hf if all(abs(h["ch"] - o["ch"]) > 3 for o in fs_null)]
        mf = inject.match(plants, hf, foff_hz, ts, nt)
        wf = inject.explain_misses(plants, hf, mf, foff_hz, ts, nt, 4.0)
        res["fastsearch"] = cell_table(plants, wf)
        zr = [hf[j]["snr_power"] / p["snr"] for p, j in zip(plants, mf) if j is not None and p["drift"] == 0]
        out["fs_zero_drift_power_ratio"] = float(np.median(zr)) if zr else None
        # 3. turboSETI
        if not a.no_turboseti:
            t0 = time.time()
            ht = ts_search(tmp, "planted", d, m, a.ts_snr)
            out["t_ts"] = time.time() - t0
            ht = [h for h in ht if all(abs(h["ch"] - o["ch"]) > 3 for o in ts0)]
            mt = inject.match(plants, ht, foff_hz, ts, nt)
            wt = inject.explain_misses(plants, ht, mt, foff_hz, ts, nt, 4.0)
            res["turboseti"] = cell_table(plants, wt)
        out["tables"] = res
        print(f"planted {n} smeared tones ({a.per_cell}/cell, seed {a.seed}); fastsearch "
              f"{out['t_fs']:.0f} s, turboSETI {out.get('t_ts', 0):.0f} s")
        for tool, tab in res.items():
            print(f"\n{tool}: fraction found with correct drift (rows drift, cols TRUE SNR)")
            print("drift  " + "".join(f"{s:>6g}" for s in SNRS))
            for dr in DRIFTS:
                print(f"{dr:+5.2f}  " + "".join(
                    f"{tab[f'{dr:+.2f}|{s:g}']['found'] / tab[f'{dr:+.2f}|{s:g}']['n']:6.2f}" for s in SNRS))
            st = {}
            for c in tab.values():
                for s in c["statuses"]:
                    st[s] = st.get(s, 0) + 1
            print("statuses:", st)
        if a.json:
            Path(a.json).write_text(json.dumps(out, indent=1))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def summarise(paths):
    """Pool several runs' JSON: 90 %-recovery TRUE SNR per |drift| for each tool."""
    runs = []
    for p in paths:
        j = json.loads(Path(p).read_text())
        runs.extend(j["runs"] if "runs" in j else [j])     # one run, or a combined results file
    print(f"{len(runs)} runs: seeds {[r['seed'] for r in runs]}; fastsearch thresholds "
          f"{[r['fs_thr'] for r in runs]}; null hits turboSETI {[r['ts_null'] for r in runs]} "
          f"fastsearch {[r['fs_null'] for r in runs]}")
    tools = sorted({t for r in runs for t in r["tables"]})
    grid = sorted({float(k.split("|")[1]) for r in runs for t in r["tables"].values() for k in t})
    print("|drift| " + "".join(f"{t:>14s}" for t in tools) + "   ratio")
    rows = {}
    for dr in (0.0, 0.1, 0.5, 1.0, 2.0, 3.5):
        vals = [snr90([r["tables"][t] for r in runs if t in r["tables"]], dr, snrs=grid)[0] for t in tools]
        rows[dr] = vals
        ratio = (vals[1] / vals[0]) if len(vals) == 2 and vals[0] and vals[1] else None
        print(f"{dr:6.2f} " + "".join(f"{(f'{v:.0f}' if v else '>max'):>14s}" for v in vals)
              + (f"   {ratio:.1f}x" if ratio else ""))
    return rows


def selftest():
    """Synthetic slab: fastsearch finds a coarse grid of plants; noise-only gives no hit
    at threshold 8; snr90 interpolation is right on a known table."""
    rng = np.random.default_rng(5)
    nt, nf, foff, ts = 16, 65536, -2.835503418452676, 17.986224128
    d = F._synthetic(rng, nt, nf, dc=False)
    m = dict(foff=foff / 1e6, tsamp=ts, fch1=1500.0)
    h0, _ = fs_search(d, m, 8.0, None)
    pl = [dict(ch=float(c), drift=dr, snr=s) for c, (dr, s) in
          zip(range(3000, 63000, 6000), [(0.0, 60), (0.5, 60), (-1, 100), (2, 150), (-3.5, 250)] * 2)]
    inject.plant(d, foff, ts, pl)
    h, _ = fs_search(d, m, 8.0, None)
    mm = inject.match(pl, h, foff, ts, nt)
    why = inject.explain_misses(pl, h, mm, foff, ts, nt, 4.0)
    tab = {"+0.10|10": dict(found=0, n=10), "+0.10|20": dict(found=8, n=10), "+0.10|40": dict(found=10, n=10),
           "-0.10|10": dict(found=0, n=10), "-0.10|20": dict(found=8, n=10), "-0.10|40": dict(found=10, n=10)}
    v, _ = snr90([tab], 0.1, snrs=[10, 20, 40])
    checks = [("noise only: no hit at true SNR 8", len(h0) == 0),
              ("all 10 plants found with correct drift", all(w["status"] == "found" for w in why)),
              ("snr90 log-interpolates (0.8 at 20, 1.0 at 40 -> ~28)", v is not None and 26 < v < 30)]
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    if sys.argv[1:2] == ["selftest"]:
        sys.exit(selftest())
    if sys.argv[1:2] == ["summary"]:
        summarise(sys.argv[2:])
        sys.exit(0)
    main()
