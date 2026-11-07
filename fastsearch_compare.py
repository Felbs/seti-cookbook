#!/usr/bin/env python3
"""fastsearch_compare.py - which turboSETI hits does fastsearch also find, and why not the rest?

A turboSETI hit counts as found if a fastsearch hit's track (start frequency,
drift) stays within TOL channels (+ the tone's per-spectrum smear) of the
turboSETI track at the first AND the last spectrum. Not one-to-one: this asks
"is every turboSETI signal also in fastsearch's list".

Every miss is then probed (--probe FILE.h5): the coarse channel is read again,
normalised exactly as fastsearch does, and fastsearch's detection statistic is
evaluated along the turboSETI track itself (best over +-3 channels, +-2 drift
steps). The reason is then one of:

  below-threshold  the statistic along that track is under fastsearch's threshold
                   (turboSETI's SNR there was inflated: its noise estimate is
                   lower than the true noise in that channel, or its own SNR is
                   biased by the bandpass shape);
  above-threshold-unmatched
                   the statistic is above threshold there but no fastsearch hit
                   matches the track: usually a different (often wider, stronger)
                   track was kept nearby, e.g. inside broadband interference, or
                   fastsearch found it a few channels off. The kept neighbour is
                   NOT identified by this tool; look in the CSV;
  capped           the coarse channel hit fastsearch's per-coarse-channel cap
                   (the CSV's `capped` column, or >= the default cap if absent);
  dc               the turboSETI hit sits on the blanked DC channel;
  other            none of the above (reported, never hidden).

    python fastsearch_compare.py TS.dat FS.fastsearch.csv [--probe FILE.h5] [--json out.json]
    python fastsearch_compare.py selftest
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fastsearch as F  # noqa: E402

TOL = 3.0


def read_ts_dat(path):
    hdr = {}
    rows = []
    for ln in open(path):
        if ln.startswith("#"):
            if "DELTAT" in ln:
                p = ln.replace("#", "").split("\t")
                hdr["tsamp"] = float(p[0].split(":")[1])
                hdr["foff_hz"] = float(p[1].split(":")[1])
                hdr["obs_length"] = float(p[3].split(":")[1])
            continue
        p = ln.split()
        if len(p) < 12:
            continue
        rows.append(dict(drift=float(p[1]), snr=float(p[2]), mhz=float(p[3]), index=int(p[5]),
                         coarse=int(p[10])))
    return hdr, rows


def read_fs_csv(path):
    out = []
    with open(path) as f:
        for r in csv.DictReader(f):
            out.append(dict(mhz=float(r["freq_mhz_t0"]), drift=float(r["drift_hz_s"]), snr=float(r["snr"]),
                            snr_power=float(r["snr_power"]), coarse=int(r["coarse_chan"]),
                            w=int(r["boxcar_w"]), n_above=int(r["n_above_in_coarse"]),
                            capped=int(r["capped"]) if r.get("capped") not in (None, "") else None))
    return out


def tracks_close(a_mhz, a_dr, b_mhz, b_dr, foff_hz, T, tol=TOL, smear=0.0):
    """Channels between the two tracks at t=0 and t=T (both ends) <= tol + smear."""
    d0 = abs(a_mhz - b_mhz) * 1e6 / abs(foff_hz)
    d1 = abs((a_mhz + a_dr * T / 1e6) - (b_mhz + b_dr * T / 1e6)) * 1e6 / abs(foff_hz)
    return max(d0, d1) <= tol + smear


def compare(ts_rows, fs_rows, foff_hz, tsamp, nt, drift_tol=None):
    T = (nt - 1) * tsamp
    fs_sorted = sorted(fs_rows, key=lambda r: r["mhz"])
    fm = np.array([r["mhz"] for r in fs_sorted])
    win = (4.0 * T + 30 * abs(foff_hz)) / 1e6
    res = []
    for h in ts_rows:
        lo, hi = np.searchsorted(fm, h["mhz"] - win), np.searchsorted(fm, h["mhz"] + win)
        smear = abs(h["drift"]) * tsamp / abs(foff_hz)
        best = None
        for j in range(lo, hi):
            r = fs_sorted[j]
            if tracks_close(h["mhz"], h["drift"], r["mhz"], r["drift"], foff_hz, T, TOL, smear):
                if best is None or r["snr"] > best["snr"]:
                    best = r
        res.append(best)
    return res


def probe(h5, misses, nfpc, max_drift, thr):
    """fastsearch's statistic along each missed turboSETI track (reads only those coarse channels)."""
    import h5py
    try:
        import hdf5plugin  # noqa: F401
    except ImportError:
        pass
    out = {}
    by_cc = {}
    for i, h in misses:
        by_cc.setdefault(h["coarse"], []).append((i, h))
    with h5py.File(h5, "r") as f:
        hdr = F.header_of(f)
        ds = f["data"]
        nt = ds.shape[0]
        foff_hz, ts = float(hdr["foff"]) * 1e6, float(hdr["tsamp"])
        step = abs(foff_hz) / ((nt - 1) * ts)
        K = int(np.ceil(max_drift / step - 1e-9))
        for cc, items in sorted(by_cc.items()):
            raw = ds[:, 0, cc * nfpc:(cc + 1) * nfpc]
            x, S = F.normalise(raw, (nfpc // 2,))
            # same calibration as the search: each width's median and spread from a real
            # drift row of that width (full resolution here)
            cal = F.calibrate_rows(x, 1, {w: [w * (nt - 1)] for w in range(1, 27)}, nt)
            P = K + 30
            xpad = np.zeros((nt, nfpc + 2 * P), np.float32)
            xpad[:, P:P + nfpc] = x
            for i, h in items:
                k0 = int(round(h["drift"] * ts * (nt - 1) / foff_hz))
                c0 = h["index"]
                bestv = -np.inf
                for k in range(max(-K, k0 - 2), min(K, k0 + 2) + 1):
                    w = max(1, int(round(abs(k) / (nt - 1))))
                    lo = max(0, c0 - 3 - w)
                    hi = min(nfpc, c0 + 4 + w)
                    D = np.zeros(hi - lo, np.float64)
                    for t in range(nt):
                        o = P + F._shift(k, t, nt, 1)
                        D += xpad[t, o + lo:o + hi]
                    _, off, sig = cal[min(w, 26)]
                    B = (F.boxcar(D.astype(np.float32), w, 1 if k >= 0 else -1) - off) / sig
                    cc_lo, cc_hi = c0 - 3 - lo, c0 + 4 - lo
                    bestv = max(bestv, float(np.max(B[max(0, cc_lo):cc_hi])))
                out[i] = dict(stat=bestv, dc=abs(c0 - nfpc // 2) <= F.DC_HALF + 1)
            del raw, x, xpad
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ts_dat")
    ap.add_argument("fs_csv")
    ap.add_argument("--nt", type=int, default=16)
    ap.add_argument("--probe", help="the .h5 file both searched, to explain misses")
    ap.add_argument("--thr", type=float, default=7.0, help="fastsearch threshold used")
    ap.add_argument("--max-drift", type=float, default=4.0)
    ap.add_argument("--nfpc", type=int, default=1033216)
    ap.add_argument("--json")
    a = ap.parse_args()
    hdr, ts_rows = read_ts_dat(a.ts_dat)
    fs_rows = read_fs_csv(a.fs_csv)
    foff_hz, tsamp = hdr["foff_hz"], hdr["tsamp"]
    T = (a.nt - 1) * tsamp
    m = compare(ts_rows, fs_rows, foff_hz, tsamp, a.nt)
    n = len(ts_rows)
    found = [i for i, r in enumerate(m) if r is not None]
    drift_good = [i for i in found if abs(m[i]["drift"] - ts_rows[i]["drift"]) <=
                  max(0.1, 0.1 * abs(ts_rows[i]["drift"]), abs(foff_hz) / T)]
    print(f"turboSETI {n} hits, fastsearch {len(fs_rows)} hits")
    print(f"turboSETI hits also found by fastsearch: {len(found)} ({100*len(found)/max(n,1):.1f} %), "
          f"of which drift agrees: {len(drift_good)}")
    from collections import Counter
    per_cc = Counter(r["coarse"] for r in fs_rows)
    capped = {r["coarse"] for r in fs_rows if r.get("capped")}
    if all(r.get("capped") is None for r in fs_rows):        # older CSV without the column
        capped = {c for c, v in per_cc.items() if v >= F.MAX_HITS_PER_COARSE}
    misses = [(i, h) for i, h in enumerate(ts_rows) if m[i] is None]
    reasons = {}
    detail = []
    pr = probe(a.probe, misses, a.nfpc, a.max_drift, a.thr) if (a.probe and misses) else {}
    for i, h in misses:
        p = pr.get(i)
        if h["coarse"] in capped:
            why = "capped"
        elif p and p["dc"]:
            why = "dc"
        elif p and p["stat"] < a.thr:
            why = "below-threshold"
        elif p:
            why = "above-threshold-unmatched"
        else:
            why = "not-probed"
        reasons[why] = reasons.get(why, 0) + 1
        detail.append(dict(mhz=h["mhz"], drift=h["drift"], ts_snr=h["snr"], coarse=h["coarse"],
                           reason=why, fs_stat=None if not p else round(p["stat"], 2)))
    print("misses by reason:", reasons)
    for d in sorted(detail, key=lambda d: -d["ts_snr"])[:25]:
        print(f"  {d['mhz']:.6f} MHz drift {d['drift']:+.4f} turboSETI SNR {d['ts_snr']:7.1f} "
              f"cc {d['coarse']:3d}: {d['reason']}" + (f" (fastsearch stat {d['fs_stat']})" if d["fs_stat"] is not None else ""))
    if pr:
        st = np.array([d["fs_stat"] for d in detail if d["fs_stat"] is not None])
        tss = np.array([d["ts_snr"] for d in detail if d["fs_stat"] is not None])
        if st.size:
            print(f"probed misses: fastsearch statistic median {np.median(st):.2f}, "
                  f"turboSETI SNR median {np.median(tss):.1f}")
    if a.json:
        Path(a.json).write_text(json.dumps(dict(n_ts=n, n_fs=len(fs_rows), found=len(found),
                                                drift_agrees=len(drift_good), reasons=reasons,
                                                misses=detail), indent=1))
    return 0


def selftest():
    foff, ts, nt = -2.835503418452676, 17.986224128, 16
    T = (nt - 1) * ts
    f0 = 1500.0
    ts_rows = [dict(mhz=f0, drift=0.5, snr=20, index=0, coarse=0),
               dict(mhz=f0 - 0.01, drift=-1.0, snr=15, index=0, coarse=0),
               dict(mhz=f0 - 0.02, drift=2.0, snr=12, index=0, coarse=0)]
    fs_rows = [dict(mhz=f0 + 1.5 * foff / 1e6, drift=0.5, snr=60, snr_power=60, coarse=0, w=3, n_above=1),
               dict(mhz=f0 - 0.01, drift=+1.0, snr=40, snr_power=40, coarse=0, w=6, n_above=1),
               dict(mhz=f0 - 0.02 + 25 * foff / 1e6, drift=2.0, snr=40, snr_power=40, coarse=0, w=6, n_above=1)]
    m = compare(ts_rows, fs_rows, foff, ts, nt)
    checks = [("same track 1.5 channels off is a match", m[0] is fs_rows[0]),
              ("opposite drift is not a match", m[1] is None),
              ("same drift 25 channels off (> 3 + smear 12.7) is not a match", m[2] is None)]
    checks.append(("track test uses both ends", not tracks_close(f0, 0.0, f0, 0.2, foff, T)))
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(selftest() if sys.argv[1:2] == ["selftest"] else main())
