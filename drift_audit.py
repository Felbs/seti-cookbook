#!/usr/bin/env python3
"""drift_audit.py - can this file's drift search tell a planet from RFI at all?

Plants tones with KNOWN drift rates into a slice of a REAL Breakthrough Listen
file, runs turboSETI on it, and prints what drift each planted tone came back
with. If every planted tone - from 0 to 15 Hz/s, either sign - comes back at
the same one or two values, then on this kind of file the reported drift rate
carries no information, and ANY veto built on it ("all hits share one drift
rate, so they are instrumental") is void.

That is exactly what happens on BL mid-resolution (.0002, 2.86 kHz x 1.07 s)
products: the drift quantum is ~5-10 Hz/s, wider than the whole range a planet
can produce (~0.05-4 Hz/s at L/S band). See SWEEP_RESULTS.md for the history.

    python drift_audit.py data/star_GJ144.h5            # Eps Eridani, mid-res
    python drift_audit.py data/star_GJ144.h5 --keep     # keep the injected slice

Needs turbo_seti (+ tools/patch_turboseti.py on numpy 2) and blimpy.
"""
import argparse
import os
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

PLANTED = [0.0, 0.05, 0.3, 1.0, 3.0, 5.2, 8.0, 15.0, -0.3, -3.0]   # Hz/s
SNR_TARGET = 60          # integrated strength, in turboSETI's own SNR units


def drift_quantum(h5):
    """The drift resolution a de-Doppler search can have on this file, Hz/s.
    One channel per observation length; turboSETI pads the time axis to a power
    of two, which is why it reports multiples of about half of this."""
    import h5py, hdf5plugin  # noqa: F401
    with h5py.File(h5, "r") as f:
        a = dict(f["data"].attrs)
        nt = f["data"].shape[0]
    return abs(a["foff"]) * 1e6 / ((nt - 1) * a["tsamp"]), a, nt


def plant(h5, out):
    import h5py, hdf5plugin  # noqa: F401
    from blimpy import Waterfall
    with h5py.File(h5, "r") as f:
        a = dict(f["data"].attrs)
        nchan = f["data"].shape[2]
    fpc = 1024 if nchan % 1024 == 0 else nchan // max(1, round(nchan / 1024))
    c0 = (nchan // 2 // fpc) * fpc                       # 4 coarse channels mid-band
    c1 = c0 + 4 * fpc
    fa, fb = a["fch1"] + a["foff"] * c0, a["fch1"] + a["foff"] * c1
    wf = Waterfall(h5, f_start=min(fa, fb), f_stop=max(fa, fb))
    d = wf.data.astype(np.float64)
    nt, _, nc = d.shape
    ts, foff_hz = wf.header["tsamp"], wf.header["foff"] * 1e6
    df = abs(foff_hz)
    t = np.arange(nt) * ts
    spacing = nc // (len(PLANTED) + 1)
    truth = []
    for i, dr in enumerate(PLANTED):
        ch0 = spacing // 2 + i * spacing
        cc = ch0 // fpc
        amp = SNR_TARGET * np.std(d[:, 0, cc * fpc:(cc + 1) * fpc].sum(0)) / nt
        pos = ch0 + dr * t / foff_hz          # dr is in Hz/s; channel axis may run downward
        for k in range(nt):
            j = int(np.floor(pos[k]))
            fr = pos[k] - j
            if 0 <= j < nc - 1:
                d[k, 0, j] += amp * (1 - fr)
                d[k, 0, j + 1] += amp * fr
        truth.append((dr, ch0))
    wf.data = d.astype(np.float32)
    fn = Path(out) / "planted.h5"
    wf.write_to_hdf5(str(fn))
    return fn, truth, wf.get_freqs(), df, nt * ts


def search(fn, **kw):
    from turbo_seti.find_doppler.find_doppler import FindDoppler
    cwd = os.getcwd()
    os.chdir(fn.parent)
    try:
        FindDoppler(fn.name, snr=10, out_dir=".", append_output=False,
                    log_level_int=40, **kw).search()
    finally:
        os.chdir(cwd)
    hits = []
    for ln in open(fn.with_suffix(".dat")):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        hits.append((float(p[1]), float(p[2]), float(p[3])))
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("h5")
    ap.add_argument("--keep", action="store_true", help="keep the planted slice")
    args = ap.parse_args()

    q, a, nt = drift_quantum(args.h5)
    print(f"file: {args.h5}")
    print(f"  {abs(a['foff'])*1e6:.2f} Hz channels x {a['tsamp']:.3f} s x {nt} spectra"
          f"  ->  drift quantum ~{q:.2f} Hz/s")
    out = tempfile.mkdtemp(prefix="drift_audit_")
    try:
        fn, truth, freqs, df, T = plant(args.h5, out)
        settings = [("as the old star_sweep ran it (max_drift=8)", dict(max_drift=8.0)),
                    ("wide open (max_drift=20, min_drift=0)",
                     dict(max_drift=20.0, min_drift=0.0))]
        seen_all, score = set(), []
        for label, kw in settings:
            hits = search(fn, **kw)
            print(f"\n{label}:")
            reported, rep_by_tone = [], []
            for dr, ch0 in truth:
                f0 = freqs[ch0]
                # turboSETI reports the track's frequency at t=0, which is where we planted it
                m = sorted((h for h in hits if abs(h[2] - f0) * 1e6 <= 2.5 * df),
                           key=lambda h: -h[1])
                rep_by_tone.append(m[0][0] if m else None)
                if m:
                    reported.append(m[0][0])
                    print(f"  planted {dr:+6.2f} Hz/s  ->  reported {m[0][0]:+8.3f} Hz/s  (SNR {m[0][1]:.0f})")
                else:
                    print(f"  planted {dr:+6.2f} Hz/s  ->  MISSED")
            distinct = sorted({round(r, 3) for r in reported})
            seen_all |= set(distinct)
            # "correct" = good enough to tell planets apart: within 0.1 Hz/s or 10 %.
            # (A tolerance in drift QUANTA would be circular: on mid-res files one
            # quantum is wider than the whole planetary range.)
            inrange = [(dr, r) for (dr, _), r in zip(truth, rep_by_tone) if abs(dr) <= kw["max_drift"]]
            ok = sum(1 for dr, r in inrange
                     if r is not None and abs(r - dr) <= max(0.1, 0.10 * abs(dr)))
            score.append((label, ok, len(inrange)))
            print(f"  {len(reported)}/{len(truth)} found; drift CORRECT for {ok}/{len(inrange)}"
                  f" in range; distinct values reported: {distinct}")
        print()
        lab, ok, n = score[0]
        print(f"old-sweep settings: drift correct for {ok}/{n} planted tones")
        if ok < 0.7 * n:
            print("VERDICT: drift is UNRESOLVED on this file. Ten different planted drifts came back")
            print("as the same few values, so 'all hits share one drift rate' says nothing about RFI.")
            print("A drift-based verdict needs fine-resolution (.0000, ~2.8 Hz) data.")
        else:
            print("VERDICT: drift is resolved on this file; drift-based vetoes are meaningful.")
            print("(Keep max_drift modest: on fine-res files a very wide max_drift aliases reports.)")
    finally:
        if args.keep:
            print(f"(planted slice kept in {out})")
        else:
            shutil.rmtree(out, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
