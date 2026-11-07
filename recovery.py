#!/usr/bin/env python3
"""recovery.py - how often does our search find a planted signal in REAL data?

Plants a grid of setigen tones (drift x strength) into one coarse channel of a
real Breakthrough Listen scan, runs turboSETI with exactly the settings of the
search being calibrated, and reports, for every grid cell, the fraction found
with the correct drift and why the rest were lost (inject.py). The same file
is searched once WITHOUT plants first (the null probe): those hits are the
data's own lines, and plants are kept away from them.

    python recovery.py data/hi_..._GJ699_0003.gpuspec.0000.npz
    python recovery.py FILE.npz --max-drift 4 --snr 10 --seed 1 --json out.json

The default settings are the GJ699 full-cadence search (R9, GJ699_CADENCE_RESULT.md):
max_drift 4, min_drift 0, snr 10.

Only the part of the slice inside ONE true coarse channel is searched
(coarse_align.py; 1,033,216 fine channels per coarse channel in the 2016 GJ699
files), with its DC spike repaired. The hi_starframe slices of 2^20 channels
straddle two coarse channels (KNOWN_ISSUES #8); --whole-slice restores the old
behaviour (turboSETI on the whole straddling slice), which the first published
numbers used.
"""
import argparse
import json
import os
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import inject  # noqa: E402

DRIFTS = [0.0, 0.1, -0.1, 0.3, -0.3, 1.0, -1.0, 2.0, -2.0, 3.5, -3.5]   # Hz/s
SNRS = [10, 20, 40, 60, 100, 200]   # TRUE radiometer SNR


def write_h5(path, data, m):
    """Minimal blimpy-readable filterbank HDF5 for one coarse channel."""
    import h5py
    with h5py.File(path, "w") as f:
        f.attrs["CLASS"] = "FILTERBANK"
        f.attrs["VERSION"] = "1.0"
        ds = f.create_dataset("data", data=data[:, None, :].astype(np.float32))
        hdr = dict(fch1=float(m["fch1"]), foff=float(m["foff"]), nchans=data.shape[1], nifs=1,
                   nbits=32, tsamp=float(m["tsamp"]), tstart=float(m["tstart"]),
                   src_raj=float(m["ra_h"]), src_dej=float(m["dec_deg"]), source_name="slice",
                   telescope_id=6, machine_id=10, data_type=1, az_start=0.0, za_start=0.0,
                   rawdatafile="slice")
        for k, v in hdr.items():
            ds.attrs[k] = v
        ds.dims[0].label = b"frequency"


def load(path):
    z = np.load(path, allow_pickle=True)
    m = {k: z[k].item() if z[k].ndim == 0 else z[k] for k in z.files if k != "data"}
    data = z["data"].astype(np.float32)
    # npz from hi_starframe: fch1/foff describe the WHOLE file; ch0 = first channel of the slice
    m["fch1"] = float(m["fch1"]) + float(m["foff"]) * int(m.get("ch0", 0))
    return data, m


def search(fn, max_drift, snr):
    from turbo_seti.find_doppler.find_doppler import FindDoppler
    cwd = os.getcwd()
    os.chdir(fn.parent)
    try:
        FindDoppler(fn.name, max_drift=max_drift, min_drift=0.0, snr=snr, out_dir=".",
                    n_coarse_chan=1, append_output=False, log_level_int=40).search()
    finally:
        os.chdir(cwd)
    return fn.with_suffix(".dat")


def layout(nc, n, avoid, guard, seed):
    """n start channels, >= guard apart, away from the DC spike and real lines."""
    rng = np.random.default_rng(seed)
    centre = nc // 2
    slots = np.arange(guard, nc - guard, guard)
    ok = [s for s in slots if abs(s - centre) > 2 * guard and all(abs(s - a) > guard for a in avoid)]
    ok = np.array(ok, dtype=float)
    rng.shuffle(ok)
    if len(ok) < n:
        raise SystemExit(f"only {len(ok)} clean slots for {n} plants")
    return ok[:n] + rng.uniform(0, 1, n)          # sub-channel start positions


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("npz")
    ap.add_argument("--max-drift", type=float, default=4.0)
    ap.add_argument("--snr", type=float, default=10.0)
    ap.add_argument("--per-cell", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--drifts", help="comma list of Hz/s (default: the standard grid)")
    ap.add_argument("--snrs", help="comma list of TRUE SNR (default: the standard grid)")
    ap.add_argument("--json")
    ap.add_argument("--whole-slice", action="store_true",
                    help="search the whole npz slice even if it straddles two coarse channels (old behaviour)")
    a = ap.parse_args()
    global DRIFTS, SNRS
    if a.drifts:
        DRIFTS = [float(v) for v in a.drifts.split(",")]
    if a.snrs:
        SNRS = [float(v) for v in a.snrs.split(",")]

    data, m = load(a.npz)
    seg_info = None
    if not a.whole_slice:
        import coarse_align
        nf_slice = data.shape[1]
        data, seg_info = coarse_align.coarse_segment(data, int(m.get("ch0", 0)), float(m["foff"]))
        m["fch1"] = float(m["fch1"]) + float(m["foff"]) * seg_info["lo"]
        print(coarse_align.describe(seg_info, nf_slice))
    nt, nc = data.shape
    foff_hz, ts = float(m["foff"]) * 1e6, float(m["tsamp"])
    T = (nt - 1) * ts
    radius = int(np.ceil(a.max_drift * T / abs(foff_hz)))
    out = Path(tempfile.mkdtemp(prefix="recovery_"))
    try:
        # 1. null probe: the file's own hits, unplanted
        f0 = out / "null.h5"
        write_h5(f0, data, m)
        own = inject.read_dat(search(f0, a.max_drift, a.snr), m["fch1"], m["foff"])
        print(f"{Path(a.npz).name}: {nt} x {nc} ch, {abs(foff_hz):.2f} Hz x {ts:.1f} s; "
              f"turboSETI max_drift {a.max_drift}, snr {a.snr}")
        print(f"null probe (no plants): {len(own)} hits of the data's own")
        k = inject.turboseti_noise(data) / inject.radiometer_noise(data)
        print(f"turboSETI noise / radiometer noise = {k:.2f}: turboSETI SNR {a.snr:g} here "
              f"means TRUE SNR ~{a.snr * k:.0f} for a channel-centred zero-drift tone")

        # 2. plant the grid
        cells = [(d, s) for d in DRIFTS for s in SNRS]
        n = len(cells) * a.per_cell
        avoid = [h["ch"] for h in own]
        if seg_info and seg_info["dc"] is not None:
            avoid.append(seg_info["dc"])          # the real DC spike (repaired, but keep plants off it)
        starts = layout(nc, n, avoid, guard=2 * radius + 50, seed=a.seed)
        plants = [dict(ch=float(starts[i]), drift=cells[i % len(cells)][0], snr=cells[i % len(cells)][1])
                  for i in range(n)]
        d = data.copy()
        inject.plant(d, foff_hz, ts, plants, seed=a.seed)
        f1 = out / "planted.h5"
        write_h5(f1, d, m)
        hits = inject.read_dat(search(f1, a.max_drift, a.snr), m["fch1"], m["foff"])
        hits = [h for h in hits if all(abs(h["ch"] - o["ch"]) > 3 for o in own)]
        matched = inject.match(plants, hits, foff_hz, ts, nt)
        why = inject.explain_misses(plants, hits, matched, foff_hz, ts, nt, a.max_drift)

        # 3. table
        print(f"planted {n} tones ({a.per_cell} per cell, seed {a.seed}); smeared (setigen)")
        print("\nfraction found WITH correct drift   (rows: drift Hz/s, columns: planted TRUE SNR)")
        print("drift  " + "".join(f"{s:>7g}" for s in SNRS) + "   misses: wrong-drift/suppressed/below")
        table = {}
        for dr in DRIFTS:
            row, reasons = [], {"found-wrong-drift": 0, "wrong-path": 0, "suppressed": 0, "below-threshold": 0}
            for s in SNRS:
                idx = [i for i, p in enumerate(plants) if p["drift"] == dr and p["snr"] == s]
                good = sum(1 for i in idx if why[i]["status"] == "found")
                for i in idx:
                    if why[i]["status"] in reasons:
                        reasons[why[i]["status"]] += 1
                row.append(good / len(idx))
                table[f"{dr:+.2f}|{s:g}"] = dict(found=good, n=len(idx),
                                               statuses=[why[i]["status"] for i in idx])
            print(f"{dr:+5.2f}  " + "".join(f"{r:7.2f}" for r in row) +
                  f"   {reasons['found-wrong-drift'] + reasons['wrong-path']}/"
                  f"{reasons['suppressed']}/{reasons['below-threshold']}")
        # zero-drift SNR calibration: does turboSETI report what we planted?
        z = [(p["snr"], hits[matched[i]]["snr"]) for i, p in enumerate(plants)
             if p["drift"] == 0 and matched[i] is not None]
        if z:
            ratio = np.median([r / p for p, r in z])
            print(f"\nzero-drift check: turboSETI-reported / planted true SNR median {ratio:.3f} "
                  f"(expect ~1/{k:.1f} = {1/k:.3f}, less when a tone straddles two channels)")
        if a.json:
            Path(a.json).write_text(json.dumps(dict(file=Path(a.npz).name, max_drift=a.max_drift,
                                                    snr=a.snr, seed=a.seed, per_cell=a.per_cell,
                                                    null_hits=len(own), noise_ratio=k, table=table,
                                                    nchan=nc, coarse_segment=seg_info), indent=1))
    finally:
        shutil.rmtree(out, ignore_errors=True)


if __name__ == "__main__":
    main()
