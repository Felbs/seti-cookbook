#!/usr/bin/env python3
"""rfidb_falseveto.py - how often does rfi_db's veto kill a signal that IS unique?

Plants Doppler-smeared tones (inject.plant, setigen) into one REAL coarse channel
of a GJ699 ON scan, searches with turboSETI exactly like the GJ699 cadence search
(max_drift 4, min_drift 0, snr 10), matches hits to plants one-to-one along the
track (inject.match), and asks the RFI database about every recovered plant,
leaving GJ699 out (its own scans cannot vouch against it). A plant is unique by
construction, so every veto is a FALSE veto.

Two placements:
  in place   the plants' real frequencies (1418.6-1421.6 MHz)
  band-wide  each recovered hit (its measured frequency error and drift kept)
             moved to a random frequency in 1.1-1.9 GHz - the false-veto rate a
             unique signal would meet anywhere in the band, given the DB's
             real density there.

    python tools/rfidb_falseveto.py data/hi_..._GJ699_0003.gpuspec.0000.npz --db DB[,DB2]

Two-machine form (setigen on one box, turboSETI + the 29 M-hit DB on another):
    python tools/rfidb_falseveto.py NPZ --plant-only OUT --own-dat SCAN.dat   # writes OUT.npz/OUT.json
    python tools/rfidb_falseveto.py OUT.npz --planted OUT.json --db ...
--own-dat takes the data's own lines from the full-file turboSETI .dat of the same
scan (same settings) instead of a null-probe search.
"""
import json
import argparse
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import inject  # noqa: E402
import recovery  # noqa: E402
import rfi_db  # noqa: E402
import recipe_api as R  # noqa: E402

DRIFTS = [0.0, 0.05, -0.05, 0.12, -0.12, 0.3, -0.3, 1.0, -1.0, 2.0, -2.0, 3.5, -3.5]
SNRS = [100, 300, 1000]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--db")
    ap.add_argument("--per-cell", type=int, default=6)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--moves", type=int, default=20, help="band-wide placements per recovered plant")
    ap.add_argument("--plant-only", metavar="OUT")
    ap.add_argument("--own-dat")
    ap.add_argument("--planted", metavar="JSON")
    a = ap.parse_args()
    if a.planted:
        return score(a)
    data, m = recovery.load(a.npz)
    nt, nc = data.shape
    foff_hz, ts = float(m["foff"]) * 1e6, float(m["tsamp"])
    T = nt * ts
    radius = int(np.ceil(4.0 * (nt - 1) * ts / abs(foff_hz)))
    out = Path(tempfile.mkdtemp(prefix="falseveto_"))
    try:
        if a.own_dat:
            own = [h for h in inject.read_dat(a.own_dat, m["fch1"], m["foff"]) if 0 <= h["ch"] < nc]
        else:
            f0 = out / "null.h5"
            recovery.write_h5(f0, data, m)
            own = inject.read_dat(recovery.search(f0, 4.0, 10.0), m["fch1"], m["foff"])
        cells = [(d, s) for d in DRIFTS for s in SNRS]
        n = len(cells) * a.per_cell
        starts = recovery.layout(nc, n, [h["ch"] for h in own], guard=2 * radius + 50, seed=a.seed)
        plants = [dict(ch=float(starts[i]), drift=cells[i % len(cells)][0], snr=cells[i % len(cells)][1])
                  for i in range(n)]
        d = data.copy()
        inject.plant(d, foff_hz, ts, plants, seed=a.seed)
        if a.plant_only:
            meta = {k: (float(v) if np.ndim(v) == 0 and not isinstance(v, str) else str(v))
                    for k, v in m.items()}
            np.savez(a.plant_only + ".npz", data=d, **{k: v for k, v in m.items()})
            Path(a.plant_only + ".json").write_text(json.dumps(dict(
                plants=plants, own=[h["ch"] for h in own], meta=meta, seed=a.seed)))
            print(f"planted {n} tones -> {a.plant_only}.npz/.json ({len(own)} own lines avoided)")
            return 0
        f1 = out / "planted.h5"
        recovery.write_h5(f1, d, m)
        hits = inject.read_dat(recovery.search(f1, 4.0, 10.0), m["fch1"], m["foff"])
        hits = [h for h in hits if all(abs(h["ch"] - o["ch"]) > 3 for o in own)]
        matched = inject.match(plants, hits, foff_hz, ts, nt)
    finally:
        shutil.rmtree(out, ignore_errors=True)
    return report(a, plants, hits, matched, m, foff_hz, ts, nt, n)


def score(a):
    """Second half of the two-machine form: search the planted file, match, query."""
    js = json.loads(Path(a.planted).read_text())
    plants, own_ch = js["plants"], js["own"]
    z = np.load(a.npz, allow_pickle=True)
    d = z["data"].astype(np.float32)
    m = {k: z[k].item() if z[k].ndim == 0 else z[k] for k in z.files if k != "data"}
    nt, nc = d.shape
    foff_hz, ts = float(m["foff"]) * 1e6, float(m["tsamp"])
    out = Path(tempfile.mkdtemp(prefix="falseveto_"))
    try:
        f1 = out / "planted.h5"
        recovery.write_h5(f1, d, m)
        hits = inject.read_dat(recovery.search(f1, 4.0, 10.0), m["fch1"], m["foff"])
    finally:
        shutil.rmtree(out, ignore_errors=True)
    hits = [h for h in hits if all(abs(h["ch"] - o) > 3 for o in own_ch)]
    matched = inject.match(plants, hits, foff_hz, ts, nt)
    return report(a, plants, hits, matched, m, foff_hz, ts, nt, len(plants))


def report(a, plants, hits, matched, m, foff_hz, ts, nt, n):
    db = rfi_db.default_db(a.db)
    print(f"DB: {db.n:,} hits, {len(db.fam_vocab):,} target families")
    T = nt * ts
    rec = [(plants[i], hits[j]) for i, j in enumerate(matched) if j is not None]
    print(f"{n} smeared plants (drifts {DRIFTS}, TRUE SNR {SNRS}); recovered {len(rec)}")
    mjd = float(m["tstart"])
    res = abs(foff_hz)
    rng = np.random.default_rng(a.seed)
    inplace = []
    moved, moved_out = [], []
    by_drift = {}
    for p, h in rec:
        v = db.query(h["mhz"], h["drift"], mjd=mjd, target="GJ699", tobs_s=T, res_hz=res)
        inplace.append(v["veto"])
        for _ in range(a.moves):
            f = rng.uniform(1100.0, 1900.0)
            vv = db.query(f, h["drift"], mjd=mjd, target="GJ699", tobs_s=T, res_hz=res)["veto"]
            moved.append(vv)
            if R.rfi_reason(f) is None:
                moved_out.append(vv)
            by_drift.setdefault(p["drift"], []).append(vv)
    print(f"false vetoes IN PLACE (1418.6-1421.6 MHz): {sum(inplace)}/{len(inplace)} "
          f"= {100 * np.mean(inplace):.2f} %")
    print(f"false vetoes BAND-WIDE 1.1-1.9 GHz: {sum(moved)}/{len(moved)} = {100 * np.mean(moved):.2f} %; "
          f"outside RFI_BANDS: {sum(moved_out)}/{len(moved_out)} = {100 * np.mean(moved_out):.2f} %")
    for dr in sorted(by_drift):
        x = by_drift[dr]
        print(f"   planted drift {dr:+5.2f} Hz/s: {100 * np.mean(x):5.2f} % of {len(x)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
