#!/usr/bin/env python3
"""fastsearch_cadence.py - one command for a whole ON/OFF cadence with fastsearch.

Given a directory holding the six .h5 files of one cadence (ON OFF ON OFF ON OFF),
it runs, in order, each stage resumable (a stage whose outputs exist is skipped):

  1. fastsearch.py on every file (workers, MemAvailable guard; a file whose .dat
     and .summary.json exist with the same settings is skipped);
  2. turboSETI's find_event_pipeline on the fastsearch .dat files, filter 1, 2, 3
     (no turboSETI installed: a small built-in ON/OFF filter of the same shape, labelled
     "builtin" everywhere - not turboSETI, used by the offline selftest); then every
     events CSV gets a `fastsearch_broadband` column (fastsearch.flag_events);
  3. cadence_forced.py run on the chosen filters' events (measured in every scan);
  4. cadence_verdict.py --forced --rfi-db on the same events; its report is saved.

Scan order and roles follow bl_archive.py's rules: files sorted by the MJD and
seconds in the guppi name (header tstart if absent); OFF = the target name + "_OFF";
one data product; the six must alternate ON OFF ON OFF ON OFF, else it stops.

    python fastsearch_cadence.py DIR --out-dir OUT --dry-run      # plan, time, disk, RAM
    python fastsearch_cadence.py DIR --out-dir OUT --workers 3 --snr 8
    python fastsearch_cadence.py selftest                           # six tiny synthetic h5 files

Outputs in OUT: fs/<stem>.dat|.fastsearch.csv|.summary.json, dat.lst, h5.lst,
events_f{1,2,3}.csv, forced_f{N}.csv (+ .json), verdict_f{N}.txt, plan.json.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fastsearch as F  # noqa: E402

# measured on the PC (FASTSEARCH_RESULT.md (b)): one 12.8 GB file, 3 workers, 11.1 min
FS_WORKER_MIN_PER_GB = 11.1 * 3 / 12.8
FORCED_S_PER_ROW = 1056 / 1029          # FORCED_CADENCE_RESULT.md, PC
NAME_RE = re.compile(r"guppi_(\d{5})_(\d{5})_(.+?)_(\d{4})\.gpuspec\.(?:\d+\.)?(\d{4})\.h5$")


# ------------------------------------------------------------------ cadence
def scan_info(p):
    m = NAME_RE.search(p.name)
    if m:
        mjd, sec, tgt, scan, prod = m.groups()
        key = (int(mjd), int(sec))
    else:
        import h5py
        with h5py.File(p, "r") as f:
            key = (float(f["data"].attrs["tstart"]), 0)
        tgt, scan, prod = p.stem, "?", "?"
    kind = "OFF" if tgt.upper().endswith("_OFF") else "ON"
    base = tgt[:-4] if kind == "OFF" else tgt
    return dict(path=p, key=key, kind=kind, target=base, scan=scan, product=prod, size=p.stat().st_size)


def order_cadence(d):
    files = sorted(Path(d).glob("*.h5"))
    if not files:
        raise SystemExit(f"no .h5 files in {d}")
    scans = sorted((scan_info(p) for p in files), key=lambda s: s["key"])
    kinds = [s["kind"] for s in scans]
    problems = []
    if len(scans) != 6:
        problems.append(f"{len(scans)} files, a cadence has 6")
    if kinds != ["ON", "OFF"] * (len(scans) // 2) or len(scans) % 2:
        problems.append(f"order is {' '.join(kinds)}, not ON OFF ON OFF ON OFF")
    if len({s["target"] for s in scans}) != 1:
        problems.append(f"more than one target: {sorted({s['target'] for s in scans})}")
    if len({s["product"] for s in scans}) != 1:
        problems.append(f"more than one product: {sorted({s['product'] for s in scans})}")
    return scans, problems


def stem_of(p):
    return p.name[:-3] if p.name.endswith(".h5") else p.stem


# ------------------------------------------------------------------ plan
def plan_text(scans, a):
    gb = sum(s["size"] for s in scans) / 1e9
    w = max(1, a.workers)
    fs_min = FS_WORKER_MIN_PER_GB * gb / w
    nfpc_mb = None
    try:
        import h5py
        import hdf5plugin  # noqa: F401
        with h5py.File(scans[0]["path"], "r") as f:
            h = F.header_of(f)
            nt = f["data"].shape[0]
            nfpc, ncoarse = F.coarse_width(h, a.n_coarse)
            cache = F.chunk_cache_bytes(f["data"], nfpc)
        per_worker = nt * nfpc * 4 * 6 + cache
        nfpc_mb = per_worker / 2 ** 20
    except Exception as e:      # an unreadable header is reported, not fatal for the plan
        nt = nfpc = ncoarse = None
        per_worker = None
        print(f"  (could not read the first header: {e})")
    lines = [f"cadence: {scans[0]['target']}, product {scans[0]['product']}, {len(scans)} files, {gb:.1f} GB"]
    for s in scans:
        lines.append(f"  {s['kind']:3s} {s['path'].name}  {s['size'] / 1e9:.2f} GB")
    if ncoarse:
        lines.append(f"geometry: {nt} spectra, {ncoarse} coarse channels of {nfpc} fine channels")
    lines.append(f"1. fastsearch, threshold {a.snr}, max drift {a.max_drift} Hz/s, {w} worker(s): "
                 f"~{fs_min:.0f} min ({FS_WORKER_MIN_PER_GB:.2f} worker-min per GB, measured on one 12.8 GB "
                 f"file on a 4+ core PC; slower on a busy or smaller machine)")
    if per_worker:
        lines.append(f"   RAM: ~{nfpc_mb:.0f} MB per worker (slab x 6 + chunk cache) -> ~{w * nfpc_mb / 1024:.1f} GB "
                     f"for {w}; workers are reduced if MemAvailable is short")
    lines.append(f"   disk: ~20-70 MB of .dat/.csv per file (~0.1-0.4 GB total); inputs are only read")
    lines.append(f"2. find_event_pipeline filters {a.filters} ({a.event_finder}): minutes; reads headers only")
    lines.append(f"3. cadence_forced on filters {a.forced_filters}: ~{FORCED_S_PER_ROW:.1f} s per events row "
                 f"(PC); the row count is known after step 2 (turboSETI gave 1,029 filter-2 rows here; "
                 f"fastsearch, more sensitive, will give more - use --forced-limit to cap)")
    lines.append(f"4. cadence_verdict --forced --rfi-db {a.rfi_db or '(seed)'}: seconds")
    return "\n".join(lines)


# ------------------------------------------------------------------ builtin ON/OFF filter
def read_dat_df(path):
    import pandas as pd
    lines = open(path).read().splitlines()
    hdr = dict(FileID=lines[1].split(":")[-1].strip(), Source=lines[3].split(":")[-1].strip().upper())
    f4 = lines[4].split("\t")
    hdr["MJD"], hdr["RA"], hdr["DEC"] = (x.split(":", 1)[-1].strip() for x in f4[:3])
    f5 = lines[5].split("\t")
    hdr["DELTAT"], hdr["DELTAF"], hdr["max_drift_rate"], hdr["obs_length"] = (
        float(x.split(":")[-1].strip()) for x in f5[:4])
    rows = []
    for ln in lines[9:]:
        p = ln.split()
        if len(p) >= 12:
            rows.append(dict(TopHitNum=int(p[0]), DriftRate=float(p[1]), SNR=float(p[2]), Freq=float(p[3]),
                             ChanIndx=int(p[5]), FreqStart=float(p[6]), FreqEnd=float(p[7]),
                             CoarseChanNum=int(p[10]), FullNumHitsInRange=int(p[11])))
    d = pd.DataFrame(rows)
    for k, v in hdr.items():
        d[k] = v
    return d


def builtin_find_events(dats, kinds, filter_threshold, snr_cut=None, tol_ch=10.0):
    """Same SHAPE as turboSETI's find_event filters (not its code): an ON hit is
    'in' another scan if a hit there lies within tol = max(tol_ch channels, |drift| x
    obs_length) of the hit's frequency extrapolated along its drift to that scan's start.
    filter 1: ON hits >= snr_cut; 2: ... in no OFF; 3: ... and in every ON (one row per ON)."""
    import pandas as pd
    tabs = [read_dat_df(p) for p in dats]
    t0 = [float(t["MJD"].iloc[0]) * 86400 if len(t) else None for t in tabs]
    out = []
    used = set()
    hid = 0
    for i, (t, k) in enumerate(zip(tabs, kinds)):
        if k != "ON" or not len(t):
            continue
        for _, r in t.iterrows():
            if snr_cut is not None and r["SNR"] < snr_cut:
                continue
            if (i, int(r["TopHitNum"])) in used:
                continue
            tol = max(tol_ch * abs(float(r["DELTAF"])) / 1e6, abs(r["DriftRate"]) * float(r["obs_length"]) / 1e6)
            inn = {}
            for j, (u, kj) in enumerate(zip(tabs, kinds)):
                if j == i or not len(u) or t0[j] is None:
                    continue
                pred = r["Freq"] + r["DriftRate"] * (t0[j] - t0[i]) / 1e6
                near = u[(u["Freq"] - pred).abs() <= tol]
                if snr_cut is not None:
                    near = near[near["SNR"] >= snr_cut]
                if len(near):
                    inn[j] = near.iloc[(near["SNR"]).values.argmax()]
            in_off = any(kinds[j] == "OFF" for j in inn)
            ons = [j for j in inn if kinds[j] == "ON"]
            if filter_threshold >= 2 and in_off:
                continue
            if filter_threshold >= 3:
                if len(ons) != kinds.count("ON") - 1:
                    continue
                hid += 1
                for j in sorted([i] + ons):
                    rr = r if j == i else inn[j]
                    used.add((j, int(rr["TopHitNum"])))
                    out.append(dict(rr, Hit_ID=f"{r['Source']}_{hid}", status=f"on_table_{j // 2 + 1}",
                                    in_n_ons=len(ons) + 1, RFI_in_range=0, event_finder="builtin"))
            else:
                out.append(dict(r, Hit_ID="", status=f"on_table_{i // 2 + 1}", in_n_ons=len(ons) + 1,
                                RFI_in_range=int(in_off), event_finder="builtin"))
    return pd.DataFrame(out)


def find_events(out, dat_lst, h5_lst, dats, kinds, flt, finder, snr_cut):
    csv = out / f"events_f{flt}.csv"
    if finder == "turboseti":
        from turbo_seti.find_event.find_event_pipeline import find_event_pipeline
        df = find_event_pipeline(str(dat_lst), str(h5_lst), filter_threshold=flt, number_in_cadence=len(dats),
                                 on_off_first="ON", saving=True, csv_name=str(csv), SNR_cut=snr_cut,
                                 user_validation=False)
        if df is None or not csv.exists():
            import pandas as pd
            pd.DataFrame(columns=["FileID", "CoarseChanNum", "Freq", "DriftRate", "SNR", "Hit_ID"]).to_csv(
                csv, index=False)
    else:
        df = builtin_find_events(dats, kinds, flt, snr_cut)
        df.to_csv(csv, index=False)
    return csv


# ------------------------------------------------------------------ run
def run(a):
    scans, problems = order_cadence(a.dir)
    out = Path(a.out_dir)
    print(plan_text(scans, a), flush=True)
    if problems:
        print("NOT A VALID CADENCE: " + "; ".join(problems), file=sys.stderr)
        if not a.allow_bad_cadence:
            return 2
    if a.dry_run:
        print("dry run: nothing executed")
        return 0
    out.mkdir(parents=True, exist_ok=True)
    fs_dir = out / "fs"
    fs_dir.mkdir(exist_ok=True)
    t_all = time.time()
    (out / "plan.json").write_text(json.dumps(dict(
        files=[str(s["path"].name) for s in scans], kinds=[s["kind"] for s in scans], snr=a.snr,
        max_drift=a.max_drift, workers=a.workers, filters=a.filters, forced_filters=a.forced_filters,
        event_finder=a.event_finder, broadband_density=a.broadband_density, broadband_keep=a.broadband_keep),
        indent=1))
    # 1. fastsearch
    dats = []
    for s in scans:
        stem = stem_of(s["path"])
        dat, summ = fs_dir / f"{stem}.dat", fs_dir / f"{stem}.summary.json"
        if dat.exists() and summ.exists():
            j = json.loads(summ.read_text())
            if j.get("threshold_snr") == a.snr and j.get("n_coarse_searched"):
                print(f"[1] have {dat.name} ({j['n_hits']} hits) - skipped", flush=True)
                dats.append(dat)
                continue
            print(f"[1] {dat.name} exists with other settings (threshold {j.get('threshold_snr')}): re-running")
        need = 1 << 30
        av = F.mem_available_bytes()
        if av is not None and av < need + (1 << 30):
            print(f"[1] waiting: MemAvailable {av / 2**30:.2f} GB", flush=True)
            F.wait_for_memory(need, timeout=3600)
        t0 = time.time()
        F.run_file(s["path"], fs_dir, a.max_drift, a.snr, a.workers, F.parse_coarse(a.coarse), a.n_coarse, "auto",
                   verbose=a.verbose, max_hits=a.max_hits, broadband_density=a.broadband_density,
                   broadband_keep=a.broadband_keep)
        print(f"[1] {s['path'].name}: {time.time() - t0:.0f} s", flush=True)
        dats.append(dat)
    kinds = [s["kind"] for s in scans]
    dat_lst, h5_lst = out / "dat.lst", out / "h5.lst"
    dat_lst.write_text("\n".join(str(p.resolve()) for p in dats) + "\n")
    h5_lst.write_text("\n".join(str(s["path"].resolve()) for s in scans) + "\n")
    # 2. events
    finder = a.event_finder
    if finder == "auto":
        try:
            import turbo_seti  # noqa: F401
            finder = "turboseti"
        except ImportError:
            finder = "builtin"
            print("[2] turbo_seti not installed: using the BUILT-IN ON/OFF filter (not turboSETI)")
    if finder == "turboseti":
        # find_event_pipeline matches each .dat to its .h5 by basename: make the names line up
        for s, d in zip(scans, dats):
            assert d.name == s["path"].name.replace(".h5", ".dat")
    for flt in a.filters:
        csv = out / f"events_f{flt}.csv"
        if csv.exists() and not a.redo_events:
            print(f"[2] have {csv.name} - skipped")
        else:
            t0 = time.time()
            find_events(out, dat_lst, h5_lst, dats, kinds, flt, finder, a.event_snr)
            print(f"[2] filter {flt} ({finder}): {time.time() - t0:.0f} s", flush=True)
        F.flag_events(csv, fs_dir)
    # 3-4. forced + verdict
    py = sys.executable
    for flt in a.forced_filters:
        csv = out / f"events_f{flt}.csv"
        import pandas as pd
        n = len(pd.read_csv(csv)) if csv.stat().st_size > 0 else 0
        if n == 0:
            print(f"[3] filter {flt}: no events - forced check and verdict skipped")
            continue
        forced = out / f"forced_f{flt}.csv"
        if forced.exists() and not a.redo_events:
            print(f"[3] have {forced.name} - skipped")
        else:
            lim = min(n, a.forced_limit) if a.forced_limit else n
            print(f"[3] cadence_forced on {lim} rows of filter {flt} (~{lim * FORCED_S_PER_ROW / 60:.0f} min at "
                  f"the PC rate)", flush=True)
            cmd = [py, str(HERE / "cadence_forced.py"), "run", "--events", str(csv), "--files",
                   *[str(s["path"]) for s in scans], "--out", str(forced), "--json", str(forced.with_suffix(".json")),
                   "--nneg", str(a.nneg), "--dats", *[str(d) for d in dats]]
            if a.forced_limit:
                cmd += ["--limit", str(a.forced_limit)]
            r = subprocess.run(cmd)
            if r.returncode not in (0, 2):
                print(f"[3] cadence_forced failed ({r.returncode})", file=sys.stderr)
                return r.returncode
        cmd = [py, str(HERE / "cadence_verdict.py"), str(csv), "--forced", str(forced), "--rfi-db"]
        if a.rfi_db:
            cmd.append(a.rfi_db)
        r = subprocess.run(cmd, capture_output=True, text=True)
        (out / f"verdict_f{flt}.txt").write_text(r.stdout + (("\n" + r.stderr) if r.returncode else ""))
        tail = [ln for ln in r.stdout.splitlines() if "unexplained" in ln][-1:] or ["(no summary line)"]
        print(f"[4] verdict filter {flt}: {tail[0]} -> {out / f'verdict_f{flt}.txt'}", flush=True)
        if r.returncode:
            print(r.stderr[-2000:], file=sys.stderr)
            return r.returncode
    print(f"done in {(time.time() - t_all) / 60:.1f} min -> {out}")
    return 0


# ------------------------------------------------------------------ selftest
def write_tiny_cadence(d, rng, nt=16, nfc=32768, ncoarse=2):
    """Six tiny GBT-like h5 files (recovery.write_h5 layout + chunking): a sky-like
    tone in the three ONs only (phase-continuous, -0.12 Hz/s), an RFI tone in all six,
    a tone in the OFFs only. Returns the planted frequencies (MHz at the first scan)."""
    import h5py
    import inject
    import recovery
    foff, ts = -2.835503418452676, 17.986224128
    fch1 = 1650.0                 # outside every listed RFI band and natural line
    starts = [48006, 48334, 48661, 48988, 49314, 49640]
    names = []
    f_sky = fch1 + foff / 1e6 * 20000.3
    f_rfi = fch1 + foff / 1e6 * 45000.6
    f_off = fch1 + foff / 1e6 * 52000.2
    t0 = 57457 * 86400 + starts[0]
    for i, s in enumerate(starts):
        kind = "ON" if i % 2 == 0 else "OFF"
        tgt = "SYNTH" if kind == "ON" else "SYNTH_OFF"
        x = np.concatenate([F._synthetic(rng, nt, nfc) for _ in range(ncoarse)], axis=1)
        dt = (57457 * 86400 + s) - t0
        pl = [dict(ch=(f_rfi - fch1) * 1e6 / foff, drift=0.0, snr=200)]
        if kind == "ON":
            pl.append(dict(ch=(f_sky + (-0.12) * dt / 1e6 - fch1) * 1e6 / foff, drift=-0.12, snr=120))
        else:
            pl.append(dict(ch=(f_off - fch1) * 1e6 / foff, drift=0.3, snr=120))
        inject.plant(x, foff, ts, pl)
        fn = Path(d) / f"spliced_blcXX_2bit_guppi_57457_{s}_{tgt}_{i + 3:04d}.gpuspec.0000.h5"
        m = dict(fch1=fch1, foff=foff / 1e6, tsamp=ts, tstart=57457 + s / 86400, ra_h=17.963, dec_deg=4.739)
        recovery.write_h5(fn, x, m)
        with h5py.File(fn, "a") as f:          # rewrite with coarse-channel-sized chunks and the real name
            data = f["data"][...]
            attrs = dict(f["data"].attrs)
            del f["data"]
            ds = f.create_dataset("data", data=data, chunks=(1, 1, nfc))
            for k, v in attrs.items():
                ds.attrs[k] = v
            ds.attrs["source_name"] = tgt
        names.append(fn)
    return names, dict(sky=f_sky, rfi=f_rfi, off=f_off)


def selftest():
    import tempfile
    import shutil
    import pandas as pd
    checks = []
    tmp = Path(tempfile.mkdtemp(prefix="fscad_selftest_"))
    try:
        cad = tmp / "cad"
        cad.mkdir()
        names, f = write_tiny_cadence(cad, np.random.default_rng(3))
        scans, problems = order_cadence(cad)
        checks.append(("cadence ordered ON OFF ON OFF ON OFF, no problems",
                       [s["kind"] for s in scans] == ["ON", "OFF"] * 3 and not problems))
        # a wrong cadence is refused
        bad = tmp / "bad"
        bad.mkdir()
        for p in names[:4]:
            os.link(p, bad / p.name)
        _, prob_bad = order_cadence(bad)
        checks.append(("a 4-file directory is refused", bool(prob_bad)))
        base = ["--out-dir", str(tmp / "out"), "--snr", "8", "--workers", "1", "--nneg", "5", "--n-coarse", "2"]
        ap = build_parser()
        rc = run(ap.parse_args([str(cad), *base, "--dry-run"]))
        checks.append(("dry run returns 0 and writes nothing", rc == 0 and not (tmp / "out").exists()))
        finder = os.environ.get("FSCAD_SELFTEST_FINDER", "auto")
        rc = run(ap.parse_args([str(cad), *base, "--event-finder", finder]))
        checks.append(("full pipeline returns 0", rc == 0))
        out = tmp / "out"
        e3 = pd.read_csv(out / "events_f3.csv")
        e2 = pd.read_csv(out / "events_f2.csv")
        near = lambda df, fm: df[(df["Freq"] - fm).abs() < 0.0005] if len(df) else df  # noqa: E731
        print(f"  filter 3: {len(e3)} rows, {e3['Hit_ID'].nunique() if len(e3) else 0} events; "
              f"filter 2: {len(e2)} rows")
        # turboSETI 2.3.2 lists this event by its first-ON row only; the builtin lists every ON
        checks.append(("filter 3 keeps the ON-only sky-like tone as exactly one event",
                       len(near(e3, f["sky"])) >= 1 and len(e3) == len(near(e3, f["sky"]))
                       and e3["Hit_ID"].nunique() == 1))
        checks.append(("the all-scan RFI tone and the OFF-only tone are not events",
                       len(near(e3, f["rfi"])) == 0 and len(near(e2, f["rfi"])) == 0 and len(near(e3, f["off"])) == 0))
        checks.append(("events CSV carries the fastsearch_broadband column", "fastsearch_broadband" in e3.columns))
        fz = pd.read_csv(out / "forced_f3.csv")
        checks.append(("forced check: the sky-like tone PASSES", set(fz["verdict"]) == {"PASSES"}))
        vt = (out / "verdict_f3.txt").read_text()
        checks.append(("verdict report written, 1 event, 1 unexplained", "1 events, 1 unexplained" in vt))
        # resumable: a second run skips every stage
        t0 = time.time()
        import io
        import contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc2 = run(ap.parse_args([str(cad), *base, "--event-finder", finder]))
        s2 = buf.getvalue()
        checks.append(("second run skips the six fastsearch files and the events",
                       rc2 == 0 and s2.count("skipped") >= 6 + 3 and "[1] spliced" not in s2.replace("have", "")))
        print(f"  second (resumed) run: {time.time() - t0:.1f} s")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dir", help="directory with the six .h5 files of one cadence")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--snr", type=float, default=8.0, help="fastsearch threshold (detection significance)")
    ap.add_argument("--max-drift", type=float, default=4.0)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--coarse", help="subset of coarse channels, e.g. 150-160 (default all)")
    ap.add_argument("--n-coarse", type=int, help="coarse channels per file (default: from the header)")
    ap.add_argument("--max-hits", type=int, default=F.MAX_HITS_PER_COARSE)
    ap.add_argument("--broadband-density", type=float, default=F.BROADBAND_DENSITY)
    ap.add_argument("--broadband-keep", type=int, default=F.BROADBAND_KEEP)
    ap.add_argument("--filters", type=lambda s: [int(v) for v in s.split(",")], default=[1, 2, 3])
    ap.add_argument("--forced-filters", type=lambda s: [int(v) for v in s.split(",")], default=[2, 3])
    ap.add_argument("--forced-limit", type=int, help="forced-check at most this many rows per filter")
    ap.add_argument("--nneg", type=int, default=300, help="cadence_forced negative controls")
    ap.add_argument("--event-snr", type=float, default=None, help="find_event SNR_cut (default: none extra)")
    ap.add_argument("--event-finder", default="auto", choices=["auto", "turboseti", "builtin"])
    ap.add_argument("--rfi-db", default="", help="cadence_verdict --rfi-db argument ('' = seed, 'full', or files)")
    ap.add_argument("--redo-events", action="store_true", help="re-run stages 2-4 even if their outputs exist")
    ap.add_argument("--allow-bad-cadence", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    return ap


def main():
    if sys.argv[1:2] == ["selftest"]:
        return selftest()
    return run(build_parser().parse_args())


if __name__ == "__main__":
    sys.exit(main())
