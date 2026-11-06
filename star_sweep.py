#!/usr/bin/env python3
"""star_sweep.py - autonomous SETI sweep over a batch of nearby stars.

For each target: query the Breakthrough Listen open archive, pull the smallest
GBT HDF5, run turboSETI's drift search, then apply the forensics that matter -
is a hit a real candidate or RFI? The RFI tell we learned on Epsilon Eridani:
hits that all share ONE drift rate are instrumental, not a planet's signal.
A target is only "interesting" if it has hits with VARIED drift rates that also
survive the sideband/structure checks. Everything is explained, not thresholded.

  python star_sweep.py                      # default nearby-star list
  python star_sweep.py GJ699 HIP54035 ...   # custom targets

Downloads land in data/ (gitignored); a summary goes to SWEEP_RESULTS.md.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
API = "http://seti.berkeley.edu/opendata/api/query-files"
PY = sys.executable

# nearby / classic SETI targets (BL archive target names)
DEFAULT = ["GJ699", "GJ411", "GJ887", "GJ273", "GJ15A", "GJ71",
           "HIP54035", "HIP57548", "GJ876", "GJ581"]


def query_smallest(target):
    out = DATA / f"_{target}.json"
    try:
        subprocess.run(["curl", "-s", "--connect-timeout", "25",
                        f"{API}?target={target}&telescopes=GBT&file-types=HDF5&limit=12",
                        "-o", str(out)], timeout=40)
        rows = json.loads(out.read_text()).get("data", [])
    except Exception as e:
        return None, f"query error {e}"
    rows = [r for r in rows if r.get("size", 0) > 1e6]
    # turboSETI needs FINE-FREQUENCY products (.0000 hi-spectral, .0002 mid);
    # the .0001 high-time-resolution files are rejected (issue #117).
    fine = [r for r in rows if ".0001." not in r["url"]
            and (".0000." in r["url"] or ".0002." in r["url"] or "gpuspec" in r["url"])]
    fine = [r for r in fine if ".0001." not in r["url"]]
    pool = fine or [r for r in rows if ".0001." not in r["url"]]
    if not pool:
        return None, "no fine-frequency GBT HDF5 (only HTR)"
    return sorted(pool, key=lambda r: r["size"])[0], None


def download(row, target):
    fn = DATA / f"star_{target}.h5"
    if fn.exists() and fn.stat().st_size > 1e7:
        return str(fn)
    r = subprocess.run(["curl", "-sL", "--connect-timeout", "40", "-o", str(fn),
                        row["url"]], timeout=900)
    return str(fn) if fn.exists() and fn.stat().st_size > 1e7 else None


def hunt(h5):
    """Drift search + forensics. Returns dict(hits, unique_drifts, verdict)."""
    import os
    from turbo_seti.find_doppler.find_doppler import FindDoppler
    dat = h5[:-3] + ".dat"
    # turboSETI on Windows needs a RELATIVE path from inside the data dir
    # (absolute out_dir -> [Errno 22]); run from DATA with basenames.
    cwd = os.getcwd(); base = os.path.basename(h5)
    try:
        os.chdir(str(DATA))
        FindDoppler(base, max_drift=8.0, snr=8, out_dir=".",
                    append_output=False).search()
    except Exception as e:
        os.chdir(cwd)
        return {"verdict": f"search error {str(e)[:40]}"}
    os.chdir(cwd)
    hits = []
    p = Path(dat)
    if p.exists():
        for ln in p.read_text().splitlines():
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.split()
            hits.append((float(f[1]), float(f[2]), float(f[3])))  # drift, snr, mhz
    drifts = sorted({round(d, 3) for d, _, _ in hits})
    n = len(hits)
    # Can this file resolve drift at all? On BL mid-res (.0002) products one
    # drift step is ~9.8 Hz/s - wider than a planet's whole range (~0.05-4
    # Hz/s) - so EVERY hit reports the same drift and "one shared drift rate =
    # RFI" is void. drift_audit.py demonstrates it with planted signals.
    import h5py, hdf5plugin  # noqa: F401
    with h5py.File(h5, "r") as f:
        a = dict(f["data"].attrs)
        nt = f["data"].shape[0]
    quantum = abs(a["foff"]) * 1e6 / ((nt - 1) * a["tsamp"])
    if quantum > 0.5:
        verdict = (f"UNRESOLVED: {n} hits, but this file's drift step is {quantum:.1f} Hz/s "
                   f"- too coarse to tell a planet from RFI. Screening only; "
                   f"re-run on fine-resolution (.0000) ON/OFF data")
        return {"hits": n, "unique_drifts": len(drifts), "verdict": verdict}
    # forensics: uniform drift => RFI; varied drift among strong hits => look closer
    if n == 0:
        verdict = "CLEAN (0 hits) - quiet sky"
    elif len(drifts) <= 2:
        verdict = f"RFI ({n} hits all at {len(drifts)} drift rate(s) = instrumental)"
    else:
        # a real candidate would have varied drift AND not be a giant RFI blast
        strong = [h for h in hits if h[1] > 20]
        verdict = (f"REVIEW: {n} hits, {len(drifts)} distinct drift rates "
                   f"({len(strong)} strong) - run sideband_pairs + find_event")
    return {"hits": n, "unique_drifts": len(drifts), "verdict": verdict}


def main():
    targets = sys.argv[1:] or [t for t in DEFAULT if t.isascii()]
    DATA.mkdir(exist_ok=True)
    results = []
    for tgt in targets:
        print(f"\n=== {tgt} ===", flush=True)
        row, err = query_smallest(tgt)
        if not row:
            print(f"  skip: {err}", flush=True)
            results.append((tgt, None, err))
            continue
        print(f"  archive: {row['size']/1e9:.2f} GB @ {row.get('center_freq'):.0f} MHz - downloading", flush=True)
        h5 = download(row, tgt)
        if not h5:
            print("  download failed", flush=True)
            results.append((tgt, None, "download failed"))
            continue
        r = hunt(h5)
        print(f"  -> {r['verdict']}", flush=True)
        results.append((tgt, r, None))
        try:
            Path(h5).unlink()          # reclaim disk after the hunt
            Path(h5[:-3] + ".dat").unlink()
        except Exception:
            pass
    # summary
    lines = ["# Nearby-star SETI sweep\n",
             f"Swept {len(results)} targets (GBT open data, turboSETI + forensics).\n",
             "| target | hits | drift rates | verdict |", "|---|---|---|---|"]
    for tgt, r, err in results:
        if r:
            lines.append(f"| {tgt} | {r.get('hits','-')} | {r.get('unique_drifts','-')} | {r['verdict']} |")
        else:
            lines.append(f"| {tgt} | - | - | {err} |")
    reviews = [t for t, r, e in results if r and r["verdict"].startswith("REVIEW")]
    lines.append(f"\n**Candidates needing review: {reviews or 'none'}**")
    (HERE / "SWEEP_RESULTS.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n" + "\n".join(lines[4:]), flush=True)


if __name__ == "__main__":
    main()
