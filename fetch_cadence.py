#!/usr/bin/env python3
"""fetch_cadence.py — download a full BL ON/OFF *cadence* for one target.

`fetch_public_data.py bl` grabs the single smallest file for a target, which is
fine for a demo and useless for the test that actually matters. Breakthrough
Listen observes ABACAD — target, throwaway sky position, target, ... — and a
signal is only interesting if it appears in the ON scans and NOT in the OFF
scans. That comparison needs the whole set.

    python fetch_cadence.py --target GJ699 --n 6

Downloads land in data/ (gitignored). Files are ~1.3 GB each, so `--n 6` is
~8 GB: check `--dry-run` first. Resumes with `curl -C -`, so an interrupted pull
costs nothing.

Credit: Breakthrough Listen open data (Lebofsky et al. 2019 PASP 131:124505;
Price et al. 2020 AJ 159:86), Berkeley SETI Research Center.
"""
import argparse
import os
import subprocess
import time
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
sys.path.insert(0, str(HERE))
import bl_archive  # noqa: E402  (exact-name, one-session, one-band cadence rules)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="GJ699")
    ap.add_argument("--n", type=int, default=6, help="how many scans of the cadence (max 6)")
    ap.add_argument("--product", default="0000",
                    choices=["0000", "0001", "0002"],
                    help="0000=fine frequency (~3 Hz, what drift "
                         "searches need), 0001=fine time, 0002=mid")
    ap.add_argument("--session", type=int, default=None,
                    help="MJD day of the session to take (default: the first valid one)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    # The archive matches target names by PREFIX (L12). bl_archive keeps only
    # exact names, real OFFs (<target>_OFF, >= 0.2 deg away), and six
    # consecutive alternating scans from ONE session, band and product.
    good = bl_archive.report(a.target, a.product)
    if a.session is not None:
        good = [c for c in good if c["day"] == a.session]
    if not good:
        print(f"no valid ON/OFF cadence for {a.target} (product {a.product}"
              f"{', session %d' % a.session if a.session else ''})", file=sys.stderr)
        return 1
    cad = good[0]
    pick = [dict(url=s["url"], name=s["name"], size=s["size"], kind=s["kind"])
            for s in cad["scans"][:a.n]]
    total = sum(f["size"] for f in pick) / 1e9
    print(f"\ntaking {len(pick)} scans of session MJD {cad['day']} "
          f"({cad['project']}, band {cad['fc']} MHz): {total:.2f} GB")
    if a.dry_run:
        return 0

    DATA.mkdir(exist_ok=True)

    # ── SINGLE-WRITER LOCK (2026-07-30) ─────────────────────────────────────
    # Learned the expensive way: two copies of this script were started by
    # accident (a background launch that looked like it had failed, then a
    # relaunch). Both ran, both `curl -C -` onto the SAME filenames, and resume
    # semantics turned two writers into interleaved garbage — all six files came
    # out 158 MB to 1.4 GB OVERSIZED and unreadable past the second time row.
    # 80 GB of download thrown away for want of six lines.
    lock = DATA / ".fetch_cadence.lock"
    if lock.exists():
        try:
            age = time.time() - lock.stat().st_mtime
        except OSError:
            age = 0
        if age < 6 * 3600:
            print(f"another fetch_cadence appears to be running "
                  f"({lock}, {age/60:.0f} min old). Refusing to run a second "
                  f"writer — delete the lock if you are sure.", file=sys.stderr)
            return 2
    lock.write_text(str(os.getpid()), encoding="utf-8")

    got, bad = [], []
    try:
        for f in pick:
            out = DATA / f["name"]
            if out.exists():
                have = out.stat().st_size
                # EXACT size check, both directions. An oversized file is
                # corrupt (this is what a double-writer produces) and resuming
                # onto it appends yet more garbage, so it is deleted, not resumed.
                if have == f["size"]:
                    print(f"[have] {f['name']}", flush=True)
                    got.append(out)
                    continue
                if have > f["size"]:
                    print(f"[bad ] {f['name']} is {(have-f['size'])/1e6:+.0f} MB "
                          f"OVERSIZED — corrupt, deleting and refetching",
                          flush=True)
                    out.unlink()
                else:
                    print(f"[part] {f['name']} {have/1e9:.2f}/"
                          f"{f['size']/1e9:.2f} GB — resuming", flush=True)
            print(f"[get ] {f['name']}  ({f['size']/1e9:.2f} GB)", flush=True)
            rc = subprocess.run(["curl", "-L", "-C", "-", "--retry", "3",
                                 "--retry-delay", "5", "-o", str(out),
                                 f["url"]]).returncode
            size_ok = out.exists() and out.stat().st_size == f["size"]
            if rc == 0 and size_ok:
                got.append(out)
            else:
                actual = out.stat().st_size if out.exists() else 0
                print(f"[fail] {f['name']}: curl rc={rc}, size {actual/1e9:.2f} "
                      f"vs {f['size']/1e9:.2f} GB", flush=True)
                bad.append(f["name"])
    finally:
        lock.unlink(missing_ok=True)

    print(f"\ndownloaded {len(got)}/{len(pick)} into {DATA}")
    for g in got:
        print(f"  {g.name}  {g.stat().st_size/1e9:.2f} GB  OK")
    for b in bad:
        print(f"  {b}  FAILED")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
