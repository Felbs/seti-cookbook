#!/usr/bin/env python3
"""cadence_verdict.py - explain every event turboSETI's ON/OFF filter leaves.

find_event_pipeline keeps signals seen in the ON scans and not the OFFs. That
is necessary, not sufficient: a satellite parked in the geostationary belt
near the target, or interference that happens to fade during the OFFs,
passes it. This tool takes the surviving events and, for each, asks:

  1. drift   - does it drift like a signal from the TARGET? A sky source carries
               at least Earth's term (earth_term.py) for that position, time
               and frequency. Drift within one quantum of zero while Earth's
               term is well away from zero = locked to Earth.
  2. band    - is it inside a known interference allocation (recipe_api.RFI_BANDS)?
  3. nature  - is it on a natural line (HI, OH, masers)?
  4. number  - is it on a round-number frequency?
Anything with no explanation is reported as rung 2 ("in every ON, not in OFF,
not yet explained") - never as a detection.

    python cadence_verdict.py data/gj699_cadence/events_f3.csv
    python cadence_verdict.py EVENTS.csv --rfi-db [DB]   # + cross-star RFI database line
        # (rfi_db.py). No DB = the checked-in seed ONLY, which is weak (0.9 % of GJ699
        # filter-2 rows); "--rfi-db full" adds Enriquez 2017 (needs the reduced 10.8 GB
        # stream and a few GB of RAM; refused when MemAvailable is too low).
"""
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")


def main(path, rfi_db_path=None):
    import pandas as pd
    from astropy.coordinates import Angle
    import recipe_api as R
    from earth_term import earth_drift

    d = pd.read_csv(path)
    db = None
    if rfi_db_path is not None:
        import rfi_db
        db = rfi_db.default_db(rfi_db_path or None)
        print(rfi_db.describe(db) + "; the event's own target is always left out")
    events = d.groupby("Hit_ID")
    print(f"{path}: {len(d)} rows = {events.ngroups} events "
          f"(each seen in {d.groupby('Hit_ID').size().max()} ON scans)\n")
    unexplained = 0
    for hid, g in events:
        r0 = g.iloc[0]
        f = float(g["Freq"].mean())
        dr = float(g["DriftRate"].median())
        snr = ", ".join(f"{s:.0f}" for s in g["SNR"])
        ra = Angle(r0["RA"]).deg
        dec = Angle(r0["DEC"]).deg
        mjd = float(g["MJD"].mean())
        e_term = earth_drift(f * 1e6, ra, dec, mjd)
        quantum = abs(float(r0["DELTAF"])) / float(r0["obs_length"])
        reasons = []
        if abs(dr) <= 1.5 * quantum and abs(e_term) > 3 * quantum:
            reasons.append(f"zero drift ({dr:+.4f} Hz/s) where a sky source must show Earth's "
                           f"{e_term:+.3f} Hz/s - locked to Earth")
        for fn in (R.rfi_reason, R.natural_reason, R.round_number_reason):
            why = fn(f)
            if why:
                reasons.append(why)
        db_line = None
        if db is not None:
            v = db.query(f, dr, mjd=mjd, target=str(r0["Source"]),
                         tobs_s=float(r0["obs_length"]), res_hz=abs(float(r0["DELTAF"])))
            why, db_line = rfi_db.explain_line(v)
            if why:
                reasons.append(why)
        if not reasons:
            unexplained += 1
        tag = "explained" if reasons else "UNEXPLAINED - rung 2, try to kill it"
        print(f"{hid:12s} {f:12.6f} MHz  drift {dr:+.4f}  SNR {snr:>14s}  -> {tag}")
        for why in reasons:
            print(f"{'':14s}- {why}")
        if db_line and not any(db_line.startswith(w[:40]) for w in reasons):
            print(f"{'':14s}  ({db_line})")
    print(f"\n{events.ngroups} events, {unexplained} unexplained")
    return 0


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="explain every event turboSETI's ON/OFF filter leaves")
    ap.add_argument("events", nargs="?", default="data/gj699_cadence/events_f3.csv")
    ap.add_argument("--rfi-db", nargs="?", const="", default=None, metavar="DB",
                    help="add a cross-star RFI database line (rfi_db.py); no DB = the weak seed only, "
                         "'full' = seed + Enriquez 2017, or a comma-separated list of DB files")
    a = ap.parse_args()
    sys.exit(main(a.events, a.rfi_db))
