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
  5. forced  - optional (--forced RESULTS.csv from cadence_forced.py run): the
               candidate MEASURED in every scan along its predicted track. On the
               track in an OFF = RFI; an OFF occupied near the track = fails the
               cadence (ambiguous: same emitter or a neighbour); missing from an ON =
               fails the cadence. NOT-RECOVERED, INCOMPLETE and UNTESTED mean the
               forced check could NOT test the event: they are printed as
               information and never count as an explanation.
Anything with no explanation is reported as rung 2 ("in every ON, not in OFF,
not yet explained") - never as a detection.

    python cadence_verdict.py data/gj699_cadence/events_f3.csv
    python cadence_verdict.py events_f3.csv --forced forced_f3.csv
    python cadence_verdict.py EVENTS.csv --rfi-db [DB]   # + cross-star RFI database line
        # (rfi_db.py). No DB = the checked-in seed ONLY, which is weak (0.9 % of GJ699
        # filter-2 rows); "--rfi-db full" adds Enriquez 2017 (needs the reduced 10.8 GB
        # stream and a few GB of RAM; refused when MemAvailable is too low).
"""
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")


UNTESTABLE = ("NOT-RECOVERED", "INCOMPLETE", "UNTESTED")


def _num(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if x == x else None


def _off_detail(r, t):
    """How an OFF scan counted as present: on-track SNR, or the free-drift look's z/drift/offset."""
    how = str(r.get("how_" + t, "track"))
    how = "track" if how in ("", "nan") else how
    if how == "nearby":
        z, d, o = _num(r.get("nearz_" + t)), _num(r.get("neard_" + t)), _num(r.get("nearoff_" + t))
        dd = f" drift {d:+.3f} Hz/s" if d is not None else ""
        oo = f" {o:+.0f} Hz off the track" if o is not None else ""
        return f"{t} free-drift look z {z if z is not None else '?'}{dd}{oo}"
    snr = _num(r.get("snr_" + t))
    tag = "on-track" if how == "track" else "on-track box, other drift"
    return f"{t} {tag} SNR {snr:.0f}" if snr is not None else f"{t} {tag}"


def forced_reasons(rows):
    """(reasons, notes) from cadence_forced.py rows of one event (one row per hit).
    Only verdicts that MEASURED the event become reasons; NOT-RECOVERED,
    INCOMPLETE and UNTESTED go to notes and never explain an event."""
    out, notes = [], []
    for _, r in rows.iterrows():
        tags = [c[2:] for c in r.index if c.startswith("z_")]
        snr = " ".join(f"{t[1:]}:{_num(r['snr_' + t]):.0f}" if _num(r['snr_' + t]) is not None
                       else f"{t[1:]}:-" for t in tags)
        v = r["verdict"]
        offs = [t for t in tags if t.endswith("OFF") and str(r["present_" + t]) in ("1", "1.0")]
        if v == "RFI":
            out.append(f"forced: RFI - same track in an OFF ({'; '.join(_off_detail(r, t) for t in offs)}) "
                       f"(true SNR {snr})")
        elif v == "OFF-OCCUPIED":
            out.append(f"forced: OFF occupied near the track ({'; '.join(_off_detail(r, t) for t in offs)}) "
                       f"- fails the cadence, but same emitter or a neighbour is ambiguous (true SNR {snr})")
        elif v in ("ONE-OFF", "PARTIAL"):
            out.append(f"forced: {v} - measured absent from some ON scans (true SNR {snr})")
        elif v in UNTESTABLE:
            extra = f"; OFF power at {', '.join(offs)}" if offs else ""
            notes.append(f"forced: {v} - the forced check could not test this hit{extra} (true SNR {snr})")
        else:
            notes.append(f"forced: {v} (true SNR {snr})")
    vs = list(rows["verdict"])
    if len(vs) > 1 and out:
        if len(set(vs)) == 1:
            return [out[0] + f" [all {len(vs)} hits of this event agree]"], []
        return out + [f"forced: hits of this event disagree ({', '.join(vs)})"], notes
    return out, notes


def main(path, forced=None, rfi_db_path=None):
    import pandas as pd
    from astropy.coordinates import Angle
    import recipe_api as R
    from earth_term import earth_drift

    d = pd.read_csv(path)
    d["_row"] = range(len(d))
    # filter-2 tables leave Hit_ID empty: treat each row as its own event
    d["Hit_ID"] = [h if isinstance(h, str) else f"row{r}" for h, r in zip(d["Hit_ID"], d["_row"])]
    fz = pd.read_csv(forced) if forced else None
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
        notes = []
        if fz is not None:
            fr = fz[fz["row"].isin(g["_row"])]
            if len(fr) == 0:
                notes.append("forced: no forced-check row for this event (not tested)")
            else:
                rs, notes = forced_reasons(fr)
                reasons += rs
        if not reasons:
            unexplained += 1
        tag = "explained" if reasons else "UNEXPLAINED - rung 2, try to kill it"
        print(f"{hid:12s} {f:12.6f} MHz  drift {dr:+.4f}  SNR {snr:>14s}  -> {tag}")
        for why in reasons:
            print(f"{'':14s}- {why}")
        if db_line and not any(db_line.startswith(w[:40]) for w in reasons):
            print(f"{'':14s}  ({db_line})")
        for n in notes:
            print(f"{'':14s}  (info) {n}")
    print(f"\n{events.ngroups} events, {unexplained} unexplained")
    return 0


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("events", nargs="?", default="data/gj699_cadence/events_f3.csv")
    ap.add_argument("--forced", help="per-scan results from cadence_forced.py run (same events file)")
    ap.add_argument("--rfi-db", nargs="?", const="", default=None, metavar="DB",
                    help="add a cross-star RFI database line (rfi_db.py); no DB = the weak seed only, "
                         "'full' = seed + Enriquez 2017, or a comma-separated list of DB files")
    a = ap.parse_args()
    sys.exit(main(a.events, a.forced, a.rfi_db))
