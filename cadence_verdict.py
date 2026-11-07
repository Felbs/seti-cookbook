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
               candidate MEASURED in every scan along its predicted track. Present
               in an OFF = RFI; missing from an ON = fails the cadence; a scan that
               could not be measured leaves the event INCOMPLETE, never clean.
Anything with no explanation is reported as rung 2 ("in every ON, not in OFF,
not yet explained") - never as a detection.

    python cadence_verdict.py data/gj699_cadence/events_f3.csv
    python cadence_verdict.py events_f3.csv --forced forced_f3.csv
"""
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")


def forced_reasons(rows):
    """Reasons from cadence_forced.py rows of one event (one row per hit)."""
    out = []
    for _, r in rows.iterrows():
        tags = [c[2:] for c in r.index if c.startswith("z_")]
        snr = " ".join(f"{t[1:]}:{r['snr_' + t]:.0f}" if r['snr_' + t] == r['snr_' + t] else f"{t[1:]}:-"
                       for t in tags)
        v = r["verdict"]
        if v == "RFI":
            offs = [t for t in tags if t.endswith("OFF") and str(r["present_" + t]) in ("1", "1.0")]
            out.append(f"forced: present in OFF scan(s) {', '.join(offs)} - local (true SNR {snr})")
        elif v in ("ONE-OFF", "PARTIAL"):
            out.append(f"forced: {v} - measured absent from some ON scans (true SNR {snr})")
        elif v == "NOT-RECOVERED":
            out.append("forced: not recovered even in its own scan (not a steady narrowband track)")
        else:
            print(f"{'':14s}  forced: {v} (true SNR {snr})")
    vs = list(rows["verdict"])
    if len(vs) > 1 and out:
        if len(set(vs)) == 1:
            return [out[0] + f" [all {len(vs)} hits of this event agree]"]
        return out + [f"forced: hits of this event disagree ({', '.join(vs)})"]
    return out


def main(path, forced=None):
    import pandas as pd
    from astropy.coordinates import Angle
    import recipe_api as R
    from earth_term import earth_drift

    d = pd.read_csv(path)
    d["_row"] = range(len(d))
    # filter-2 tables leave Hit_ID empty: treat each row as its own event
    d["Hit_ID"] = [h if isinstance(h, str) else f"row{r}" for h, r in zip(d["Hit_ID"], d["_row"])]
    fz = pd.read_csv(forced) if forced else None
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
        if fz is not None:
            reasons += forced_reasons(fz[fz["row"].isin(g["_row"])])
        if not reasons:
            unexplained += 1
        tag = "explained" if reasons else "UNEXPLAINED - rung 2, try to kill it"
        print(f"{hid:12s} {f:12.6f} MHz  drift {dr:+.4f}  SNR {snr:>14s}  -> {tag}")
        for why in reasons:
            print(f"{'':14s}- {why}")
    print(f"\n{events.ngroups} events, {unexplained} unexplained")
    return 0


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("events", nargs="?", default="data/gj699_cadence/events_f3.csv")
    ap.add_argument("--forced", help="per-scan results from cadence_forced.py run (same events file)")
    a = ap.parse_args()
    sys.exit(main(a.events, a.forced))
