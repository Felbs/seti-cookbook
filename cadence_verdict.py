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
"""
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")


def main(path):
    import pandas as pd
    from astropy.coordinates import Angle
    import recipe_api as R
    from earth_term import earth_drift

    d = pd.read_csv(path)
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
        if not reasons:
            unexplained += 1
        tag = "explained" if reasons else "UNEXPLAINED - rung 2, try to kill it"
        print(f"{hid:12s} {f:12.6f} MHz  drift {dr:+.4f}  SNR {snr:>14s}  -> {tag}")
        for why in reasons:
            print(f"{'':14s}- {why}")
    print(f"\n{events.ngroups} events, {unexplained} unexplained")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "data/gj699_cadence/events_f3.csv"))
