#!/usr/bin/env python3
"""earth_term.py - the drift EARTH puts on any signal from the sky.

A transmitter on another world sends a tone; by the time it reaches a
telescope on a spinning, orbiting Earth, the tone slides in frequency. Most of
that slide is OURS: the telescope's line-of-sight acceleration from Earth's
rotation (up to ~0.034 m/s^2 at the equator) and orbit (~0.006 m/s^2). This
module predicts it for any sky position, time and site:

    v_los(t)  = (observer's barycentric velocity) . (unit vector to source)
    f_obs     = f_emit * (1 + v_los / c)
    drift     = f * d(v_los)/dt / c                 [Hz/s]

Why it matters: a real sky signal's drift = Earth term + the transmitter's own
term. Subtract the Earth term and what is left is a property of the alien
world (its spin, its orbit) - or, for local interference, nonsense. It is the
shared piece behind the hydrogen-line star-frame search, phase-connecting a
signal across ON scans, and geometry-locked drift.

GATE: Voyager 1 (Green Bank, 2016-09-19). The spacecraft's own acceleration is
~1e-7 m/s^2, so its measured drift (-0.3776 Hz/s at 8421.39 MHz, re-measured
by our own loader) must be almost entirely Earth term. We must predict it to
within 5 %.

    python earth_term.py gate
    python earth_term.py drift --ra 269.4486 --dec 4.6933 --mjd 60000.5 --f 1420.4e6
"""
import argparse
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

C = 299_792_458.0
GBT = dict(lon_deg=-79.8398, lat_deg=38.4331, height_m=807.0)   # same site as natural_signals.py

VOYAGER = dict(file="data/Voyager1.single_coarse.fine_res.h5", f_hz=8421.39e6,
               measured_drift=-0.3776)


def _setup():
    from astropy.utils import iers
    iers.conf.auto_max_age = None          # never refuse to run for stale Earth-rotation tables
    return iers


def v_los(ra_deg, dec_deg, mjd, site=GBT):
    """Observer's barycentric velocity along the line of sight, m/s (positive =
    moving toward the source). Accepts an array of MJDs."""
    _setup()
    import astropy.units as u
    from astropy.coordinates import (EarthLocation, SkyCoord,
                                     get_body_barycentric_posvel)
    from astropy.time import Time
    t = Time(np.atleast_1d(mjd), format="mjd", scale="utc")
    loc = EarthLocation.from_geodetic(lon=site["lon_deg"] * u.deg,
                                      lat=site["lat_deg"] * u.deg,
                                      height=site["height_m"] * u.m)
    _, v_site = loc.get_gcrs_posvel(t)                    # site about Earth's centre
    _, v_earth = get_body_barycentric_posvel("earth", t)  # Earth about the barycentre
    v = (v_site.xyz + v_earth.xyz).to(u.m / u.s).value    # (3, n)
    n = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs").cartesian.xyz.value
    return n @ v


def earth_drift(f_hz, ra_deg, dec_deg, mjd, site=GBT, dt_s=30.0):
    """Drift rate (Hz/s) that Earth's motion alone imposes at frequency f_hz."""
    d = dt_s / 86400.0
    vp, vm = v_los(ra_deg, dec_deg, [mjd + d, mjd - d], site)
    return f_hz * (vp - vm) / (2 * dt_s) / C


def earth_curvature(f_hz, ra_deg, dec_deg, mjd, site=GBT, dt_s=120.0):
    """Second derivative of the Earth-imposed frequency, Hz/s^2."""
    d = dt_s / 86400.0
    vp, v0, vm = v_los(ra_deg, dec_deg, [mjd + d, mjd, mjd - d], site)
    return f_hz * (vp - 2 * v0 + vm) / dt_s ** 2 / C


def cmd_gate():
    import h5py, hdf5plugin  # noqa: F401
    with h5py.File(VOYAGER["file"], "r") as f:
        a = dict(f["data"].attrs)
        nt = f["data"].shape[0]
    ra = float(a["src_raj"]) * 15.0           # header RA is in hours
    dec = float(a["src_dej"])
    mid = float(a["tstart"]) + 0.5 * nt * float(a["tsamp"]) / 86400.0
    pred = earth_drift(VOYAGER["f_hz"], ra, dec, mid)
    curv = earth_curvature(VOYAGER["f_hz"], ra, dec, mid)
    meas = VOYAGER["measured_drift"]
    err = abs(pred / meas - 1)
    ok = err < 0.05
    print(f"Voyager 1 @ RA {ra:.4f} Dec {dec:.4f}, MJD {mid:.5f}, {VOYAGER['f_hz']/1e6:.2f} MHz")
    print(f"  Earth-term drift predicted : {pred:+.4f} Hz/s")
    print(f"  drift measured in the file : {meas:+.4f} Hz/s")
    print(f"  difference                 : {100*err:.1f} %  (gate < 5 %)")
    print(f"  Earth-term curvature       : {curv:+.2e} Hz/s^2 "
          f"(bends the track {0.5*abs(curv)*1300**2:.1f} Hz over a 1300 s ON-OFF-ON-OFF-ON span)")
    print("GATE:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("gate")
    p = sub.add_parser("drift")
    p.add_argument("--ra", type=float, required=True, help="deg, ICRS")
    p.add_argument("--dec", type=float, required=True, help="deg, ICRS")
    p.add_argument("--mjd", type=float, required=True)
    p.add_argument("--f", type=float, default=1420.405751768e6, help="Hz")
    args = ap.parse_args()
    if args.cmd == "gate":
        return cmd_gate()
    print(f"Earth-term drift  {earth_drift(args.f, args.ra, args.dec, args.mjd):+.4f} Hz/s")
    print(f"Earth-term curve  {earth_curvature(args.f, args.ra, args.dec, args.mjd):+.3e} Hz/s^2")
    print(f"line-of-sight v   {v_los(args.ra, args.dec, args.mjd)[0]/1000:+.3f} km/s (toward source)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
