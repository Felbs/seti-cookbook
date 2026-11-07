#!/usr/bin/env python3
"""xepoch.py - is the same tone there on another DATE, in the star's own frame?

A beacon held at a fixed frequency in its home star's rest frame keeps that
frequency from one observing session to the next, once Earth's motion and the
star's own motion are taken out. Local interference does not: it is fixed in
OUR (topocentric) frame, so in the star's frame it moves by Earth's change of
line-of-sight velocity between the dates (km/s, i.e. kHz at 1.4 GHz). This
tool turns every hit of every scan into the STAR frame and asks which tones
come back on different dates.

THE FRAME (per scan, at MID-scan time, never the file start):
    f_star = f_obs * gamma_s (1 + rv(t)/c) / (gamma_o (1 + v_los(t)/c))
  v_los   observer's barycentric velocity toward the star (earth_term.py's
          Voyager-gated ephemeris and Green Bank site), projected on the star's
          direction AT THAT DATE: the position is propagated with proper motion,
          parallax and RV (astropy apply_space_motion). For Barnard's Star the
          J2000 direction is 2.8' off by 2016, which mis-projects Earth's
          velocity by up to ~24 m/s (~110 Hz at 1.42 GHz).
  rv(t)   the star's radial velocity AT THAT DATE: includes the secular
          (perspective) acceleration mu^2 d. Barnard's Star: +4.55 m/s/yr =
          21.5 Hz/yr at 1420 MHz (selftest checks it against mu^2 d).
  gammas  special-relativistic transverse terms (Earth ~7 Hz, Barnard ~64 Hz at
          1.42 GHz, nearly constant between dates). The Sun's gravitational
          redshift at Earth (varies ~0.3 Hz over a year at 1.42 GHz) is NOT
          modelled; it is far below the tolerance.
  drift   the Earth term is removed BEFORE drifts are compared:
          drift_star = (drift_obs - f_obs * a_los(t)/c) * f_star/f_obs.

MATCHING (no fixed buckets: a pair 3-9 Hz apart that straddles a bucket edge is
not lost). ON hits that survive the vetoes are clustered within each epoch
(same star-frame frequency within --tol in >= --min-ons ON scans), then
epoch clusters are matched across epochs by a sorted merge: every cluster is
an anchor, and the nearest cluster of each other epoch within +-tol joins it.
An EPOCH is one session: scans less than 6 h apart. Never a scan index.

VETOES
  topocentric OFF veto  an OFF hit of the same epoch within --tol-topo of the ON
                        hit's observed mid-scan frequency, either where it sits
                        or where its drift would carry it by the OFF's time
                        (RFI is fixed, or drifts on, in OUR frame).
  star-frame OFF veto   an OFF hit of ANY epoch within --tol-star-off (default
                        3 x tol) in the star frame (OFF hits converted with the
                        ON star's direction): the tone is there when we look away.
  drift consistency     the members' Earth-removed drifts agree within --drift-tol.
  degenerate pair       two dates whose star-frame/topocentric offset differs by
                        less than 2 x tol-topo cannot tell a star from RFI; such
                        pairs are reported and never count as confirmation.

FALSE ALARMS. With n_i surviving clusters in epoch i spread over the common
searched bandwidth B (star frame), the chance that a random 2*tol cell holds
one is p_i = n_i * 2 tol / B; the expected number of chance k-of-N coincidences
is (B / 2 tol) * P(>= k of N), P from the Poisson-binomial of the p_i. (Not the
number of OCCUPIED cells, which is meaningless as a denominator.) The selftest
checks this against Monte Carlo.

WHAT IT TESTS: a transmitter whose frequency is stabilised in its star's rest
frame (a deliberate beacon; intrinsic drift ~0). A transmitter on a planet that
does not correct for its own orbit and spin moves by kHz-MHz between dates and
is NOT tested here.

    python xepoch.py selftest                # synthetic epochs + setigen end-to-end, offline (~1 min)
    python xepoch.py inject                  # real GJ699 noise, three PSEUDO-epochs (needs data/hi_*GJ699*.npz)
    python xepoch.py run --target GJ699 A.dat B.dat ...   # real turboSETI/fastsearch .dat files
    python xepoch.py fetch --star GJ725 --days 57587 57667 57804   # one coarse channel per scan, HTTP range reads
    python xepoch.py real --star GJ725       # search the fetched coarse channels and match across dates
"""
import argparse
import json
import math
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
DATA = HERE / "data"
C = 299_792_458.0
F_HI = 1420.405751768e6
SESSION_GAP_D = 0.25            # scans more than 6 h apart are different epochs

DEFAULTS = dict(tol_hz=10.0, tol_topo_hz=500.0, tol_star_off_hz=30.0, drift_tol=0.05,
                min_ons=2, min_epochs=2)


def _setup():
    from astropy.utils import iers
    iers.conf.auto_max_age = None
    iers.conf.auto_download = False      # historical dates: the bundled IERS-B table is enough
    return iers


# ================================================================ star model
class Star:
    """A star's J2000/ICRS astrometry and space motion. pm in mas/yr (pmra is
    mu_alpha* = cos(dec) d(alpha)/dt), parallax in mas, rv in km/s."""

    def __init__(self, name, ra, dec, pmra=None, pmdec=None, plx=None, rv=None):
        self.name, self.ra, self.dec = name, float(ra), float(dec)
        self.pmra, self.pmdec, self.plx = pmra, pmdec, plx
        self.rv = 0.0 if rv is None else float(rv)
        self.moving = pmra is not None and pmdec is not None and plx is not None and plx > 0

    @classmethod
    def from_targets(cls, name):
        d = json.loads((HERE / "targets.json").read_text())
        for s in d["stars"]:
            if name in (s["star"], s.get("bl_token")) or name.upper() in (s["star"].upper(),
                                                                          (s.get("bl_token") or "").upper()):
                rv = s.get("rv_km_s")
                if rv is not None and abs(rv) > 600:
                    raise SystemExit(f"{name}: SIMBAD RV {rv} km/s is not a stellar RV; pass --rv")
                return cls(s["star"], s["ra"], s["dec"], s.get("pmra_mas_yr"), s.get("pmdec_mas_yr"),
                           s.get("parallax_mas"), rv)
        raise SystemExit(f"{name} not in targets.json")

    def state(self, mjd):
        """Unit vectors (3, n), radial velocity (n,) m/s and total speed (n,) m/s at each MJD."""
        _setup()
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from astropy.time import Time
        mjd = np.atleast_1d(np.asarray(mjd, dtype=float))
        if self.moving:
            c = SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg, distance=(1000.0 / self.plx) * u.pc,
                         pm_ra_cosdec=self.pmra * u.mas / u.yr, pm_dec=self.pmdec * u.mas / u.yr,
                         radial_velocity=self.rv * u.km / u.s, obstime=Time("J2000"))
            n = c.apply_space_motion(new_obstime=Time(mjd, format="mjd"))
            xyz = n.cartesian.xyz.value
            xyz = xyz / np.linalg.norm(xyz, axis=0)
            rv = np.atleast_1d(n.radial_velocity.to(u.m / u.s).value) * np.ones(mjd.size)
            mu = math.radians(math.hypot(self.pmra, self.pmdec) / 3.6e6) / (365.25 * 86400)
            vt = mu * (1000.0 / self.plx) * 3.0856775814913673e16
            speed = np.sqrt(rv ** 2 + vt ** 2)
        else:
            c = SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg, frame="icrs")
            xyz = np.repeat(c.cartesian.xyz.value[:, None], mjd.size, axis=1)
            rv = np.full(mjd.size, self.rv * 1e3)
            speed = np.abs(rv)
        return xyz, rv, speed

    def secular_accel(self):
        """mu^2 d in m/s^2 (0 without space motion)."""
        if not self.moving:
            return 0.0
        mu = math.radians(math.hypot(self.pmra, self.pmdec) / 3.6e6) / (365.25 * 86400)
        return mu ** 2 * (1000.0 / self.plx) * 3.0856775814913673e16


def observer_velocity(mjd, site=None):
    """Green Bank's barycentric velocity (3, n) m/s - the same calculation as earth_term.v_los."""
    _setup()
    import astropy.units as u
    from astropy.coordinates import EarthLocation, get_body_barycentric_posvel
    from astropy.time import Time
    from earth_term import GBT
    site = site or GBT
    t = Time(np.atleast_1d(np.asarray(mjd, dtype=float)), format="mjd", scale="utc")
    loc = EarthLocation.from_geodetic(lon=site["lon_deg"] * u.deg, lat=site["lat_deg"] * u.deg,
                                      height=site["height_m"] * u.m)
    _, v_site = loc.get_gcrs_posvel(t)
    _, v_earth = get_body_barycentric_posvel("earth", t)
    return (v_site.xyz + v_earth.xyz).to(u.m / u.s).value


def frame(star, mjd, site=None, model=True, dt_s=30.0):
    """For each MJD: K (f_star = K f_obs), a_los (m/s^2, d v_los/dt with the
    direction held fixed), v_los (m/s toward the star) and rv (m/s).
    model=False: J2000 direction and constant RV (what a naive pipeline does)."""
    mjd = np.atleast_1d(np.asarray(mjd, dtype=float))
    if model:
        n, rv, speed = star.state(mjd)
    else:
        n, rv, speed = Star(star.name, star.ra, star.dec, rv=star.rv).state(mjd)
    d = dt_s / 86400.0
    v0 = observer_velocity(mjd, site)
    vp = observer_velocity(mjd + d, site)
    vm = observer_velocity(mjd - d, site)
    vlos = np.einsum("ij,ij->j", n, v0)
    alos = (np.einsum("ij,ij->j", n, vp) - np.einsum("ij,ij->j", n, vm)) / (2 * dt_s)
    g_o = 1.0 / np.sqrt(1.0 - (v0 ** 2).sum(0) / C ** 2)
    g_s = 1.0 / np.sqrt(1.0 - speed ** 2 / C ** 2)
    K = g_s * (1.0 + rv / C) / (g_o * (1.0 + vlos / C))
    return K, alos, vlos, rv


# ================================================================ hits
def read_dat(path, target=None):
    """A turboSETI or fastsearch .dat -> scan dict with hit arrays. kind is ON if
    the header's Source equals target (case-insensitive), else OFF."""
    hdr, f0, dr, snr = {}, [], [], []
    for ln in open(path):
        if ln.startswith("#"):
            for k, pat in (("source", r"Source:\s*(\S+)"), ("mjd", r"MJD:\s*([0-9.]+)"),
                           ("tsamp", r"DELTAT:\s*([0-9.eE+-]+)"), ("foff_hz", r"DELTAF\(Hz\):\s*([0-9.eE+-]+)"),
                           ("obs_length", r"obs_length:\s*([0-9.eE+-]+)"), ("ra", r"RA:\s*(\S+)"),
                           ("dec", r"DEC:\s*(\S+)")):
                m = re.search(pat, ln)
                if m and k not in hdr:
                    hdr[k] = m.group(1)
            continue
        p = ln.split()
        if len(p) < 4:
            continue
        dr.append(float(p[1]))
        snr.append(float(p[2]))
        f0.append(float(p[3]) * 1e6)
    tsamp = float(hdr["tsamp"])
    nt = max(1, int(round(float(hdr["obs_length"]) / tsamp)))
    src = hdr.get("source", "?")
    kind = "ON" if target and src.upper() == target.upper() else "OFF"
    return make_scan(Path(path).name, kind, float(hdr["mjd"]), tsamp, nt, float(hdr["foff_hz"]),
                     f0, dr, snr, source=src)


def make_scan(name, kind, tstart_mjd, tsamp, nt, foff_hz, f0_hz, drift, snr, coverage_hz=None, source=None):
    """f0_hz: hit frequency at the FIRST spectrum; drift Hz/s (observed)."""
    f0 = np.asarray(f0_hz, dtype=float)
    dr = np.asarray(drift, dtype=float)
    span = (nt - 1) * tsamp
    return dict(name=name, kind=kind, source=source, tstart=float(tstart_mjd), tsamp=float(tsamp), nt=int(nt),
                foff_hz=float(foff_hz), t_mid=float(tstart_mjd) + nt * tsamp / 2 / 86400.0,
                f_mid=f0 + dr * span / 2, drift=dr, snr=np.asarray(snr, dtype=float),
                coverage=None if coverage_hz is None else tuple(sorted(coverage_hz)))


def group_epochs(scans, gap_d=SESSION_GAP_D):
    """Assign scan['epoch'] = 'MJD<day>' by session: a new epoch starts after a gap > gap_d."""
    order = sorted(range(len(scans)), key=lambda i: scans[i]["t_mid"])
    labels, last, cur = {}, None, None
    for i in order:
        t = scans[i]["t_mid"]
        if last is None or t - last > gap_d:
            cur = f"MJD{int(math.floor(scans[i]['tstart']))}"
            while cur in labels.values() and any(labels[j] == cur and abs(scans[j]["t_mid"] - t) > gap_d
                                                 for j in labels):
                cur += "b"
        labels[i] = cur
        last = t
    for i, s in enumerate(scans):
        s["epoch"] = labels[i]
    return sorted(set(labels.values()), key=lambda e: min(scans[i]["t_mid"] for i in labels if labels[i] == e))


def to_star_frame(star, scans, site=None, model=True):
    for s in scans:
        K, a, v, rv = frame(star, [s["t_mid"]], site, model)
        s["K"], s["a_los"], s["v_los"], s["rv"] = float(K[0]), float(a[0]), float(v[0]), float(rv[0])
        s["f_star"] = s["f_mid"] * s["K"]
        s["drift_star"] = (s["drift"] - s["f_mid"] * s["a_los"] / C) * s["K"]
    return scans


def _near(sorted_vals, x, tol):
    """Boolean: some value of sorted_vals within tol of each x."""
    if sorted_vals.size == 0:
        return np.zeros(np.shape(x), bool)
    j = np.searchsorted(sorted_vals, x)
    lo = np.abs(x - sorted_vals[np.clip(j - 1, 0, sorted_vals.size - 1)])
    hi = np.abs(sorted_vals[np.clip(j, 0, sorted_vals.size - 1)] - x)
    return np.minimum(lo, hi) <= tol


def cross_merge(cf, ce, tol_hz, min_epochs):
    """Sorted merge across epochs. cf: sorted star-frame frequencies of the
    epoch clusters, ce: their epoch labels. Every cluster is an anchor; the
    nearest cluster of each other epoch within +-tol joins it. Returns member
    index sets (frozensets) with >= min_epochs epochs, strict subsets dropped."""
    groups, seen = [], set()
    for i in range(cf.size):
        a, b = np.searchsorted(cf, cf[i] - tol_hz), np.searchsorted(cf, cf[i] + tol_hz, side="right")
        pick = {ce[i]: i}
        for j in range(a, b):
            e = ce[j]
            if e == ce[i]:
                continue
            if e not in pick or abs(cf[j] - cf[i]) < abs(cf[pick[e]] - cf[i]):
                pick[e] = j
        if len(pick) < min_epochs:
            continue
        key = frozenset(pick.values())
        if key not in seen:
            seen.add(key)
            groups.append(key)
    big = [g for g in groups if len(g) > min_epochs]
    return [g for g in groups if not any(g < h for h in big)]


_GEOM = {}


def geom_factor(m, n=400_000):
    """How much larger the anchor rule's acceptance volume is than (2 tol)^(m-1):
    m points form a group when SOME member has all others within +-tol. Exact
    values: m=2 -> 1, m=3 -> 1.5; others by a fixed-seed Monte Carlo."""
    if m <= 2:
        return 1.0
    if m not in _GEOM:
        r = np.random.default_rng(12345)
        x = np.concatenate([np.zeros((n, 1)), r.uniform(-2, 2, (n, m - 1))], axis=1)
        ok = np.zeros(n, bool)
        for i in range(m):
            ok |= np.all(np.abs(x - x[:, i:i + 1]) <= 1.0, axis=1)
        _GEOM[m] = float(ok.mean() * 4.0 ** (m - 1) / 2.0 ** (m - 1))
    return _GEOM[m]


def chance_groups(p, k, n_cells):
    """Expected chance groups spanning >= k epochs: n_cells * sum_m P(exactly m) * geom(m)."""
    dist = np.zeros(len(p) + 1)
    dist[0] = 1.0
    for pi in p:
        dist[1:] = dist[1:] * (1 - pi) + dist[:-1] * pi
        dist[0] *= (1 - pi)
    return float(n_cells * sum(dist[m] * geom_factor(m) for m in range(k, len(p) + 1)))


def shift_test(cf, ce, cd, lo, hi, tol_hz, drift_tol, min_epochs, n_epochs, rounds=40, seed=3):
    """Empirical chance rate: shift every epoch's clusters by its own random
    offset (cyclic inside the common coverage), keep each epoch's real
    structure (combs, crowding), and count groups again. Returns the mean
    number of groups with >= k epochs, before and after the drift test."""
    if rounds <= 0 or cf.size == 0 or hi <= lo:
        return None
    rng = np.random.default_rng(seed)
    labs = sorted(set(ce.tolist()))
    acc = {k: [] for k in range(min_epochs, n_epochs + 1)}
    acc_d = {k: [] for k in range(min_epochs, n_epochs + 1)}
    for _ in range(rounds):
        f = cf.copy()
        for e in labs:
            m = ce == e
            f[m] = lo + np.mod(f[m] - lo + rng.uniform(0, hi - lo), hi - lo)
        o = np.argsort(f)
        fs, es, ds = f[o], ce[o], cd[o]
        g = cross_merge(fs, es, tol_hz, min_epochs)
        for k in acc:
            gk = [x for x in g if len(x) >= k]
            acc[k].append(len(gk))
            acc_d[k].append(sum(1 for x in gk if np.ptp(ds[list(x)]) <= drift_tol))
    return dict(rounds=rounds, mean_groups={str(k): float(np.mean(v)) for k, v in acc.items()},
                mean_groups_drift_ok={str(k): float(np.mean(v)) for k, v in acc_d.items()})


# ================================================================ analysis
def analyse(star, scans, tol_hz=10.0, tol_topo_hz=500.0, tol_star_off_hz=None, drift_tol=0.05,
            min_ons=2, min_epochs=2, band_hz=None, site=None, model=True, verbose=True, frame_fn=None,
            shift_rounds=40):
    """frame_fn(star, scans) replaces the physical frame (selftest only).
    tol_star_off_hz defaults to 3 x tol_hz."""
    tol_star_off_hz = 3 * tol_hz if tol_star_off_hz is None else tol_star_off_hz
    epochs = group_epochs(scans)
    (frame_fn or (lambda st, sc: to_star_frame(st, sc, site, model)))(star, scans)
    ons = [s for s in scans if s["kind"] == "ON"]
    offs = [s for s in scans if s["kind"] == "OFF"]
    # ---- common searched coverage, star frame
    cov_note = "per-scan coverage"
    for s in scans:
        if s["coverage"] is None and band_hz is not None:
            s["coverage"] = tuple(sorted(band_hz))
    if any(s["coverage"] is None for s in ons):
        allf = np.concatenate([s["f_mid"] for s in scans if s["f_mid"].size]) if scans else np.array([0.0])
        for s in scans:
            if s["coverage"] is None:
                s["coverage"] = (float(allf.min()), float(allf.max()))
        cov_note = "NO coverage given: hit span used (B and false-alarm numbers are approximate)"
    ep_cov = {}
    for e in epochs:
        lo = max(s["coverage"][0] * s["K"] for s in ons if s["epoch"] == e) if any(
            s["epoch"] == e for s in ons) else None
        hi = min(s["coverage"][1] * s["K"] for s in ons if s["epoch"] == e) if lo is not None else None
        ep_cov[e] = (lo, hi)
    on_epochs = [e for e in epochs if ep_cov[e][0] is not None]
    if len(on_epochs) < 2 and verbose:
        print(f"note: only {len(on_epochs)} epoch (session) with ON scans - cross-epoch matching needs two dates; "
              "the vetoes and clusters below are single-session only")
    lo = max(ep_cov[e][0] for e in on_epochs)
    hi = min(ep_cov[e][1] for e in on_epochs)
    if hi <= lo:
        raise SystemExit("the epochs share no star-frame frequency coverage: refusing to match")
    B = hi - lo
    # ---- vetoes
    off_star = np.sort(np.concatenate([s["f_star"] for s in offs])) if offs else np.array([])
    rows = []
    stats = dict(on_hits=0, topo_vetoed=0, star_off_vetoed=0, outside_coverage=0, survivors=0)
    for s in ons:
        n = s["f_mid"].size
        stats["on_hits"] += n
        if n == 0:
            continue
        topo = np.zeros(n, bool)
        for o in offs:
            if o["epoch"] != s["epoch"] or o["f_mid"].size == 0:
                continue
            of = np.sort(o["f_mid"])
            dt = (o["t_mid"] - s["t_mid"]) * 86400.0
            topo |= _near(of, s["f_mid"], tol_topo_hz)
            topo |= _near(of, s["f_mid"] + s["drift"] * dt, tol_topo_hz)
        staroff = _near(off_star, s["f_star"], tol_star_off_hz) & ~topo
        inside = (s["f_star"] >= lo) & (s["f_star"] <= hi)
        keep = ~topo & ~staroff & inside
        stats["topo_vetoed"] += int(topo.sum())
        stats["star_off_vetoed"] += int(staroff.sum())
        stats["outside_coverage"] += int((~topo & ~staroff & ~inside).sum())
        stats["survivors"] += int(keep.sum())
        for i in np.flatnonzero(keep):
            rows.append(dict(epoch=s["epoch"], scan=s["name"], t_mid=s["t_mid"], f_mid=float(s["f_mid"][i]),
                             drift=float(s["drift"][i]), snr=float(s["snr"][i]), f_star=float(s["f_star"][i]),
                             drift_star=float(s["drift_star"][i])))
    # chance that a tone that IS unique is vetoed anyway (crowded OFFs): random frequencies,
    # with drifts drawn from the scan's own hits, through the same two vetoes
    rng_v = np.random.default_rng(1)
    fv = []
    for s in ons:
        fr_ = rng_v.uniform(s["coverage"][0], s["coverage"][1], 2000)
        dr_ = rng_v.choice(s["drift"], 2000) if s["drift"].size else np.zeros(2000)
        v = np.zeros(2000, bool)
        for o in offs:
            if o["epoch"] != s["epoch"] or o["f_mid"].size == 0:
                continue
            of = np.sort(o["f_mid"])
            v |= _near(of, fr_, tol_topo_hz) | _near(of, fr_ + dr_ * (o["t_mid"] - s["t_mid"]) * 86400.0, tol_topo_hz)
        v |= _near(off_star, fr_ * s["K"], tol_star_off_hz)
        fv.append(float(v.mean()))
    stats["chance_veto_per_on_hit"] = float(np.mean(fv)) if fv else 0.0
    # ---- per-epoch clusters (same star-frame frequency in >= min_ons ON scans)
    clusters = []
    for e in epochs:
        r = sorted([x for x in rows if x["epoch"] == e], key=lambda x: x["f_star"])
        if not r:
            continue
        f = np.array([x["f_star"] for x in r])
        used = np.zeros(len(r), bool)
        for i in sorted(range(len(r)), key=lambda i: -r[i]["snr"]):
            if used[i]:
                continue
            a, b = np.searchsorted(f, f[i] - tol_hz), np.searchsorted(f, f[i] + tol_hz, side="right")
            mem = [j for j in range(a, b) if not used[j]]
            used[mem] = True
            if len({r[j]["scan"] for j in mem}) < min_ons:
                continue
            w = np.array([max(r[j]["snr"], 1e-3) for j in mem])
            clusters.append(dict(epoch=e, f_star=float(np.average(f[mem], weights=w)),
                                 drift_star=float(np.median([r[j]["drift_star"] for j in mem])),
                                 snr=float(w.max()), n_ons=len({r[j]["scan"] for j in mem}),
                                 members=[r[j] for j in mem]))
    n_by_ep = {e: sum(1 for c in clusters if c["epoch"] == e) for e in on_epochs}
    # ---- cross-epoch sorted merge
    clusters.sort(key=lambda c: c["f_star"])
    cf = np.array([c["f_star"] for c in clusters])
    ce = np.array([c["epoch"] for c in clusters])
    groups = cross_merge(cf, ce, tol_hz, min_epochs)
    # ---- per-group tests
    K_ep = {e: np.mean([s["K"] for s in ons if s["epoch"] == e]) for e in on_epochs}
    out = []
    for g in groups:
        mem = sorted((clusters[j] for j in g), key=lambda c: c["epoch"])
        fs = [c["f_star"] for c in mem]
        ds = [c["drift_star"] for c in mem]
        f_obs = [np.mean([m["f_mid"] for m in c["members"]]) for c in mem]
        eps = [c["epoch"] for c in mem]
        f_ref = np.mean(fs)
        degenerate = []
        for x in range(len(eps)):
            for y in range(x + 1, len(eps)):
                if abs(K_ep[eps[x]] - K_ep[eps[y]]) * f_ref / max(K_ep[eps[x]], 1e-9) < 2 * tol_topo_hz:
                    degenerate.append((eps[x], eps[y]))
        topo_same = (max(f_obs) - min(f_obs)) <= tol_topo_hz
        drift_ok = (max(ds) - min(ds)) <= drift_tol
        verdict = ("PASS" if drift_ok and not degenerate else
                   "FAIL drift" if not drift_ok else "UNTESTABLE (dates too alike in Earth velocity)")
        out.append(dict(f_star_mhz=f_ref / 1e6, spread_hz=max(fs) - min(fs), n_epochs=len(eps), epochs=eps,
                        drift_star=ds, drift_spread=max(ds) - min(ds), drift_ok=drift_ok,
                        topo_freqs_mhz=[x / 1e6 for x in f_obs], same_topocentric=topo_same,
                        degenerate_pairs=degenerate, verdict=verdict,
                        snr=[c["snr"] for c in mem], n_ons=[c["n_ons"] for c in mem]))
    out.sort(key=lambda g: (g["verdict"] != "PASS", -g["n_epochs"], -min(g["snr"])))
    # ---- false-alarm estimate: analytic (occupancy) and empirical (shift test)
    p = [min(1.0, n_by_ep[e] * 2 * tol_hz / B) for e in on_epochs]
    n_cells = B / (2 * tol_hz)
    fap = {k: chance_groups(p, k, n_cells) for k in range(2, len(on_epochs) + 1)}
    shift = shift_test(cf, ce, np.array([c["drift_star"] for c in clusters]), lo, hi, tol_hz, drift_tol,
                       min_epochs, len(on_epochs), rounds=shift_rounds)
    # drift-consistency survival of random pairs (empirical, from the clusters' own drifts)
    q = None
    if len(on_epochs) >= 2:
        d_by = [np.array([c["drift_star"] for c in clusters if c["epoch"] == e]) for e in on_epochs]
        prs = []
        rng = np.random.default_rng(0)
        for x in range(len(d_by)):
            for y in range(x + 1, len(d_by)):
                if d_by[x].size and d_by[y].size:
                    a_ = rng.choice(d_by[x], 4000)
                    b_ = rng.choice(d_by[y], 4000)
                    prs.append(np.mean(np.abs(a_ - b_) <= drift_tol))
        q = float(np.mean(prs)) if prs else None
    res = dict(target=star.name, params=dict(tol_hz=tol_hz, tol_topo_hz=tol_topo_hz,
                                              tol_star_off_hz=tol_star_off_hz, drift_tol=drift_tol,
                                              min_ons=min_ons, min_epochs=min_epochs,
                                              star_model="propagated position + rv(t) incl. secular"
                                              if model else "J2000 direction, constant RV"),
               epochs=[dict(epoch=e, scans=[s["name"] for s in scans if s["epoch"] == e],
                            n_on=sum(1 for s in ons if s["epoch"] == e),
                            n_off=sum(1 for s in offs if s["epoch"] == e),
                            v_los_kms=float(np.mean([s["v_los"] for s in ons if s["epoch"] == e])) / 1e3
                            if e in on_epochs else None,
                            K_minus_1_ppm=(K_ep[e] - 1) * 1e6 if e in on_epochs else None,
                            clusters=n_by_ep.get(e, 0)) for e in epochs],
               coverage_star_mhz=[lo / 1e6, hi / 1e6], B_hz=B, coverage_note=cov_note, counts=stats,
               trials=dict(cells=n_cells, p_i=p),
               false_alarm=dict(expected_chance_groups={str(k): v for k, v in fap.items()},
                                shift_test=shift,
                                drift_pair_survival=q,
                                expected_after_drift={str(k): (v * q ** (k - 1) if q is not None else None)
                                                      for k, v in fap.items()}),
               groups=out, n_pass=sum(1 for g in out if g["verdict"] == "PASS"), _survivors=rows)
    if verbose:
        report(res)
    return res


def report(res, n=12):
    p = res["params"]
    print(f"{res['target']}: {len(res['epochs'])} epochs; tol {p['tol_hz']} Hz star frame, topo veto "
          f"{p['tol_topo_hz']} Hz, star-frame OFF veto {p['tol_star_off_hz']} Hz, drift tol {p['drift_tol']} Hz/s, "
          f">= {p['min_ons']} ONs per epoch, >= {p['min_epochs']} epochs; model: {p['star_model']}")
    for e in res["epochs"]:
        print(f"  {e['epoch']}: {e['n_on']} ON + {e['n_off']} OFF, v_los {e['v_los_kms'] or 0:+.3f} km/s, "
              f"{e['clusters']} surviving clusters")
    c = res["counts"]
    print(f"  ON hits {c['on_hits']}: topo-vetoed {c['topo_vetoed']}, star-OFF-vetoed {c['star_off_vetoed']}, "
          f"outside common coverage {c['outside_coverage']}, survivors {c['survivors']}; a unique tone is "
          f"vetoed by chance with p = {c['chance_veto_per_on_hit']:.3f} per ON hit")
    print(f"  common coverage {res['coverage_star_mhz'][0]:.6f}-{res['coverage_star_mhz'][1]:.6f} MHz (star frame), "
          f"B = {res['B_hz']/1e3:.1f} kHz  [{res['coverage_note']}]")
    fa = res["false_alarm"]
    print("  expected chance groups: " + ", ".join(f">={k} epochs {v:.3g}" for k, v in
                                                  fa["expected_chance_groups"].items())
          + (f"; x drift-pair survival {fa['drift_pair_survival']:.3f}" if fa["drift_pair_survival"] is not None
             else ""))
    print(f"  groups found {len(res['groups'])}, PASS {res['n_pass']}")
    for g in res["groups"][:n]:
        print(f"    {g['verdict']:<8} {g['f_star_mhz']:.6f} MHz  {g['n_epochs']} epochs, spread {g['spread_hz']:.1f} Hz, "
              f"drift spread {g['drift_spread']:.3f}, snr {', '.join(f'{x:.0f}' for x in g['snr'])}")


# ================================================================ data helpers
def load_npz(path):
    z = np.load(path)
    m = {k: (z[k].item() if z[k].ndim == 0 else z[k]) for k in z.files if k != "data"}
    nf = z["data"].shape[1]
    m["freqs_hz"] = (m["fch1"] + m["foff"] * (m["ch0"] + np.arange(nf))) * 1e6
    return z["data"], m


def search_scan(x, freqs_hz, foff_hz, tsamp, max_drift=0.5, thr=10.0, dc_idx=()):
    """fastsearch.search_slab on one waterfall -> (f0_hz, drift, snr) arrays."""
    import fastsearch as F
    hits, info = F.search_slab(x, foff_hz, tsamp, max_drift=max_drift, thr=thr, dc_idx=dc_idx)
    ch = np.array([h["ch"] for h in hits])
    f0 = np.interp(ch, np.arange(freqs_hz.size), freqs_hz) if ch.size else np.array([])
    return f0, np.array([h["drift"] for h in hits]), np.array([h["snr"] for h in hits]), info


def plant_track(star, t_mid, tstart, tsamp, nt, freqs_hz, f_star, drift_intr=0.0, site=None):
    """Where a star-frame tone sits in one scan: (channel at t=0 (float), observed drift)."""
    K, a, _, _ = frame(star, [t_mid], site)
    f_mid = f_star / K[0]
    drift = f_mid * a[0] / C + drift_intr / K[0]
    f0 = f_mid - drift * (nt - 1) * tsamp / 2
    df = freqs_hz[1] - freqs_hz[0]
    return (f0 - freqs_hz[0]) / df, drift


# ================================================================ selftest
def _mc_check(rng, B=2e6, n=(3000, 2500, 4000), tol=10.0, trials=60):
    """Monte Carlo: uniform random cluster frequencies through the SAME matcher
    -> mean chance groups, vs chance_groups()."""
    counts = {2: [], 3: []}
    for _ in range(trials):
        f = np.concatenate([rng.uniform(0, B, k) for k in n])
        e = np.concatenate([np.full(k, i) for i, k in enumerate(n)])
        o = np.argsort(f)
        g = cross_merge(f[o], e[o], tol, 2)
        counts[2].append(sum(1 for x in g if len(x) >= 2))
        counts[3].append(sum(1 for x in g if len(x) >= 3))
    p = [k * 2 * tol / B for k in n]
    pred = {k: chance_groups(p, k, B / (2 * tol)) for k in (2, 3)}
    return {k: float(np.mean(v)) for k, v in counts.items()}, pred


def selftest():
    t_start = time.time()
    rng = np.random.default_rng(7)
    checks = []
    barnard = Star("GJ699", 269.4520769586187, 4.693364966576667, -801.551, 10362.394, 546.9759, -110.11)
    # 1. secular acceleration: astropy's propagated RV vs mu^2 d
    _, rv, _ = barnard.state([51544.5, 51544.5 + 365.25])
    drv = rv[1] - rv[0]
    sec = barnard.secular_accel() * 365.25 * 86400
    checks.append((f"Barnard's secular RV change {drv:.3f} m/s/yr = mu^2 d {sec:.3f} "
                   f"({F_HI*drv/C:.2f} Hz/yr at 1420 MHz)", abs(drv - sec) < 0.02 and 21.0 < F_HI * drv / C < 22.0))
    # 2. our v_los equals earth_term's for a fixed direction
    import earth_term as E
    star_fixed = Star("X", 83.6, 22.0)
    K, a, v, _ = frame(star_fixed, [57500.3])
    v_et = E.v_los(83.6, 22.0, 57500.3)[0]
    d_et = E.earth_drift(1.42e9, 83.6, 22.0, 57500.3)
    checks.append(("v_los and Earth drift identical to earth_term (Voyager-gated)",
                   abs(v[0] - v_et) < 1e-3 and abs(1.42e9 * a[0] / C - d_et) < 1e-6))
    # 3. star-frame matching of a beacon two years apart: needs the full model
    days = [57457.56, 57457.56 + 200.0, 57457.56 + 700.0]     # not whole years: Earth's velocity must differ
    f_star0 = F_HI
    scans = []
    for k, t in enumerate(days):
        Kt, at, _, _ = frame(barnard, [t])
        fm = f_star0 / Kt[0]
        def away(n_):                    # random lines, kept 2 kHz from the beacon (chance vetoes are
            v = rng.uniform(fm - 5e5, fm + 5e5, 4 * n_)    # measured by `inject`, not here)
            return v[np.abs(v - fm) > 2e3][:n_]
        noise = away(200)
        for j, kind in enumerate(["ON", "OFF", "ON", "OFF", "ON", "OFF"]):
            tm = t + j * 5.5 / 1440
            Kj, aj, _, _ = frame(barnard, [tm])
            fmj = f_star0 / Kj[0]
            dj = fmj * aj[0] / C
            f0 = [fmj - dj * 15 * 18.25 / 2] if kind == "ON" else []
            fr = list(f0) + list(away(60))
            dr = [dj] * len(f0) + list(rng.uniform(-0.5, 0.5, 60))
            fr += list(noise[:40])                       # RFI: in ON and OFF at the same topocentric freq
            dr += [0.0] * 40
            scans.append(make_scan(f"e{k}s{j}", kind, tm - 16 * 18.25 / 2 / 86400, 18.25, 16, -2.79,
                                   fr, dr, [30.0] * len(fr), coverage_hz=(fm - 6e5, fm + 6e5)))
    res = analyse(barnard, scans, min_epochs=3, verbose=False)
    ok_full = any(g["verdict"] == "PASS" and abs(g["f_star_mhz"] * 1e6 - f_star0) < 10 and g["n_epochs"] == 3
                  for g in res["groups"])
    res_n = analyse(barnard, [dict(s) for s in scans], min_epochs=2, model=False, verbose=False)
    lost_naive = not any(g["n_epochs"] >= 2 for g in res_n["groups"])
    checks.append(("beacon on 3 dates over 1.9 yr matched in the star frame (10 Hz tol)", ok_full))
    checks.append(("...and LOST with J2000 direction + constant RV (the naive frame)", lost_naive))
    checks.append(("RFI in ON and OFF never survives", res["counts"]["topo_vetoed"] >= 3 * 3 * 40))
    # 4. bucket-edge: two tones 6 Hz apart straddling a 10 Hz bucket boundary are matched
    sc = [make_scan("a", "ON", 57000.0, 18.25, 16, -2.79, [1.42e9 + 9.0e3 - 3], [0.0], [20], (1.419e9, 1.421e9)),
          make_scan("a2", "ON", 57000.002, 18.25, 16, -2.79, [1.42e9 + 9.0e3 - 3], [0.0], [20], (1.419e9, 1.421e9)),
          make_scan("b", "ON", 57100.0, 18.25, 16, -2.79, [1.42e9 + 9.0e3 + 3], [0.0], [20], (1.419e9, 1.421e9)),
          make_scan("b2", "ON", 57100.002, 18.25, 16, -2.79, [1.42e9 + 9.0e3 + 3], [0.0], [20], (1.419e9, 1.421e9))]
    unit = Star("U", 0.0, 0.0)

    def flat(star, s_list):              # identity frame: K = 1, no Earth term
        for s in s_list:
            s["K"], s["a_los"], s["v_los"], s["rv"] = 1.0, 0.0, 0.0, 0.0
            s["f_star"], s["drift_star"] = s["f_mid"].copy(), s["drift"].copy()
        return s_list
    r = analyse(unit, sc, tol_hz=10.0, min_ons=2, verbose=False, frame_fn=flat)
    r_epochs = len(r["epochs"])
    # 5. scans of one session are ONE epoch, never two
    same_day = [make_scan(f"s{i}", "ON", 57000.0 + i * 0.004, 18.25, 16, -2.79, [1.42e9], [0.0], [20],
                          (1.419e9, 1.421e9)) for i in range(3)]
    r1 = analyse(unit, same_day, min_ons=1, verbose=False, frame_fn=flat)
    checks.append(("tones 6 Hz apart across a 10-Hz bucket edge are matched (no buckets)",
                   len(r["groups"]) == 1 and r["groups"][0]["spread_hz"] < 7))
    checks.append(("an epoch is a session: 3 scans 6 min apart = 1 epoch, no self-match",
                   len(r1["epochs"]) == 1 and r1["groups"] == [] and r_epochs == 2))
    # 6. false-alarm formula vs Monte Carlo
    mc, pred = _mc_check(rng)
    ok_mc = abs(mc[2] / pred[2] - 1) < 0.05 and abs(mc[3] / pred[3] - 1) < 0.2
    checks.append((f"chance groups, Monte Carlo vs formula: >=2 {mc[2]:.1f} vs {pred[2]:.1f}, "
                   f">=3 {mc[3]:.3f} vs {pred[3]:.3f}", ok_mc))
    # 7. end-to-end on synthetic noise with setigen plants (fastsearch, numpy ok)
    e2e = synthetic_injection(rng)
    checks.append((f"setigen end-to-end: star-frame tone recovered on all 3 dates "
                   f"({e2e['recovered']}/{e2e['planted']})", e2e["recovered"] == e2e["planted"]))
    checks.append((f"topocentric tone in ON+OFF vetoed ({e2e['topo_vetoed']}/{e2e['topo_planted']}), "
                   f"ON-only topocentric tone not matched across dates ({e2e['ononly_matched']} matched)",
                   e2e["topo_vetoed"] == e2e["topo_planted"] and e2e["ononly_matched"] == 0))
    checks.append((f"no other PASS group in noise ({e2e['other_pass']})", e2e["other_pass"] == 0))
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print(f"selftest {'PASS' if ok else 'FAIL'} ({time.time()-t_start:.0f} s)")
    return 0 if ok else 1


def synthetic_injection(rng, nf=1 << 16, days=(57457.55, 57470.40, 57489.60)):
    """Three synthetic dates (1 ON + 1 OFF each), chi-square noise, setigen plants:
    a star-frame tone (all ONs), a topocentric tone (all ONs and OFFs) and an
    ON-only topocentric tone. Barnard's Star geometry."""
    import inject
    barnard = Star("GJ699", 269.4520769586187, 4.693364966576667, -801.551, 10362.394, 546.9759, -110.11)
    tsamp, nt, foff = 18.253611008, 16, -2.7939677238464355
    f_star0 = 1421.05e6
    Ks = [frame(barnard, [t])[0][0] for t in days]
    f_mids = [f_star0 / k for k in Ks]
    top = 0.5 * (max(f_mids) + min(f_mids)) + 0.5 * nf * abs(foff)    # shared topocentric coverage
    freqs = top + foff * np.arange(nf)
    f_topo = top + foff * nf * 0.25          # topocentric tone (ON+OFF)
    f_onl = top + foff * nf * 0.80           # ON-only topocentric tone
    scans = []
    for k, t in enumerate(days):
        for j, kind in enumerate(("ON", "OFF")):
            tm = t + j * 5.5 / 1440
            ts = tm - nt * tsamp / 2 / 86400
            x = rng.chisquare(64, (nt, nf)).astype(np.float32)
            plants = []
            if kind == "ON":
                ch, dr = plant_track(barnard, tm, ts, tsamp, nt, freqs, f_star0)
                plants.append(dict(ch=ch, drift=dr, snr=40))
                plants.append(dict(ch=(f_onl - freqs[0]) / foff, drift=0.0, snr=40))
            plants.append(dict(ch=(f_topo - freqs[0]) / foff, drift=-0.05, snr=40))
            inject.plant(x, foff, tsamp, plants)
            f0, d, s, _ = search_scan(x, freqs, foff, tsamp, max_drift=0.4, thr=10.0)
            scans.append(make_scan(f"d{k}{kind}", kind, ts, tsamp, nt, foff, f0, d, s,
                                   coverage_hz=(freqs.min(), freqs.max())))
    res = analyse(barnard, scans, min_ons=1, min_epochs=3, verbose=False)
    rec = sum(1 for g in res["groups"] if g["verdict"] == "PASS" and abs(g["f_star_mhz"] * 1e6 - f_star0) < 10)
    on_scans = [s for s in scans if s["kind"] == "ON"]
    f_topo_mid = f_topo - 0.05 * (nt - 1) * tsamp / 2
    on_only = sum(1 for g in res["groups"] if any(abs(f * 1e6 - f_onl) < 100 for f in g["topo_freqs_mhz"]))
    other = sum(1 for g in res["groups"] if g["verdict"] == "PASS" and abs(g["f_star_mhz"] * 1e6 - f_star0) >= 10)
    # topocentric tone: found in the ON scan, and NOT among the survivors of the vetoes
    vetoed = 0
    for s in on_scans:
        found = np.any(np.abs(s["f_mid"] - f_topo_mid) < 50)
        left = any(x["scan"] == s["name"] and abs(x["f_mid"] - f_topo_mid) < 50 for x in res["_survivors"])
        vetoed += int(found and not left)
    return dict(planted=1, recovered=min(rec, 1), topo_planted=len(on_scans), topo_vetoed=vetoed,
                ononly_matched=on_only, other_pass=other, result=res)


# ================================================================ real data: GJ699 pseudo-epochs
PSEUDO_SHIFT_D = [0.0, 200.0, 400.0]


def cmd_inject(args):
    """Real GJ699 noise and RFI, three PSEUDO-epochs: (ON 0003, OFF 0004), (ON 0005,
    OFF 0006), (ON 0007, OFF 0008) relabelled +0, +200, +400 days. The data are
    one session; the dates are labels that give the frame physics a year to act
    over. Plants: star-frame tones (all three ONs), topocentric tones (all ONs and
    OFFs), ON-only topocentric tones."""
    import fastsearch_recovery as FR
    files = sorted(DATA.glob("hi_*GJ699*.npz"))
    if len(files) < 6:
        raise SystemExit("needs the six data/hi_*GJ699*.npz slices (python hi_starframe.py fetch)")
    star = Star.from_targets("GJ699")
    segs = []
    for p in files:
        d, m = load_npz(p)
        seg, lo, dc = FR.coarse_segment(np.asarray(d, dtype=np.float32), m)
        fr = m["freqs_hz"][lo:lo + seg.shape[1]]
        kind = "OFF" if "_OFF_" in p.name else "ON"
        segs.append(dict(name=p.name, kind=kind, x=seg, freqs=fr, dc=dc, tstart=float(m["tstart"]),
                         tsamp=float(m["tsamp"]), foff=float(m["foff"]) * 1e6))
        del d
    nt = segs[0]["x"].shape[0]
    pairs = [(0, 1), (2, 3), (4, 5)]
    for (a, b), sh in zip(pairs, PSEUDO_SHIFT_D):
        for i in (a, b):
            segs[i]["tstart_label"] = segs[i]["tstart"] + sh
            segs[i]["t_mid_label"] = segs[i]["tstart_label"] + nt * segs[i]["tsamp"] / 2 / 86400
    print(f"GJ699 pseudo-epochs: {len(segs)} scans, {segs[0]['x'].shape[1]} channels each "
          f"({segs[0]['freqs'].min()/1e6:.4f}-{segs[0]['freqs'].max()/1e6:.4f} MHz), labels +0/+200/+400 d")
    return inject_run(star, segs, args, min_ons=1,
                      kind="pseudo-epochs: one real GJ699 session (2016-03-10), relabelled +0/+200/+400 d")


def inject_run(star, segs, args, min_ons, kind):
    """Plant star-frame tones (every ON), topocentric tones (every ON and OFF) and
    ON-only topocentric tones into real waterfalls `segs` (dicts with x, freqs,
    foff, tsamp, dc, kind, tstart_label, t_mid_label), search each with
    fastsearch, and score the whole cross-epoch chain at each SNR level."""
    import inject
    rng = np.random.default_rng(args.seed)
    nt = segs[0]["x"].shape[0]
    # null search (no plants): the data's own lines, to keep plants away from them
    null = []
    t0 = time.time()
    for s in segs:
        f0, d, sn, info = search_scan(s["x"], s["freqs"], s["foff"], s["tsamp"], args.max_drift, args.thr,
                                      dc_idx=[s["dc"]] if s["dc"] is not None else ())
        null.append(make_scan(s["name"], s["kind"], s["tstart_label"], s["tsamp"], nt, s["foff"], f0, d, sn,
                              coverage_hz=(s["freqs"].min(), s["freqs"].max())))
        print(f"  null {s['kind']:3s} {s['name'][-40:]}: {f0.size} hits ({time.time()-t0:.0f} s)")
    res0 = analyse(star, null, min_ons=min_ons, min_epochs=2, verbose=True)
    n_ep = len({x["epoch"] for x in null if x["kind"] == "ON"})
    n_on = sum(1 for x in segs if x["kind"] == "ON")
    # plant positions: star-frame tones whose three topocentric positions avoid every null hit
    lo_f = max(s["freqs"].min() for s in segs) + 30e3
    hi_f = min(s["freqs"].max() for s in segs) - 30e3
    allhits = np.sort(np.concatenate([s["f_mid"] for s in null]))
    Ks = [frame(star, [x["t_mid_label"]])[0][0] for x in segs if x["kind"] == "ON"]
    n_star, n_topo, n_onl = args.n_star, args.n_topo, args.n_topo
    star_f, topo_f, onl_f = [], [], []
    tries = 0
    while len(star_f) < n_star and tries < 20000:
        tries += 1
        fs = rng.uniform(lo_f * max(Ks), hi_f * min(Ks))
        tops = [fs / k for k in Ks]
        if any(not (lo_f < t < hi_f) for t in tops) or any(_near(allhits, np.array([t]), 2e3)[0] for t in tops):
            continue
        if any(abs(fs - x) < 8e3 for x in star_f) or any(abs(t - u_) < 8e3 for t in tops for u_ in topo_f + onl_f):
            continue
        star_f.append(fs)
    for lst, n in ((topo_f, n_topo), (onl_f, n_onl)):
        while len(lst) < n and tries < 40000:
            tries += 1
            ft = rng.uniform(lo_f, hi_f)
            if _near(allhits, np.array([ft]), 2e3)[0]:
                continue
            if any(abs(ft - x) < 8e3 for x in topo_f + onl_f) or any(abs(ft - fs / k) < 8e3
                                                                       for fs in star_f for k in Ks):
                continue
            lst.append(ft)
    levels = args.snr
    out = dict(kind=kind, n_epochs=n_ep, min_ons=min_ons,
               files=[s["name"] for s in segs], null=dict(groups=len(res0["groups"]), n_pass=res0["n_pass"],
                                                          expected=res0["false_alarm"],
                                                          counts=res0["counts"]),
               n_star=len(star_f), n_topo=len(topo_f), n_ononly=len(onl_f), levels={})
    for snr in levels:
        scans = []
        for s in segs:
            x = s["x"].astype(np.float32, copy=True)
            df = s["freqs"][1] - s["freqs"][0]
            pl = []
            if s["kind"] == "ON":
                for fs in star_f:
                    ch, dr = plant_track(star, s["t_mid_label"], s["tstart_label"], s["tsamp"], nt, s["freqs"], fs)
                    pl.append(dict(ch=ch, drift=dr, snr=snr))
                for ft in onl_f:
                    pl.append(dict(ch=(ft - s["freqs"][0]) / df, drift=-0.02, snr=snr))
            for ft in topo_f:
                pl.append(dict(ch=(ft - s["freqs"][0]) / df, drift=-0.02, snr=snr))
            inject.plant(x, s["foff"], s["tsamp"], pl)
            f0, d, sn, _ = search_scan(x, s["freqs"], s["foff"], s["tsamp"], args.max_drift, args.thr,
                                       dc_idx=[s["dc"]] if s["dc"] is not None else ())
            scans.append(make_scan(s["name"], s["kind"], s["tstart_label"], s["tsamp"], nt, s["foff"], f0, d, sn,
                                   coverage_hz=(s["freqs"].min(), s["freqs"].max())))
            del x
        r = analyse(star, scans, min_ons=min_ons, min_epochs=args.min_epochs, verbose=False)
        passes = [g for g in r["groups"] if g["verdict"] == "PASS"]
        rec = sum(1 for fs in star_f if any(abs(g["f_star_mhz"] * 1e6 - fs) <= r["params"]["tol_hz"]
                                            and g["n_epochs"] == n_ep for g in passes))
        rec2 = sum(1 for fs in star_f if any(abs(g["f_star_mhz"] * 1e6 - fs) <= r["params"]["tol_hz"]
                                             for g in passes))
        # each star-frame plant found per ON scan (before matching), to separate search from matching
        found = 0
        for sc in scans:
            if sc["kind"] != "ON":
                continue
            for fs in star_f:
                found += int(np.any(np.abs(sc["f_star"] - fs) <= 10.0))
        topo_survive = sum(1 for g in r["groups"] if any(abs(f * 1e6 - ft) < 300 for f in g["topo_freqs_mhz"]
                                                         for ft in topo_f))
        onl_match = sum(1 for g in r["groups"] if any(abs(f * 1e6 - ft) < 300 for f in g["topo_freqs_mhz"]
                                                      for ft in onl_f))
        other = [g for g in passes if not any(abs(g["f_star_mhz"] * 1e6 - fs) <= 10 for fs in star_f)]
        out["levels"][str(snr)] = dict(recovered_all=rec, recovered_any=rec2, of=len(star_f),
                                       found_per_on_scan=found, of_on=n_on * len(star_f),
                                       topo_groups=topo_survive, ononly_groups=onl_match,
                                       other_pass=len(other), counts=r["counts"],
                                       expected=r["false_alarm"]["expected_chance_groups"])
        print(f"  true SNR {snr:>4}: star-frame tones matched on all {n_ep} dates {rec}/{len(star_f)} "
              f"(>= {args.min_epochs}: {rec2}); found per ON scan {found}/{n_on*len(star_f)}; "
              f"topocentric groups {topo_survive}, ON-only topocentric groups {onl_match}, other PASS {len(other)}"
              f"  ({time.time()-t0:.0f} s)")
    if getattr(args, "json", None) and args.cmd == "inject":
        Path(args.json).write_text(json.dumps(out, indent=1, default=float))
        print(f"wrote {args.json}")
    return out


# ================================================================ real data: other epochs
def census_row(star, band="L"):
    d = json.loads((HERE / "EPOCH_CENSUS.json").read_text())
    for r in d["rows"]:
        if (r["star"] == star or r["token"] == star) and r["band"] == band:
            return r
    raise SystemExit(f"{star} {band}: not in EPOCH_CENSUS.json")


def cmd_fetch(args):
    """One coarse channel (the chunk holding --at-mhz, default the star-frame HI
    line) of every scan of the chosen sessions, by HTTP range read."""
    import fsspec
    import h5py
    import hdf5plugin  # noqa: F401
    row = census_row(args.star, args.band)
    star = Star.from_targets(args.star)
    out_dir = DATA / f"xepoch_{row['token']}"
    out_dir.mkdir(parents=True, exist_ok=True)
    sess = [s for s in row["sessions"] if s["mjd_day"] in args.days] if args.days else row["sessions"]
    total = 0
    for s in sess:
        if s.get("pointing_ok") is False:
            print(f"skip MJD {s['mjd_day']}: pointing {s['pointing_offset_arcmin']}' off the star")
            continue
        for name, kind in zip(s["files"], s["kinds"]):
            out = out_dir / name.replace(".h5", ".npz")
            if out.exists():
                print(f"have {out.name}")
                continue
            url = s["url_dir"] + "/" + name
            with fsspec.open(url, block_size=8 * 2 ** 20).open() as fh:
                with h5py.File(fh, "r") as h:
                    dset = h["data"]
                    a = dict(dset.attrs)
                    tm = float(a["tstart"]) + dset.shape[0] * float(a["tsamp"]) / 2 / 86400
                    at = args.at_mhz * 1e6 if args.at_mhz else F_HI / frame(star, [tm])[0][0]
                    fpc = dset.chunks[2]
                    i = int((at / 1e6 - a["fch1"]) / a["foff"])
                    c0 = (i // fpc) * fpc
                    x = dset[:, 0, c0:c0 + fpc].astype(np.float32)
            np.savez(out, data=x, ch0=c0, fch1=a["fch1"], foff=a["foff"], tsamp=a["tsamp"], tstart=a["tstart"],
                     ra_h=a["src_raj"], dec_deg=a["src_dej"], kind="ON" if kind == "A" else "OFF",
                     source=str(a.get("source_name", "")), nfpc=fpc)
            total += x.nbytes
            print(f"{'ON ' if kind == 'A' else 'OFF'} {name[-48:]}: {x.nbytes/1e6:.0f} MB, ch {c0}.. "
                  f"(total {total/1e9:.2f} GB)")
            time.sleep(0.5)


def cmd_real(args):
    star = Star.from_targets(args.star)
    d = DATA / f"xepoch_{args.star}"
    files = sorted(d.glob("*.npz"))
    if not files:
        raise SystemExit(f"nothing in {d}: run `xepoch.py fetch --star {args.star}` first")
    if args.inject:
        segs = []
        for p in files:
            x, m = load_npz(p)
            nf = x.shape[1]
            segs.append(dict(name=p.name, kind=str(m["kind"]), x=np.asarray(x, np.float32), freqs=m["freqs_hz"],
                             dc=nf // 2 if int(m.get("nfpc", nf)) == nf else None, tstart=float(m["tstart"]),
                             tstart_label=float(m["tstart"]), tsamp=float(m["tsamp"]), foff=float(m["foff"]) * 1e6,
                             t_mid_label=float(m["tstart"]) + x.shape[0] * float(m["tsamp"]) / 2 / 86400))
        print(f"{args.star}: {len(segs)} real scans, injection through the full chain")
        out = inject_run(star, segs, args, min_ons=args.min_ons, kind=f"{args.star}: real dates")
        if args.json:
            Path(args.json).write_text(json.dumps(out, indent=1, default=float))
        return out
    scans = []
    t0 = time.time()
    for p in files:
        x, m = load_npz(p)
        nf = x.shape[1]
        dc = [nf // 2] if int(m.get("nfpc", nf)) == nf else []
        f0, dr, sn, info = search_scan(np.asarray(x, np.float32), m["freqs_hz"], float(m["foff"]) * 1e6,
                                       float(m["tsamp"]), args.max_drift, args.thr, dc_idx=dc)
        kind = str(m["kind"])
        scans.append(make_scan(p.name, kind, float(m["tstart"]), float(m["tsamp"]), x.shape[0],
                               float(m["foff"]) * 1e6, f0, dr, sn,
                               coverage_hz=(m["freqs_hz"].min(), m["freqs_hz"].max())))
        print(f"  {kind:3s} {p.name[-52:]}: {f0.size} hits ({time.time()-t0:.0f} s)")
        del x
    res = analyse(star, scans, min_ons=args.min_ons, min_epochs=args.min_epochs)
    if args.json:
        Path(args.json).write_text(json.dumps({k: v for k, v in res.items() if not k.startswith("_")}, indent=1, default=float))
        print(f"wrote {args.json}")


def cmd_run(args):
    star = Star.from_targets(args.target) if args.rv is None else Star(
        args.target, args.ra, args.dec, args.pmra, args.pmdec, args.plx, args.rv)
    scans = [read_dat(p, args.target) for p in args.dat]
    band = tuple(x * 1e6 for x in args.band) if args.band else None
    res = analyse(star, scans, args.tol, args.tol_topo, args.tol_star_off, args.drift_tol, args.min_ons,
                  args.min_epochs, band)
    if args.json:
        Path(args.json).write_text(json.dumps({k: v for k, v in res.items() if not k.startswith("_")}, indent=1, default=float))
        print(f"wrote {args.json}")


def main():
    if sys.argv[1:2] in (["selftest"], []):
        return selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run")
    r.add_argument("dat", nargs="+")
    r.add_argument("--target", required=True)
    for k, v in (("--tol", 10.0), ("--tol-topo", 500.0), ("--tol-star-off", 30.0), ("--drift-tol", 0.05)):
        r.add_argument(k, type=float, default=v)
    r.add_argument("--min-ons", type=int, default=2)
    r.add_argument("--min-epochs", type=int, default=2)
    r.add_argument("--band", type=float, nargs=2, metavar=("LO_MHZ", "HI_MHZ"),
                   help="searched band of every scan (topocentric); without it the hit span is used")
    for k in ("--ra", "--dec", "--pmra", "--pmdec", "--plx", "--rv"):
        r.add_argument(k, type=float)
    r.add_argument("--json")
    i = sp.add_parser("inject")
    i.add_argument("--seed", type=int, default=1)
    i.add_argument("--snr", type=float, nargs="+", default=[10, 12, 15, 20, 30])
    i.add_argument("--n-star", type=int, default=20)
    i.add_argument("--n-topo", type=int, default=10)
    i.add_argument("--thr", type=float, default=10.0)
    i.add_argument("--max-drift", type=float, default=0.5)
    i.add_argument("--min-epochs", type=int, default=2)
    i.add_argument("--json")
    f = sp.add_parser("fetch")
    f.add_argument("--star", required=True)
    f.add_argument("--band", default="L")
    f.add_argument("--days", type=int, nargs="*")
    f.add_argument("--at-mhz", type=float, help="default: the star-frame HI line")
    e = sp.add_parser("real")
    e.add_argument("--star", required=True)
    e.add_argument("--thr", type=float, default=10.0)
    e.add_argument("--max-drift", type=float, default=1.0)
    e.add_argument("--min-ons", type=int, default=2)
    e.add_argument("--min-epochs", type=int, default=2)
    e.add_argument("--json")
    e.add_argument("--inject", action="store_true", help="plant star-frame and topocentric tones, score the chain")
    e.add_argument("--seed", type=int, default=1)
    e.add_argument("--snr", type=float, nargs="+", default=[10, 12, 15, 20, 30])
    e.add_argument("--n-star", type=int, default=20)
    e.add_argument("--n-topo", type=int, default=10)
    a = ap.parse_args()
    return {"run": cmd_run, "inject": cmd_inject, "fetch": cmd_fetch, "real": cmd_real}[a.cmd](a) and 0


if __name__ == "__main__":
    sys.exit(main())
