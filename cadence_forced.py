#!/usr/bin/env python3
"""cadence_forced.py - measure a candidate in EVERY scan of a cadence, along its predicted track.

turboSETI's ON/OFF filter asks "was there a HIT in this OFF scan?". A missing
hit is not a measurement: the signal may simply have stayed under turboSETI's
threshold, been suppressed by a neighbour, or drifted a little differently
(the 1674.9977 MHz track of GJ699 had no OFF hits yet sits at SNR 31 and 33 in
both later OFFs). This tool does the forced measurement instead:

  1. take the candidate's frequency and drift where it was found, refine both
     with a fine drift-and-offset search in that scan;
  2. predict its channel track in every scan from the elapsed time
     (f + drift * dt; optional Earth-term curvature from earth_term.py);
  3. read ONLY that window of each scan (h5py hyperslab; HTTP range reads via
     fsspec when a URL is given and fsspec is installed);
  4. sum power ALONG the predicted track, with a box as wide as the tone's own
     Doppler smear (|drift| * tsamp per spectrum) plus 1.5 channels;
  5. take the best offset inside a small search window that grows only with
     the drift error x elapsed time (+ curvature bend when it is not modelled),
     and score it against the robust spread (median, 1.4826 MAD) of the SAME
     track shape at other offsets in the same window.

Per scan it reports z (significance against those null tracks) and the TRUE
radiometer SNR of the power found (the planted-tone units of inject.py:
time-summed power over the channel-to-channel noise of a single-channel sum),
plus the true SNR the bar corresponds to. "Absent" therefore means "measured,
and below this bar" - and a scan that could not be measured is None, never
absent.

Verdicts (8), in order of precedence:
               INCOMPLETE   the scan it was found in cannot be measured, or (after the
                            OFF tests below) some other scan cannot be measured
               NOT-RECOVERED  measured, but not seen even in the scan it was found in:
                            the check could not TEST it (not "explained")
               UNTESTED     no OFF scan in the cadence
               RFI          present ON THE TRACK in an OFF, at a drift consistent with
                            the track (same-emitter evidence)
               OFF-OCCUPIED an OFF has power near the track, but only through the
                            free-drift look or at an inconsistent drift: something is
                            there, maybe the same emitter, maybe a neighbour (ambiguous)
               PASSES       present in every ON, absent from every OFF, all measured
               ONE-OFF      present in exactly one ON (the one it was found in)
               PARTIAL      present in some, not all, ONs
PASSES is rung 2 at most: "in every ON, not in any OFF, not yet explained".
In crowded bands an OFF often holds an unrelated carrier near the track, so a
sky-like ON-only tone there can come back RFI or OFF-OCCUPIED (measured rates in
FORCED_CADENCE_RESULT.md); in those bands an RFI call is E1, not E2.

    python cadence_forced.py selftest
    python cadence_forced.py run --events events_f2.csv --files ON1.h5 OFF1.h5 ... --out forced_f2.csv
    python cadence_forced.py control --files data/hi_*GJ699*.npz      # plants in real data
    python cadence_forced.py control --files h5s.lst --near-events events_f2.csv --where ON
                                         # plants 3-20 kHz from real hits (crowded bands)

Idea credited to seti-hunter's weirdness_engine.py (no licence; reimplemented
from a description, no code copied).
"""
import argparse
import json
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

VERDICTS = ["RFI", "OFF-OCCUPIED", "PASSES", "ONE-OFF", "PARTIAL", "INCOMPLETE", "NOT-RECOVERED",
            "UNTESTED"]

DEFAULTS = dict(
    z_bar=6.0,          # presence: z >= z_bar against the null tracks
    base_ch=2.0,        # search half-width at zero elapsed time, channels
    drift_err_q=1.0,    # drift error floor (quanta, foff/((nt-1)tsamp)); see drift_error()
    box_pad=1.5,        # box width = smear + box_pad channels
    n_null=512,         # null offsets each side, channels
    min_null=128,       # fewer usable null offsets -> unmeasured
    scan_drift_q=1.0,   # per-scan drift wiggle +-this many quanta (3 values)
    refine=True,        # refine freq+drift in the found scan first
    refine_q=4.0, refine_ch=3.0,
    curvature=False,    # model Earth-term curvature in the prediction
    med_width=129,      # running-median bandpass removal along frequency
    nearby_hz=100.0,    # OFF scans only: also look +-this far from the track mid ...
    nearby_drift=0.75,  # ... at ANY drift within +-this (Hz/s) of the track's: RFI wanders
    same_drift_k=2.0,   # on-track OFF presence counts as RFI only if its best free drift in the
                        # on-track window is within same_drift_k * drift_error() of the track's
)


def robust(x):
    x = np.asarray(x, dtype=np.float64)
    med = float(np.median(x))
    return med, float(1.4826 * np.median(np.abs(x - med)))


# --------------------------------------------------------------------- scans
class Scan:
    """One scan of a cadence. Times in MJD, foff/fch1 in Hz (fch1 = frequency of
    data index 0), data indexed [spectrum, channel]."""
    kind = "ON"
    name = "scan"
    dc_period = None      # coarse-channel length; its centre channel is the DC spike
    dc_origin = 0

    def times_s(self):
        """Centre time of each spectrum, seconds since MJD 0."""
        return self.tstart * 86400.0 + (np.arange(self.nt) + 0.5) * self.tsamp

    def read(self, c0, c1):
        raise NotImplementedError


class ArrayScan(Scan):
    def __init__(self, data, fch1_hz, foff_hz, tstart, tsamp, kind="ON", name="array",
                 dc_period=None, dc_origin=0):
        self.data = data
        self.nt, self.nchan = data.shape
        self.fch1, self.foff, self.tstart, self.tsamp = float(fch1_hz), float(foff_hz), float(tstart), float(tsamp)
        self.kind, self.name, self.dc_period, self.dc_origin = kind, name, dc_period, dc_origin

    def read(self, c0, c1):
        return np.asarray(self.data[:, c0:c1], dtype=np.float64)


class H5Scan(Scan):
    """Filterbank HDF5 (local path or http(s) URL). Only the hyperslab is read;
    a large chunk cache keeps a coarse channel decompressed between candidates."""
    def __init__(self, path, kind=None, cache_mb=144):
        import h5py
        import hdf5plugin  # noqa: F401
        self.path = str(path)
        self.name = Path(self.path).name
        self.kind = kind or ("OFF" if "_OFF" in self.name.upper() else "ON")
        kw = dict(rdcc_nbytes=cache_mb * 2 ** 20, rdcc_nslots=100003)
        if self.path.startswith(("http://", "https://")):
            try:
                import fsspec
            except ImportError:
                raise SystemExit("HTTP range reads need fsspec (pip install fsspec aiohttp)")
            self._fh = fsspec.open(self.path, block_size=4 * 2 ** 20).open()
            self.h = h5py.File(self._fh, "r", **kw)
        else:
            self.h = h5py.File(self.path, "r", **kw)
        d = self.h["data"]
        a = dict(d.attrs)
        self.d = d
        self.nt, self.nchan = d.shape[0], d.shape[-1]
        self.fch1, self.foff = float(a["fch1"]) * 1e6, float(a["foff"]) * 1e6
        self.tstart, self.tsamp = float(a["tstart"]), float(a["tsamp"])
        self.ra_deg, self.dec_deg = float(a.get("src_raj", np.nan)) * 15.0, float(a.get("src_dej", np.nan))
        self.dc_period = d.chunks[-1] if d.chunks and d.chunks[-1] > 4096 else None

    def read(self, c0, c1):
        x = self.d[:, 0, c0:c1] if self.d.ndim == 3 else self.d[:, c0:c1]
        return np.asarray(x, dtype=np.float64)


def npz_scan(path):
    """One coarse channel saved by hi_starframe.py fetch (fch1 of the WHOLE file + ch0)."""
    z = np.load(path, allow_pickle=True)
    g = {k: (z[k].item() if z[k].ndim == 0 else z[k]) for k in z.files if k != "data"}
    foff = float(g["foff"]) * 1e6
    s = ArrayScan(z["data"].astype(np.float32), float(g["fch1"]) * 1e6 + foff * int(g["ch0"]), foff,
                  float(g["tstart"]), float(g["tsamp"]), kind=str(g["kind"]), name=Path(path).name,
                  dc_period=z["data"].shape[1], dc_origin=0)
    s.ra_deg, s.dec_deg = float(g["ra_h"]) * 15.0, float(g["dec_deg"])
    return s


def open_scan(p):
    return npz_scan(p) if str(p).endswith(".npz") else H5Scan(p)


# --------------------------------------------------------------- measurement
def quantum(scan):
    return abs(scan.foff) / ((scan.nt - 1) * scan.tsamp)


def track(scan, f_ref, t_ref, drift, curv=0.0, t_curv=None):
    """Fractional channel of the track in each spectrum of `scan`, and dt (s).
    Curvature is pivoted at t_curv (the found scan), default t_ref."""
    tk = scan.times_s()
    dt = tk - t_ref
    dc = tk - (t_ref if t_curv is None else t_curv)
    f = f_ref + drift * dt + 0.5 * curv * dc ** 2
    return (f - scan.fch1) / scan.foff, dt


def _box_sums(cum, pos, width):
    """Sum of each spectrum over [pos-w/2, pos+w/2] (fractional channels,
    exact integral of the piecewise-constant spectrum), summed over spectra.
    cum: (nt, nw+1) cumulative sums with a leading 0; pos: (nt, n_o) local positions."""
    grid = np.arange(cum.shape[1], dtype=np.float64)
    out = np.zeros(pos.shape[1])
    for k in range(cum.shape[0]):
        out += np.interp(pos[k] + 0.5 + width / 2, grid, cum[k]) - np.interp(pos[k] + 0.5 - width / 2, grid, cum[k])
    return out


def _normalise(x, scan, c0, med_width):
    """Bandpass removal ALONG FREQUENCY per spectrum (a time-median would erase
    zero-drift tones); masks the DC spike of each coarse channel."""
    from scipy.ndimage import median_filter
    bp = median_filter(x, size=(1, med_width), mode="nearest")
    good = bp > 0
    y = np.where(good, x / np.where(good, bp, 1.0) - 1.0, 0.0)
    y[~np.isfinite(y)] = 0.0
    if scan.dc_period:
        ch = np.arange(c0, c0 + x.shape[1]) - scan.dc_origin
        dc = np.abs((ch % scan.dc_period) - scan.dc_period // 2) <= 1
        y[:, dc] = 0.0
    return y, good


def measure(scan, f_ref, t_ref, drift, opts, curv_model=0.0, curv_tol=0.0, half_ch=None,
            drift_grid=None, t_curv=None):
    """Forced measurement of one track in one scan. Returns a dict; measured=False
    with a reason when the scan cannot say present or absent."""
    o = opts
    q = quantum(scan)
    p, dt = track(scan, f_ref, t_ref, drift, curv_model, t_curv)
    dtmax = float(np.max(np.abs(dt)))
    smear = abs(drift) * scan.tsamp / abs(scan.foff)
    width = smear + o["box_pad"]
    if half_ch is None:
        half_ch = o["base_ch"] + (o["drift_err_q"] * q * dtmax + 0.5 * abs(curv_tol) * dtmax ** 2) / abs(scan.foff)
    guard = width + 2.0
    wig = o["scan_drift_q"] * q
    if drift_grid is None:
        drift_grid = [drift - wig, drift, drift + wig] if wig > 0 else [drift]
    tmid = scan.times_s().mean()
    extra = max(abs(dd - drift) for dd in drift_grid) * scan.nt * scan.tsamp / abs(scan.foff)
    lo_need, hi_need = p.min() - half_ch - width - extra, p.max() + half_ch + width + extra
    res = dict(scan=scan.name, kind=scan.kind, measured=False, present=None, z=None, snr=None,
               bar_snr=None, offset_ch=None, half_ch=round(float(half_ch), 2), dt_s=round(float(dt.mean()), 1),
               reason="")
    if lo_need < 0 or hi_need > scan.nchan - 1:
        res["reason"] = "track leaves the band"
        return res
    pad = half_ch + guard + o["n_null"] + width + extra + o["med_width"] // 2 + 2
    c0 = int(max(0, math.floor(p.min() - pad)))
    c1 = int(min(scan.nchan, math.ceil(p.max() + pad) + 1))
    x = scan.read(c0, c1)
    y, good = _normalise(x, scan, c0, o["med_width"])
    pl = p - c0
    # data present under the track?
    ki = np.clip(np.round(pl).astype(int), 0, y.shape[1] - 1)
    if (~good[np.arange(scan.nt), ki]).mean() > 0.2 or np.all(x[:, max(0, int(pl.min()) - 2):int(pl.max()) + 3] == 0):
        res["reason"] = "no data under the track (blanked)"
        return res
    step = 0.5
    lim_lo = (width / 2 + extra + 1) - pl.min()
    lim_hi = (y.shape[1] - 1 - width / 2 - extra - 1) - pl.max()
    offs = np.arange(-(half_ch + guard + o["n_null"]), half_ch + guard + o["n_null"] + step, step)
    offs = offs[(offs >= lim_lo) & (offs <= lim_hi)]
    search = np.abs(offs) <= half_ch
    null = np.abs(offs) >= half_ch + guard
    if null.sum() < o["min_null"] / step or not search.any():
        res["reason"] = "too few null tracks in the window"
        return res
    best = None
    cum = np.concatenate((np.zeros((y.shape[0], 1)), np.cumsum(y, axis=1)), axis=1)
    for dd in drift_grid:
        pk = pl + (dd - drift) * (scan.times_s() - tmid) / scan.foff
        s = _box_sums(cum, pk[:, None] + offs[None, :], abs(dd) * scan.tsamp / abs(scan.foff) + o["box_pad"])
        med, sig = robust(s[null])
        if sig <= 0:
            continue
        z = (s - med) / sig
        j = int(np.flatnonzero(search)[np.argmax(z[search])])
        if best is None or z[j] > best[0]:
            best = (float(z[j]), j, dd, s[j] - med, sig, pk)
    if best is None:
        res["reason"] = "zero noise (constant data)"
        return res
    zb, j, dd, power, sig, pk = best
    # radiometer unit: integer-aligned single-channel sums along the same track at null offsets
    ints = np.round(offs[null][::2]).astype(int)
    rows = np.round(pk).astype(int)[:, None] + ints[None, :]
    okc = (rows >= 0).all(0) & (rows < y.shape[1]).all(0)
    s1 = y[np.arange(scan.nt)[:, None], rows[:, okc]].sum(0)
    _, sig1 = robust(s1)
    sig1 = sig1 if sig1 > 0 else sig
    res.update(measured=True, present=bool(zb >= o["z_bar"]), z=round(zb, 2),
               snr=round(float(power / sig1), 1), bar_snr=round(float(o["z_bar"] * sig / sig1), 1),
               offset_ch=float(offs[j]), drift_used=round(float(dd), 5), reason="")
    return res


def refine(scan, f_ref, t_ref, drift, opts):
    """Best freq+drift near the reported ones in the found scan (fine grid)."""
    q = quantum(scan)
    grid = drift + q * np.arange(-opts["refine_q"], opts["refine_q"] + 1e-9, 0.25)
    r = measure(scan, f_ref, t_ref, drift, opts, half_ch=opts["refine_ch"], drift_grid=list(grid))
    if not r["measured"]:
        return f_ref, drift, r
    d2 = r["drift_used"]
    # the offset is at the scan mid (pivot) = t_ref, so it moves f_ref directly
    return f_ref + r["offset_ch"] * scan.foff, d2, r


def verdict(kinds, present, found=None, how=None):
    """kinds: 'ON'/'OFF' per scan; present: True / False / None (unmeasured);
    how (optional, per scan): 'track' (on-track, consistent drift), 'crossing'
    (on-track box, inconsistent drift) or 'nearby' (free-drift look only).
    Without `how`, every OFF presence counts as 'track'."""
    if found is not None and present[found] is None:
        return "INCOMPLETE"
    if found is not None and present[found] is not True:
        return "NOT-RECOVERED"
    offs = [i for i, k in enumerate(kinds) if k == "OFF"]
    ons = [i for i, k in enumerate(kinds) if k == "ON"]
    if not offs:
        return "UNTESTED"
    how = how or ["track"] * len(kinds)
    if any(present[i] is True and how[i] == "track" for i in offs):
        return "RFI"
    if any(present[i] is True for i in offs):
        return "OFF-OCCUPIED"
    if any(p is None for p in present):
        return "INCOMPLETE"
    n_on = sum(present[i] is True for i in ons)
    if n_on == len(ons):
        return "PASSES"
    if n_on == 1:
        return "ONE-OFF"
    return "PARTIAL"


def drift_error(scan, drift, z, opts):
    """Drift error (Hz/s) of a single-scan estimate. Measured on synthetic
    smeared plants after refine() (see selftest/README): ~0.5 quantum for slow
    tones at z ~ 20-100, up to ~5 quanta at 3 Hz/s, where the smear box is
    nearly blind to small slope changes. Generous by design."""
    smear = abs(drift) * scan.tsamp / abs(scan.foff)
    return quantum(scan) * (opts["drift_err_q"] + 0.15 * smear + 20.0 / max(z or 1.0, 1.0))


def forced(cand, scans, opts, curv_per_hz=0.0):
    """cand: dict(f_hz (at the centre of spectrum 0 of the found scan, as
    turboSETI reports it), drift, found=index into scans). Returns a dict.

    Scans are measured nearest-first. Every scan where the track is PRESENT
    becomes an anchor; with two or more anchors the drift is re-fitted across
    them (phase connection) and its error shrinks with the lever arm, so later
    windows stay small. ON and OFF anchors are treated the same."""
    i0 = cand["found"]
    fs = scans[i0]
    t_ref = fs.times_s().mean()
    f_ref = cand["f_hz"] + cand["drift"] * (t_ref - fs.times_s()[0])
    drift = cand["drift"]
    curv = curv_per_hz * f_ref
    cm = curv if opts["curvature"] else 0.0
    ct = 0.0 if opts["curvature"] else curv
    ref = None
    if opts["refine"]:
        f_ref, drift, ref = refine(fs, f_ref, t_ref, drift, opts)
    per = [None] * len(scans)
    per[i0] = measure(fs, f_ref, t_ref, drift, opts, curv_model=cm, t_curv=t_ref)
    z0 = per[i0]["z"]
    dd1 = drift_error(fs, drift, z0, opts)
    T = fs.nt * fs.tsamp
    anchors = [(t_ref, f_ref)]
    if per[i0]["present"]:
        anchors = [(t_ref, f_ref + per[i0]["offset_ch"] * fs.foff)]
    order = sorted((i for i in range(len(scans)) if i != i0),
                   key=lambda i: abs(scans[i].times_s().mean() - t_ref))
    for i in order:
        s = scans[i]
        ta = np.array([a[0] for a in anchors])
        fa = np.array([a[1] for a in anchors]) - 0.5 * cm * (ta - t_ref) ** 2
        if len(anchors) >= 2:
            d_use, a0 = np.polyfit(ta - ta.mean(), fa, 1)
            t_lin, f_lin = ta.mean(), a0
            dd = dd1 * T / (ta.max() - ta.min())
        else:
            d_use, t_lin, f_lin, dd = drift, ta[0], fa[0], dd1
        tk = s.times_s()
        dt_near = float(np.max(np.min(np.abs(tk[:, None] - ta[None, :]), axis=1)))
        dt_c = float(np.max(np.abs(tk - t_ref)))
        half = opts["base_ch"] + (dd * dt_near + 0.5 * abs(ct) * dt_c ** 2) / abs(s.foff)
        m = measure(s, f_lin, t_lin, d_use, opts, curv_model=cm, t_curv=t_ref, half_ch=half)
        m["drift_model"] = round(float(d_use), 5)
        on_track = bool(m["present"])
        m["how"] = "track" if on_track else ""
        if s.kind == "OFF" and on_track:
            # Same emitter, or an unrelated carrier crossing the box? Free drift, SAME window:
            # the presence counts as RFI only if the best drift there agrees with the track.
            q = quantum(s)
            nd = opts["nearby_drift"]
            grid = list(d_use + q * np.arange(-round(nd / q), round(nd / q) + 1))
            m3 = measure(s, f_lin, t_lin, d_use, opts, curv_model=cm, t_curv=t_ref, half_ch=half,
                         drift_grid=grid)
            tol = opts["same_drift_k"] * drift_error(s, d_use, m["z"], opts)
            if m3["measured"]:
                m["free_drift"] = m3["drift_used"]
                if abs(m3["drift_used"] - d_use) > tol:
                    m["how"] = "crossing"
                    m["reason"] = (f"on-track box, but best drift there {m3['drift_used']:+.3f} Hz/s "
                                   f"(track {d_use:+.3f} +- {tol:.3f}): maybe a neighbour")
        if s.kind == "OFF" and opts["nearby_hz"] > 0 and m["measured"]:
            # Local interference does not keep one drift (GJ699 1674.9977 MHz: -0.33, -0.53,
            # -0.20, +0.21 Hz/s in successive scans). Absence from an OFF must survive a
            # free-drift look near the track, or a wandering transmitter could PASS.
            q = quantum(s)
            nd = opts["nearby_drift"]
            grid = list(d_use + q * np.arange(-round(nd / q), round(nd / q) + 1))
            m2 = measure(s, f_lin, t_lin, d_use, opts, curv_model=cm, t_curv=t_ref,
                         half_ch=max(half, opts["nearby_hz"] / abs(s.foff)), drift_grid=grid)
            m["nearby_z"] = m2["z"]
            if m2["measured"]:
                m["nearby_drift"] = m2["drift_used"]
                m["nearby_offset_hz"] = round(m2["offset_ch"] * abs(s.foff), 1)
                if m2["present"] and not m["present"]:
                    m["present"] = True
                    m["how"] = "nearby"
                    m["reason"] = (f"nearby: z {m2['z']} at drift {m2['drift_used']:+.3f} Hz/s, "
                                   f"{m2['offset_ch'] * abs(s.foff):+.0f} Hz off the track")
        per[i] = m
        if on_track and m["how"] == "track":   # only consistent ON-TRACK presence phase-connects
            tm = tk.mean()
            f_at = f_lin + d_use * (tm - t_lin) + 0.5 * cm * (tm - t_ref) ** 2 + m["offset_ch"] * s.foff
            anchors.append((tm, f_at))
    v = verdict([s.kind for s in scans], [r["present"] for r in per], i0,
                [r.get("how", "track" if r["present"] else "") for r in per])
    return dict(verdict=v, f_ref_hz=f_ref, drift=drift, scans=per,
                refine_z=None if ref is None else ref["z"], drift_err=dd1)


def curvature_per_hz(scans):
    """Earth-term curvature (Hz/s^2 per Hz of sky frequency) at the cadence centre."""
    s = scans[0]
    if not np.isfinite(getattr(s, "ra_deg", np.nan)):
        return 0.0
    from earth_term import earth_curvature
    mid = 0.5 * (scans[0].tstart + scans[-1].tstart + scans[-1].nt * scans[-1].tsamp / 86400)
    return float(earth_curvature(1.0, s.ra_deg, s.dec_deg, mid))


# ------------------------------------------------------------------ controls
def negative_controls(scans, opts, n, rng, chan_pool, drift_pool, curv_per_hz=0.0, avoid=None):
    """Random positions (same coarse channels, same drift mix, away from known
    hits) through the same machinery. Returns rates and verdict counts."""
    ons = [i for i, s in enumerate(scans) if s.kind == "ON"]
    pres = {"ON": [0, 0], "OFF": [0, 0]}
    found_pres = [0, 0]
    counts = {}
    avoid = np.sort(np.asarray(avoid if avoid is not None else [], dtype=float))
    tries = 0
    done = 0
    while done < n and tries < 20 * n:
        tries += 1
        i = int(rng.choice(ons))
        s = scans[i]
        jit = min(30000.0, s.nchan / 2)
        c = float(rng.choice(chan_pool)) + rng.uniform(-jit, jit)
        f = s.fch1 + s.foff * c
        if avoid.size:
            k = np.searchsorted(avoid, f)
            near = [avoid[m] for m in (k - 1, k) if 0 <= m < avoid.size]
            if near and min(abs(f - a) for a in near) < 2000.0:
                continue
        if not (2000 < c < s.nchan - 2000):
            continue
        r = forced(dict(f_hz=f, drift=float(rng.choice(drift_pool)), found=i), scans, opts, curv_per_hz)
        done += 1
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
        for j, m in enumerate(r["scans"]):
            if not m["measured"]:
                continue
            if j == i:
                found_pres[0] += m["present"]
                found_pres[1] += 1
            else:
                pres[m["kind"]][0] += m["present"]
                pres[m["kind"]][1] += 1
    rate = lambda a: round(a[0] / a[1], 4) if a[1] else None  # noqa: E731
    return dict(n=done, false_presence_found_scan=rate(found_pres),
                false_presence_other_on=rate(pres["ON"]), false_presence_off=rate(pres["OFF"]),
                verdicts=counts)


# ----------------------------------------------------------------------- run
def load_events(path):
    import pandas as pd
    d = pd.read_csv(path)
    return d


def cmd_run(a):
    import pandas as pd
    opts = dict(DEFAULTS, z_bar=a.zbar, curvature=a.curvature, refine=not a.no_refine,
                drift_err_q=a.drift_err_q, nearby_hz=a.nearby_hz)
    files = a.files
    if len(files) == 1 and files[0].endswith(".lst"):
        files = [ln.strip() for ln in open(files[0]) if ln.strip()]
    scans = [open_scan(p) for p in files]
    print("cadence:", " ".join(f"{s.kind}" for s in scans))
    names = [s.name for s in scans]
    curv = curvature_per_hz(scans)
    f_mid = 1.5e9
    print(f"Earth-term curvature {curv*f_mid:+.2e} Hz/s^2 at 1.5 GHz; bend over the cadence "
          f"{0.5*abs(curv*f_mid)*((scans[-1].tstart-scans[0].tstart)*86400+288)**2/abs(scans[0].foff):.2f} ch "
          f"({'modelled' if opts['curvature'] else 'absorbed by the search window'})")
    d = load_events(a.events)
    if a.limit:
        d = d.iloc[:a.limit]
    d = d.assign(_row=np.arange(len(d)))
    gch = (d["Freq"] * 1e6 - scans[0].fch1) / scans[0].foff
    d = d.assign(_ch=gch).sort_values("_ch")
    rows = []
    skipped = []
    t0 = time.time()
    for n, (_, r) in enumerate(d.iterrows()):
        fid = Path(str(r["FileID"])).name
        found = names.index(fid) if fid in names else None
        if found is None:
            # never drop a row silently: it is written as INCOMPLETE and counted
            skipped.append(fid)
            rows.append(dict(row=int(r["_row"]), hit_id=r.get("Hit_ID"), freq_mhz=float(r["Freq"]),
                             drift=float(r["DriftRate"]), turbo_snr=float(r["SNR"]), found="",
                             verdict="INCOMPLETE", note=f"scan file not given: {fid}"))
            continue
        res = forced(dict(f_hz=float(r["Freq"]) * 1e6, drift=float(r["DriftRate"]), found=found),
                     scans, opts, curv)
        row = dict(row=int(r["_row"]), hit_id=r.get("Hit_ID"), freq_mhz=float(r["Freq"]),
                   drift=float(r["DriftRate"]), turbo_snr=float(r["SNR"]), found=found,
                   verdict=res["verdict"], drift_refined=round(res["drift"], 5))
        for j, m in enumerate(res["scans"]):
            tag = f"{j}{m['kind']}"
            row[f"z_{tag}"] = m["z"]
            row[f"snr_{tag}"] = m["snr"]
            row[f"bar_{tag}"] = m["bar_snr"]
            row[f"present_{tag}"] = {True: 1, False: 0, None: ""}[m["present"]]
            row[f"why_{tag}"] = m["reason"]
            if m["kind"] == "OFF":
                row[f"how_{tag}"] = m.get("how", "")
                row[f"nearz_{tag}"] = m.get("nearby_z")
                row[f"neard_{tag}"] = m.get("nearby_drift")
                row[f"nearoff_{tag}"] = m.get("nearby_offset_hz")
        row["off_present"] = " ".join(f"{j}{m['kind']}" for j, m in enumerate(res["scans"])
                                      if m["kind"] == "OFF" and m["present"])
        rows.append(row)
        if (n + 1) % 100 == 0:
            print(f"  {n+1}/{len(d)}  {time.time()-t0:.0f} s", flush=True)
    out = pd.DataFrame(rows).sort_values("row")
    out.to_csv(a.out, index=False)
    counts = out["verdict"].value_counts().to_dict()
    print(f"\n{len(out)} candidates in {time.time()-t0:.0f} s -> {a.out}")
    if skipped:
        print(f"  WARNING: {len(skipped)} rows name a scan not in --files (written as INCOMPLETE): "
              f"{sorted(set(skipped))[:5]}")
    for v in VERDICTS:
        print(f"  {v:14s} {counts.get(v, 0)}")
    rng = np.random.default_rng(a.seed)
    neg = None
    if a.nneg:
        t1 = time.time()
        avoid = list(d["Freq"].values * 1e6)
        for fn in (a.dats or []):
            avoid += [float(ln.split()[3]) * 1e6 for ln in open(fn)
                      if ln.strip() and not ln.startswith("#")]
        neg = negative_controls(scans, opts, a.nneg, rng, gch.values, d["DriftRate"].values, curv,
                                avoid=avoid)
        print(f"\nnegative controls ({neg['n']} random positions, {time.time()-t1:.0f} s): "
              f"false presence found-scan {neg['false_presence_found_scan']}, other ON "
              f"{neg['false_presence_other_on']}, OFF {neg['false_presence_off']}; verdicts {neg['verdicts']}")
    if a.json:
        Path(a.json).write_text(json.dumps(dict(events=a.events, files=names, opts=opts,
                                                counts=counts, negative_controls=neg,
                                                rows_scan_not_given=len(skipped),
                                                curvature_hz_s2_at_1p5GHz=curv * f_mid), indent=1))
    passes = out[out["verdict"] == "PASSES"]
    if len(passes):
        print("\nPASSES rows (not yet explained - try to kill them):")
        for _, r in passes.iterrows():
            print(f"  {r['freq_mhz']:.6f} MHz drift {r['drift']:+.4f} found in scan {r['found']}")
    return 2 if skipped else 0


# ------------------------------------------------------------------- control
def box_plant(data, foff_hz, tsamp, plants):
    """setigen-free injector (fallback, and a cross-check of inject.plant): a
    1-channel box tone whose centre moves ch + drift*t/foff, integrated exactly
    over each spectrum (256 sub-steps, each split linearly between its two
    nearest channels). Same units as inject.plant (inject.level_for), same
    convention: ch = channel at the START of spectrum 0, data index order."""
    import inject
    nt, nc = data.shape
    for p in plants:
        lvl = p.get("level") or inject.level_for(data, p["snr"])
        u = (np.arange(256) + 0.5) / 256
        for k in range(nt):
            c = p["ch"] + p["drift"] * (k + u) * tsamp / foff_hz
            j = np.floor(c).astype(int)
            w = c - j
            for jj, ww in ((j, 1 - w), (j + 1, w)):
                ok = (jj >= 0) & (jj < nc)
                np.add.at(data[k], jj[ok], (lvl * ww[ok] / u.size).astype(data.dtype))
    return data


def get_injector(name="auto"):
    """'setigen' (inject.plant), 'box' (box_plant) or 'auto' (setigen if importable)."""
    if name in ("auto", "setigen"):
        try:
            import setigen  # noqa: F401
            import inject
            return "setigen", inject.plant
        except ImportError:
            if name == "setigen":
                raise SystemExit("setigen is not installed; use --injector box")
    return "box", lambda data, foff, ts, pl: box_plant(data, foff, ts, pl)


def _scan_key(p):
    """Sort scans by (MJD, seconds) in a GUPPI-style name; fall back to the name."""
    import re
    m = re.search(r"_(\d{5})_(\d{5})_", Path(str(p)).name)
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0, Path(str(p)).name)
def cmd_control(a):
    """Positive + negative controls on REAL data: plant a smeared tone along one
    phase-connected track into the ONs (expect PASSES) or into every scan
    (expect RFI), at several strengths and drifts."""
    import pandas as pd
    opts = dict(DEFAULTS, z_bar=a.zbar, curvature=a.curvature, nearby_hz=a.nearby_hz)
    files = a.files
    if len(files) == 1 and files[0].endswith(".lst"):
        files = [ln.strip() for ln in open(files[0]) if ln.strip()]
    inj_name, plant = get_injector(a.injector)
    print(f"injector: {inj_name}")
    raw = [open_scan(p) for p in sorted(files, key=_scan_key)]
    scans0 = raw
    near = None
    if a.near_events:
        near = pd.read_csv(a.near_events)["Freq"].values * 1e6
        print(f"plants 3-20 kHz from a random hit of {a.near_events} ({len(near)} rows)")
    print("cadence:", " ".join(s.kind for s in scans0), f"({scans0[0].nchan} channels each)")
    curv = curvature_per_hz(scans0)
    rng = np.random.default_rng(a.seed)
    nc = scans0[0].nchan
    centre = nc // 2
    results = []
    cases = [(d, snr, where) for d in a.drifts for snr in a.snrs for where in a.where]
    dcp = scans0[0].dc_period
    for d, snr, where in cases:
        for rep in range(a.reps):
            while True:
                if near is not None:
                    # crowded control: a few kHz from a real hit, in the same coarse channel
                    f = float(rng.choice(near)) + rng.choice([-1, 1]) * rng.uniform(3e3, 20e3)
                    c = (f - scans0[0].fch1) / scans0[0].foff
                    dc_far = dcp is None or abs((c % dcp) - dcp // 2) > 200
                    if 14000 < c < nc - 14000 and dc_far:
                        break
                    continue
                # start channel well inside the slice, away from DC and the HI hump
                c = float(rng.uniform(40000, nc - 40000))
                f = scans0[0].fch1 + scans0[0].foff * c
                if abs(c - centre) > 20000 and not (1420.0e6 < f < 1420.8e6):
                    break
            lo, hi = int(c) - 12000, int(c) + 12000
            sc = []
            fs = scans0[0]
            t_ref = fs.times_s()[0]
            for s in scans0:
                x = s.read(lo, hi).astype(np.float32)
                ns = ArrayScan(x, s.fch1 + s.foff * lo, s.foff, s.tstart, s.tsamp, s.kind, s.name)
                if where == "ALL" or s.kind == "ON":
                    t0 = s.tstart * 86400.0
                    dt0 = t0 - t_ref
                    f0 = f + d * dt0 + 0.5 * curv * f * dt0 ** 2
                    plant(ns.data, s.foff, s.tsamp, [dict(ch=(f0 - ns.fch1) / s.foff, drift=d, snr=snr)])
                sc.append(ns)
            # the 'hit' as a search would report it: freq at spectrum-0 centre (= f, since t_ref is
            # that centre), drift off by 0.4 quantum
            f_hit = f
            r = forced(dict(f_hz=f_hit, drift=d + 0.4 * quantum(fs), found=0), sc, opts, curv)
            want = "PASSES" if where == "ON" else "RFI"
            results.append(dict(freq_mhz=round(f / 1e6, 6), drift=d, snr=snr, where=where,
                                verdict=r["verdict"], ok=r["verdict"] == want,
                                snr_read=[m["snr"] for m in r["scans"]], z=[m["z"] for m in r["scans"]]))
            print(f"  plant drift {d:+.2f} true SNR {snr:4.0f} in {where:3s} -> {r['verdict']:13s} "
                  f"{'ok' if r['verdict']==want else 'MISS'}  read SNR "
                  + " ".join(f"{m['snr'] if m['snr'] is not None else '-'}" for m in r["scans"]))
    # negative controls on the unplanted slice
    neg = negative_controls(scans0, opts, a.nneg, rng, [nc / 2 + 200000, nc / 2 - 200000],
                            [0.0105, -0.12, -0.3, 1.0, -2.0], curv) if a.nneg else None
    if neg:
        print(f"negative controls ({neg['n']}): false presence found-scan {neg['false_presence_found_scan']}, "
              f"other ON {neg['false_presence_other_on']}, OFF {neg['false_presence_off']}; {neg['verdicts']}")
    summ = {}
    for r in results:
        k = f"{r['where']} drift {r['drift']:+.2f} snr {r['snr']:.0f}"
        summ.setdefault(k, [0, 0])
        summ[k][0] += r["ok"]
        summ[k][1] += 1
    print("\nsummary (correct / planted):")
    for k, v in summ.items():
        print(f"  {k:30s} {v[0]}/{v[1]}")
    for where in a.where:
        vc = {}
        for r in results:
            if r["where"] == where:
                vc[r["verdict"]] = vc.get(r["verdict"], 0) + 1
        print(f"  verdicts of plants in {where}: {vc}")
    if a.json:
        Path(a.json).write_text(json.dumps(dict(results=results, negative_controls=neg, opts=opts,
                                                injector=inj_name, near_events=a.near_events), indent=1))
    return 0


# ------------------------------------------------------------------ selftest
def _synthetic(rng, nchan=8192, kinds=("ON", "OFF", "ON", "OFF", "ON", "OFF"), ripple=0.4):
    nt, foff, ts = 16, -2.835503418452676, 17.986224128
    starts = [48006, 48334, 48661, 48988, 49314, 49640][:len(kinds)]
    bp = 1.0 + ripple * np.cos(np.linspace(0, 3 * np.pi, nchan)) ** 2     # bandpass-like ripple
    out = []
    for k, s in zip(kinds, starts):
        x = (rng.chisquare(64, (nt, nchan)) / 64 * bp).astype(np.float32)
        out.append(ArrayScan(x, 1.5e9, foff, 57457.0 + s / 86400, ts, kind=k, name=f"{k}{s}"))
    return out


def _plant_track(scans, f, d, snr, which, off_drift=None, plant=None):
    """Plant along one straight track; with off_drift, the OFF scans get a tone
    at the track's mid-scan position but with a DIFFERENT drift (wandering RFI)."""
    if plant is None:
        import inject
        plant = inject.plant
    t_ref = scans[0].times_s()[0] - 0.5 * scans[0].tsamp
    for i, s in enumerate(scans):
        if i not in which:
            continue
        f0 = f + d * (s.tstart * 86400 - t_ref)
        dd = d
        if off_drift is not None and s.kind == "OFF":
            half = 0.5 * s.nt * s.tsamp
            f0, dd = f0 + (d - off_drift) * half, off_drift
        ch = (f0 - s.fch1) / s.foff
        if -50 < ch < s.nchan + 50:
            plant(s.data, s.foff, s.tsamp, [dict(ch=ch, drift=dd, snr=snr)])


def selftest():
    rng = np.random.default_rng(7)
    checks = []
    opts = dict(DEFAULTS)
    on3 = [0, 2, 4]

    def case(name, which, d, snr, want, ch=3000.0, kinds=None, ripple=0.4, off_drift=None, o=None,
             plant=None, extra=None):
        sc = _synthetic(rng, kinds=kinds, ripple=ripple) if kinds else _synthetic(rng, ripple=ripple)
        f = sc[0].fch1 + sc[0].foff * ch
        _plant_track(sc, f, d, snr, which, off_drift, plant)
        if extra:
            extra(sc, f)
        # reported hit: centre of spectrum 0, drift off by 0.4 quantum, freq off by 0.3 channel
        f_hit = f + d * 0.5 * sc[0].tsamp + 0.3 * sc[0].foff
        r = forced(dict(f_hz=f_hit, drift=d + 0.4 * quantum(sc[0]), found=0), sc, o or opts)
        ok = r["verdict"] == want
        snrs = " ".join(f"{m['snr'] if m['snr'] is not None else '-'}" for m in r["scans"])
        checks.append((f"{name}: {r['verdict']} (want {want}); per-scan true SNR {snrs}", ok))
        return r

    case("sky-like tone in ONs only", on3, -0.12, 40, "PASSES")
    r = case("sky-like tone in ONs only, flat bandpass", on3, -0.12, 40, "PASSES", ripple=0.0)
    rd = [m["snr"] for m in r["scans"][::2]]
    checks.append((f"true SNR read back {rd} for planted 40 (each within 25 %)",
                   all(x is not None and abs(x / 40 - 1) < 0.25 for x in rd)))
    case("fast sky-like tone in ONs only (smeared)", on3, -1.5, 80, "PASSES")
    case("tone in every scan", range(6), -0.30, 40, "RFI")
    case("tone in every scan, zero drift", range(6), 0.0, 40, "RFI")
    case("wandering RFI: ONs at -0.30, OFFs at +0.20 Hz/s (seen only by the free-drift look)",
         range(6), -0.30, 40, "OFF-OCCUPIED", off_drift=0.20)
    case("same, free-drift OFF look switched off (shows why it exists)", range(6), -0.30, 40, "PASSES",
         off_drift=0.20, o=dict(opts, nearby_hz=0.0))
    case("sky-like tone in ONs only, planted with the setigen-free box injector", on3, -0.12, 40, "PASSES",
         plant=box_plant)

    def crossing(sc, f):
        # an unrelated strong carrier in every OFF, drift +0.40 Hz/s, crossing the track mid-scan
        for s in sc:
            if s.kind != "OFF":
                continue
            t_mid = s.tstart * 86400 + 0.5 * s.nt * s.tsamp
            f_track = f + (-0.12) * (t_mid - (sc[0].tstart * 86400))
            f0 = f_track - 0.40 * 0.5 * s.nt * s.tsamp
            box_plant(s.data, s.foff, s.tsamp, [dict(ch=(f0 - s.fch1) / s.foff, drift=0.40, snr=300)])
    case("ON-only sky tone + unrelated carrier crossing the track in every OFF (not same-emitter)",
         on3, -0.12, 40, "OFF-OCCUPIED", extra=crossing)
    case("tone in first ON only", [0], -0.12, 40, "ONE-OFF")
    case("tone in ON1+ON2 only", [0, 2], -0.12, 40, "PARTIAL")
    case("track leaves the band", on3, -2.0, 60, "INCOMPLETE", ch=7300.0)
    case("nothing planted at the hit", [], -0.12, 40, "NOT-RECOVERED")
    case("no OFF scans", [0, 1, 2], -0.12, 40, "UNTESTED", kinds=("ON", "ON", "ON"))
    # verdict logic edge cases
    v1 = verdict(["ON", "OFF", "ON", "OFF"], [True, None, True, None], 0)
    checks.append((f"ON/OFF/ON/OFF = T/None/T/None -> {v1} (must not PASS)", v1 != "PASSES" and v1 == "INCOMPLETE"))
    v2 = verdict(["ON", "OFF", "ON", "OFF"], [True, False, None, False], 0)
    checks.append((f"T/F/None/F -> {v2} (must not be ONE-OFF)", v2 != "ONE-OFF" and v2 == "INCOMPLETE"))
    v3 = verdict(["ON", "OFF", "ON", "OFF"], [True, True, None, None], 0)
    checks.append((f"present in an OFF beats unmeasured scans -> {v3}", v3 == "RFI"))
    v4 = verdict(["ON", "OFF", "ON", "OFF"], [None, True, True, False], 0)
    checks.append((f"found scan unmeasurable -> {v4} (INCOMPLETE, not NOT-RECOVERED/RFI)", v4 == "INCOMPLETE"))
    v5 = verdict(["ON", "OFF", "ON", "OFF"], [True, True, True, True], 0, ["track", "crossing", "", "nearby"])
    checks.append((f"OFF presence only crossing/nearby -> {v5} (not RFI)", v5 == "OFF-OCCUPIED"))
    # cadence_verdict --forced: an untestable verdict must never explain an event
    try:
        import pandas as pd
        import cadence_verdict as CV
        base = {"z_0ON": 3.0, "snr_0ON": 4.0, "present_0ON": 0, "z_1OFF": 2.0, "snr_1OFF": 3.0,
                "present_1OFF": 0, "how_1OFF": ""}
        for v, want in (("NOT-RECOVERED", 0), ("INCOMPLETE", 0), ("UNTESTED", 0), ("RFI", 1),
                        ("OFF-OCCUPIED", 1)):
            rs, nt_ = CV.forced_reasons(pd.DataFrame([dict(base, verdict=v)]))
            checks.append((f"cadence_verdict: {v} gives {len(rs)} reason(s), {len(nt_)} note(s) "
                           f"(want {want} reason)", len(rs) == want))
    except ImportError as e:
        checks.append((f"cadence_verdict import failed: {e}", False))
    # the box injector against inject.plant (setigen): same power, same place
    try:
        import inject
        a0 = np.zeros((16, 400), np.float32) + 1.0
        a1, a2 = a0.copy(), a0.copy()
        pl = [dict(ch=200.3, drift=-0.7, snr=50, level=10.0)]
        inject.plant(a1, -2.835503418452676, 17.986224128, [dict(pl[0])])
        box_plant(a2, -2.835503418452676, 17.986224128, [dict(pl[0])])
        d1, d2 = a1 - 1.0, a2 - 1.0
        cen = lambda d: (d * np.arange(d.shape[1])).sum(1) / d.sum(1)  # noqa: E731
        dp = float(np.max(np.abs(d1.sum(1) / d2.sum(1) - 1)))
        dc = float(np.max(np.abs(cen(d1) - cen(d2))))
        # setigen samples the 1-channel box at channel centres, so its centroid snaps by up to
        # 0.5 ch; the box injector integrates exactly. Same power, same place to < 0.5 ch.
        checks.append((f"box injector vs inject.plant: power per spectrum within {dp:.3f}, "
                       f"centroid within {dc:.3f} ch (want <= 2 %, < 0.5 ch)", dp <= 0.02 and dc < 0.5))
    except ImportError:
        print("  (setigen not installed: box injector cross-check skipped)")
    # negative controls: random positions in plain noise must be absent
    sc = _synthetic(rng)
    neg = negative_controls(sc, opts, 30, rng, [4096.0], [0.0, -0.12, 0.5, -2.0])
    fp = max(neg["false_presence_found_scan"], neg["false_presence_other_on"], neg["false_presence_off"])
    checks.append((f"negative controls: false presence {fp} over {neg['n']} random tracks (<= 2 %); "
                   f"verdicts {neg['verdicts']}", fp <= 0.02 and neg["verdicts"].get("NOT-RECOVERED", 0) >= 29))
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd")
    sub.add_parser("selftest")
    p = sub.add_parser("run", help="forced check of every row of a turboSETI events CSV")
    p.add_argument("--events", required=True)
    p.add_argument("--files", nargs="+", required=True, help="scans in cadence order (or one .lst)")
    p.add_argument("--out", default="forced.csv")
    p.add_argument("--json")
    p.add_argument("--nneg", type=int, default=300, help="random negative-control positions")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--dats", nargs="*", help="turboSETI .dat files: keep negative controls 2 kHz from every hit")
    p.add_argument("--limit", type=int)
    p.add_argument("--zbar", type=float, default=DEFAULTS["z_bar"])
    p.add_argument("--drift-err-q", type=float, default=DEFAULTS["drift_err_q"])
    p.add_argument("--curvature", action="store_true", help="model Earth-term curvature")
    p.add_argument("--nearby-hz", type=float, default=DEFAULTS["nearby_hz"],
                   help="OFF free-drift look half-width, Hz (0 = off)")
    p.add_argument("--no-refine", action="store_true")
    c = sub.add_parser("control", help="plant smeared tones in real scans: ONs -> PASSES, all -> RFI")
    c.add_argument("--files", nargs="+", required=True)
    c.add_argument("--drifts", type=float, nargs="+", default=[-0.12, -0.5, 1.0])
    c.add_argument("--snrs", type=float, nargs="+", default=[40, 80])
    c.add_argument("--reps", type=int, default=3)
    c.add_argument("--nneg", type=int, default=100)
    c.add_argument("--seed", type=int, default=1)
    c.add_argument("--zbar", type=float, default=DEFAULTS["z_bar"])
    c.add_argument("--curvature", action="store_true")
    c.add_argument("--json")
    c.add_argument("--where", nargs="+", default=["ON", "ALL"], choices=["ON", "ALL"],
                   help="ON = plant in the ONs only (want PASSES); ALL = every scan (want RFI)")
    c.add_argument("--near-events", help="events CSV: plant 3-20 kHz from its hits (crowded bands)")
    c.add_argument("--injector", default="auto", choices=["auto", "setigen", "box"],
                   help="auto = inject.plant (setigen) if installed, else the setigen-free box tone")
    c.add_argument("--nearby-hz", type=float, default=DEFAULTS["nearby_hz"])
    a = ap.parse_args()
    if a.cmd in (None, "selftest"):
        return selftest()
    return cmd_run(a) if a.cmd == "run" else cmd_control(a)


if __name__ == "__main__":
    sys.exit(main())
