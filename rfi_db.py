#!/usr/bin/env python3
"""rfi_db.py - a cross-star RFI database: has anyone ELSE's beam seen this hit?

A signal from a star should appear only when the telescope points at that star.
If the same frequency, with the same drift, shows up in the hit lists of OTHER
stars - on the same night or on another night - it is almost certainly made on
Earth. This tool keeps those hit lists in one table and answers, for any hit:

    how many OTHER targets, and in how many sessions, show a matching hit?
    is that more than chance, given how crowded that part of the band is?

Matching is keyed on FREQUENCY, DRIFT and TIME, not frequency alone, so a broad
persistent interference band does not mask everything near it, and a drifting
signal that merely CROSSES a steady RFI line is not killed for crossing it.

The rule (documented, measured in RFI_DB_RESULT.md)
---------------------------------------------------
A query hit q (freq f_q at its own t=0, drift d_q, observed for T_q seconds with
channel width r_q) matches a database hit h when BOTH hold:

  1. the frequency TRACKS overlap: the interval swept by q, [f_q, f_q + d_q*T_q],
     and the one swept by h, each widened by
         slack = K_RES * max(r_q, r_h) + f * V_SLACK / c
     (K_RES = 2 channels; V_SLACK = 10 m/s by default = 3.3e-8 * f: 33 Hz at
     1 GHz, 47 Hz at 1.42 GHz. All hit tables here are TOPOCENTRIC, uncorrected
     frequencies, so Earth-made RFI keeps its frequency from one pointing to the
     next; the slack covers oscillator wander and a geostationary satellite's
     station-keeping Doppler, not a fast-moving satellite's kHz Doppler.)
  2. the drifts agree: |d_q - d_h| <= DRIFT_ABS + DRIFT_REL*max(|d_q|,|d_h|)
     + r_q/T_q + r_h/T_h       (0.02 Hz/s + 20 % + one drift quantum each)

Hits whose drift was never resolved (mid-resolution products, drift quantum
~9.8 Hz/s, KNOWN_ISSUES #7) carry drift_valid=0: they are counted on a separate
frequency-only line and NEVER contribute to a veto.

Only OTHER targets count: hits whose target has the same base name as the
query (GJ699, GJ699_OFF, ProxCen_S, ...) are left out, because a signal seen
again from the SAME star on another date is rung 4 evidence FOR it, not against
it, and because a star's own OFF scans are part of its own cadence test.

    n_t     = number of distinct other targets with a matching hit
    n_s     = number of distinct sessions (telescope + MJD day) those come from
    lambda  = matches expected by CHANCE: drift-compatible other-target hits in a
              +-200 kHz neighbourhood (M of them) x (match width / neighbourhood
              width), with M counted as at least 1: a rate measured from zero
              background hits is small, not zero
    p       = Poisson P(N >= n_t | lambda)   (counting targets against a hit-rate
              lambda is conservative: lambda >= the chance rate of targets)
    score   = -log10(p), capped at 10;  VETO when score >= 3 and n_t >= 1, or
              n_t >= 2 when M == 0 (one other target in an empty neighbourhood
              never vetoes on its own).

Intermodulation relatives (optional, `relatives()`): blc1 was explained by
lookalikes at integer multiples of a 2.000004 MHz clock from itself (Sheikh et
al. 2021). For a list of clock spacings, it counts SAME-session, time-overlapping hits at
f + k*spacing (|k| <= KMAX) within 100 Hz + 4 Hz*|k| and compatible drift. The
chance expectation is measured in a window AT each relative position (so it
stays honest when the hit tables cover only narrow windows, as for blc1).
Off by default.

Schema (one row per hit; rfi_seed.csv.gz is this table, gzip CSV):
    freq_mhz     topocentric frequency at the hit's t=0 (turboSETI start)
    drift_hz_s   Hz/s, positive = rising frequency
    snr          as published; snr_units says which: "turboseti" (turboSETI's own,
                 ~4x below radiometer SNR on raw GBT data, KNOWN_ISSUES #1) or
                 "true" (radiometer SNR, inject.radiometer_noise definition)
    target, session (telescope:MJD-day), mjd (scan start), telescope,
    product, res_hz, tobs_s, drift_valid (0/1), source (file or paper)

    python rfi_db.py selftest                 # synthetic, offline, no downloads
    python rfi_db.py info [DB]                # what is in the database
    python rfi_db.py query FREQ DRIFT [--mjd M --target T --db DB|full]
    python rfi_db.py build OUT.npz SOURCES... # combine seed + extra sources
    python rfi_db.py explain EVENTS.csv [--db DB] [--exclude-session S]
    python rfi_db.py coverage [--db DB]       # how much band the veto removes
    python rfi_db.py blc1 [--db DB]           # the blc1 end-to-end control

DB: by default ONLY rfi_seed.csv.gz next to this file (29,202 hits; weak - it
explains 0.9 % of GJ699 filter-2 rows). `--db full` adds data/rfi_db_enriquez.npz
(Enriquez et al. 2017, 28.9 M hits; `tools/rfidb_reduce_enriquez.py`), which
explains 40.7 %; loading it peaks at 1.76 GB resident (measured) and is refused
when MemAvailable is too low. The seed is DERIVED DATA from published tables,
not covered by the MIT code licence: see RFI_DB_RESULT.md, "Data sources".
"""
import argparse
import csv
import gzip
import io
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SEED = HERE / "rfi_seed.csv.gz"
ENRIQUEZ_NPZ = HERE / "data" / "rfi_db_enriquez.npz"
C_M_S = 299792458.0

K_RES = 2.0
V_SLACK = 10.0          # m/s
DRIFT_ABS = 0.02        # Hz/s
DRIFT_REL = 0.20
NEIGH_HZ = 200e3        # chance-rate neighbourhood (+-)
SCORE_VETO = 3.0        # p <= 1e-3
CLOCKS_MHZ = (2.000004,)    # Sheikh et al. 2021 Parkes clock
KMAX = 200

COLUMNS = ["freq_mhz", "drift_hz_s", "snr", "snr_units", "target", "session", "mjd",
           "telescope", "product", "res_hz", "tobs_s", "drift_valid", "source"]
NUM = {"freq_mhz": np.float64, "drift_hz_s": np.float32, "snr": np.float32,
       "mjd": np.float64, "res_hz": np.float32, "tobs_s": np.float32, "drift_valid": np.int8}
STR = [c for c in COLUMNS if c not in NUM]

_SUFFIX = re.compile(r"(_(off|on|s|r|nd|a|b|c))+$", re.I)


def family(target):
    """Base name of a target: GJ699_OFF -> gj699, ProxCen_S -> proxcen,
    1421-490_ND_S -> 1421-490."""
    t = str(target).strip().lower().replace(" ", "")
    return _SUFFIX.sub("", t)


# ----------------------------------------------------------------------------- table
def _code_dtype(n):
    """Smallest signed integer type that can index a vocabulary of n strings."""
    return np.int8 if n < 127 else np.int16 if n < 32767 else np.int32


class RFIDB:
    """Column store, sorted by frequency. String columns are held as small integer
    codes into a vocabulary. Memory: ~43 bytes per hit (29 M hits ~1.25 GB of
    columns); measured peak while loading the full DB is in RFI_DB_RESULT.md.

    The constructor CONSUMES `cols` (pops each column) so that a large table is
    not held twice while it is converted and sorted."""

    def __init__(self, cols):
        n = len(cols["freq_mhz"])
        self.n = n
        self.c = {}
        for k, t in NUM.items():
            self.c[k] = np.asarray(cols.pop(k), dtype=t)
        self.vocab = {}
        for k in STR:
            v = cols.pop(k)
            if isinstance(v, tuple):                 # (codes, vocab)
                voc = [str(x) for x in v[1]]
                self.c[k], self.vocab[k] = np.asarray(v[0]).astype(_code_dtype(len(voc)), copy=False), voc
            else:
                u, inv = np.unique(np.asarray(v, dtype=str), return_inverse=True)
                self.c[k], self.vocab[k] = inv.astype(_code_dtype(len(u))), [str(x) for x in u]
            del v
        fam = [family(t) for t in self.vocab["target"]]
        fu, finv = np.unique(np.asarray(fam, dtype=str), return_inverse=True)
        self.fam_vocab = list(fu)
        self.c["family"] = finv.astype(_code_dtype(len(fu)))[self.c["target"]]
        f = self.c["freq_mhz"]
        if n > 1 and not bool(np.all(f[1:] >= f[:-1])):      # saved DBs are already sorted
            o = np.argsort(f, kind="stable")
            del f
            for k in list(self.c):                           # one column at a time
                self.c[k] = self.c[k][o]
            del o
        # widest frequency track of any hit (MHz); tracks themselves are computed per slice
        self.max_track = float((np.abs(self.c["drift_hz_s"]) * self.c["tobs_s"]).max()) * 1e-6 if n else 0.0

    def _track(self, sl):
        """(lo, hi) MHz of the frequency tracks [f, f + drift*T] of rows sl."""
        fh = self.c["freq_mhz"][sl]
        tr = self.c["drift_hz_s"][sl].astype(np.float64) * self.c["tobs_s"][sl] * 1e-6
        return fh + np.minimum(tr, 0), fh + np.maximum(tr, 0)

    # --- construction
    @classmethod
    def from_rows(cls, rows):
        cols = {k: [r[k] for r in rows] for k in COLUMNS}
        return cls(cols)

    @classmethod
    def concat(cls, dbs):
        """One DB from several. CONSUMES the inputs (their columns are freed as
        they are copied), so the big table is never held twice."""
        dbs = [d for d in dbs if d is not None and d.n]
        cols = {}
        for k in NUM:
            cols[k] = np.concatenate([d.c.pop(k) for d in dbs])
        for k in STR:
            vocab = sorted({str(v) for d in dbs for v in d.vocab[k]})
            idx = {v: i for i, v in enumerate(vocab)}
            dt = _code_dtype(len(vocab))
            cols[k] = (np.concatenate([np.array([idx[str(v)] for v in d.vocab[k]], dt)[d.c.pop(k)]
                                       for d in dbs]), vocab)
        for d in dbs:
            d.c.clear()
            d.n = 0
        return cls(cols)

    @classmethod
    def load(cls, path):
        path = Path(path)
        if path.suffix == ".npz":
            z = np.load(path, allow_pickle=False)
            cols = {k: z[k] for k in NUM}
            for k in STR:
                cols[k] = (z[k], [str(s) for s in z[k + "__vocab"]])
            return cls(cols)
        op = gzip.open if path.suffix == ".gz" else open
        with op(path, "rt", newline="") as f:
            rows = list(csv.DictReader(f))
        cols = {k: [r[k] for r in rows] for k in COLUMNS}
        for k in NUM:
            cols[k] = np.array(cols[k], dtype=np.float64).astype(NUM[k])
        return cls(cols)

    def save(self, path):
        path = Path(path)
        if path.suffix == ".npz":
            out = {k: self.c[k] for k in NUM}
            for k in STR:
                out[k] = self.c[k]
                out[k + "__vocab"] = np.array(self.vocab[k], dtype=str)
            np.savez_compressed(path, **out)
            return
        buf = io.StringIO()
        w = csv.writer(buf, lineterminator="\n")
        w.writerow(COLUMNS)
        for i in range(self.n):
            row = []
            for k in COLUMNS:
                if k in NUM:
                    v = self.c[k][i]
                    row.append(f"{v:.6f}" if k == "freq_mhz" else f"{v:.6f}" if k == "mjd"
                               else f"{v:.5f}" if k == "drift_hz_s" else f"{v:.4g}" if k != "drift_valid"
                               else str(int(v)))
                else:
                    row.append(self.vocab[k][self.c[k][i]])
            w.writerow(row)
        data = buf.getvalue().encode()
        if path.suffix == ".gz":
            with open(path, "wb") as f:                     # mtime=0: diffable, reproducible
                with gzip.GzipFile(fileobj=f, mode="wb", mtime=0, filename="") as g:
                    g.write(data)
        else:
            path.write_bytes(data)

    def str_of(self, k, i):
        return self.vocab[k][self.c[k][i]]

    # --- matching
    def candidates(self, f_mhz, pad_mhz):
        a = np.searchsorted(self.c["freq_mhz"], f_mhz - pad_mhz - self.max_track, "left")
        b = np.searchsorted(self.c["freq_mhz"], f_mhz + pad_mhz + self.max_track, "right")
        return a, b

    def _drift_ok(self, sl, d_q, res_q, tobs_q):
        dh = self.c["drift_hz_s"][sl].astype(np.float64)
        tol = (DRIFT_ABS + DRIFT_REL * np.maximum(abs(d_q), np.abs(dh))
               + res_q / tobs_q + self.c["res_hz"][sl] / self.c["tobs_s"][sl])
        return np.abs(dh - d_q) <= tol

    def query(self, f_mhz, drift, mjd=None, target=None, tobs_s=300.0, res_hz=2.8,
              exclude_sessions=(), only_sessions=None, v_slack=V_SLACK, k_res=K_RES):
        """Verdict dict for one hit (see module docstring for the rule)."""
        fam = family(target) if target is not None else None
        famcode = self.fam_vocab.index(fam) if fam in self.fam_vocab else -1
        slack = (k_res * res_hz + f_mhz * 1e6 * v_slack / C_M_S) * 1e-6     # MHz (q side)
        tr = drift * tobs_s * 1e-6
        qlo, qhi = f_mhz + min(tr, 0.0), f_mhz + max(tr, 0.0)
        a, b = self.candidates(f_mhz, NEIGH_HZ * 1e-6 + abs(tr))
        sl = slice(a, b)
        other = self.c["family"][sl] != famcode
        if exclude_sessions or only_sessions is not None:
            ses = np.array(self.vocab["session"], dtype=object)[self.c["session"][sl]]
            if exclude_sessions:
                other &= ~np.isin(ses, list(exclude_sessions))
            if only_sessions is not None:
                other &= np.isin(ses, list(only_sessions))
        valid = self.c["drift_valid"][sl] == 1
        dok = self._drift_ok(sl, drift, res_hz, tobs_s)
        sl_h = (k_res * self.c["res_hz"][sl] * 1e-6)
        lo_h, hi_h = self._track(sl)
        overlap = (lo_h - sl_h <= qhi + slack) & (hi_h + sl_h >= qlo - slack)
        m = other & valid & dok & overlap
        fo = other & ~valid & overlap      # unresolved-drift hits: information only
        # chance: drift-compatible other-target hits in the +-NEIGH neighbourhood
        # (the matching hits themselves are left out: the background is the annulus)
        neigh = other & valid & dok & ~m & (np.abs(self.c["freq_mhz"][sl] - f_mhz) <= NEIGH_HZ * 1e-6)
        M = int(neigh.sum())
        if M:
            mtrack, msl = float(np.mean(hi_h[neigh] - lo_h[neigh])), float(np.mean(sl_h[neigh]))
        else:                       # no background measured: assume hits like the query
            mtrack, msl = qhi - qlo, k_res * res_hz * 1e-6
        width = (qhi - qlo) + mtrack + 2 * slack + 2 * msl
        # rate floor: a rate measured from M = 0 background hits is NOT zero; count it as 1
        lam = max(M, 1) * width / max(2 * NEIGH_HZ * 1e-6 - width, width)
        idx = np.nonzero(m)[0] + a
        tcodes = np.unique(self.c["target"][idx])
        fams = np.unique(self.c["family"][idx])
        scodes = np.unique(self.c["session"][idx])
        n_t = int(len(fams))
        p = _poisson_sf(n_t, lam) if n_t else 1.0
        score = min(10.0, -math.log10(max(p, 1e-10))) if n_t else 0.0
        # with no measured background, one other target is never enough on its own
        veto = bool(n_t >= (2 if M == 0 else 1) and score >= SCORE_VETO)
        return dict(n_hits=int(len(idx)), n_targets=n_t, n_sessions=int(len(scodes)),
                    targets=[self.vocab["target"][t] for t in tcodes][:8],
                    sessions=[self.vocab["session"][s] for s in scodes][:8],
                    sources=sorted({self.vocab["source"][self.c["source"][i]].split(" ")[0]
                                    for i in idx[:200]}),
                    n_freq_only=int(len(np.unique(self.c["family"][np.nonzero(fo)[0] + a]))),
                    lam=lam, n_background=M, p=p, score=score, veto=veto, idx=idx)

    def relatives(self, f_mhz, drift, session, tobs_s=300.0, res_hz=2.8,
                  spacings=CLOCKS_MHZ, kmax=KMAX, tol0_hz=100.0, tol_k_hz=4.0, mjd=None,
                  neigh_hz=NEIGH_HZ, ks=None, offset=0.0):
        """Same-session hits at f + k*spacing (any target, k != 0) with compatible
        drift: evidence of an intermodulation family (Sheikh et al. 2021). With
        mjd given, only hits from scans that OVERLAP the query scan in time count
        (an intermodulation family exists at one moment, in one receiver chain).

        Chance rate, SELECTION-AWARE: around each position f + k*spacing the
        density of eligible (same session, time-overlapping, drift-compatible)
        hits is measured in a +-neigh_hz window AT THAT POSITION, and
            lambda = sum_k  M_k * 2 tol_k / (2 neigh_hz)
        where M_k counts every eligible hit in the window, the matching one
        included (so a lone match in an otherwise empty window costs
        2 tol_k / 2 neigh_hz, never zero). If the hit tables only cover narrow
        windows around the relatives (blc1: +-5 kHz fetched), pass neigh_hz =
        the covered half-width; a band-wide density would be far too low.

        ks restricts to given k values; offset shifts every position by
        offset*spacing (offset 0.5 = half-clock CONTROL positions, which are no
        relative of f at all)."""
        if not 0 < neigh_hz < 0.5e6 * min(spacings):
            raise ValueError("neigh_hz must be below half the clock spacing")
        ses_code = self.vocab["session"].index(session) if session in self.vocab["session"] else -1
        found = []
        lam = 0.0
        for sp in spacings:
            span = (kmax + abs(offset)) * sp
            a, b = self.candidates(f_mhz, span + 1.0)
            sl = slice(a, b)
            same = (self.c["session"][sl] == ses_code) & (self.c["drift_valid"][sl] == 1)
            if mjd is not None:
                t0 = self.c["mjd"][sl]
                same &= (t0 <= mjd + tobs_s / 86400 + 1e-4) & \
                    (t0 + self.c["tobs_s"][sl] / 86400 >= mjd - 1e-4)
            dok = self._drift_ok(sl, drift, res_hz, tobs_s)
            pool = same & dok
            fr = self.c["freq_mhz"][sl]
            k = np.round((fr - f_mhz) / sp - offset)
            dist = np.abs(fr - f_mhz - (k + offset) * sp)
            tol = (tol0_hz + tol_k_hz * np.abs(k + offset)) * 1e-6
            ok_k = ((k + offset) != 0) & (np.abs(k) <= kmax)
            if ks is not None:
                ok_k &= np.isin(k, list(ks))
            hit = pool & ok_k & (dist <= tol)
            for i in np.nonzero(hit)[0]:
                found.append((sp, int(k[i]), float(fr[i]), float(self.c["drift_hz_s"][a + i]),
                              self.vocab["target"][self.c["target"][a + i]]))
            inwin = pool & ok_k & (dist <= neigh_hz * 1e-6)
            lam += float(np.sum(2 * tol[inwin] / (2 * neigh_hz * 1e-6)))
        kset = sorted({(sp, kk) for sp, kk, *_ in found})
        p = _poisson_sf(len(kset), lam) if kset else 1.0
        return dict(n_k=len(kset), found=found, lam=lam, p=p,
                    score=min(10.0, -math.log10(max(p, 1e-10))) if kset else 0.0)


def _poisson_sf(n, lam):
    """P(N >= n) for N ~ Poisson(lam)."""
    if n <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    from scipy.stats import poisson
    return float(poisson.sf(n - 1, lam))


def explain_line(v):
    """One plain-language line for cadence_verdict.py."""
    if v["n_targets"] == 0:
        extra = f" ({v['n_freq_only']} other target(s) at this frequency with unresolved drift)" \
            if v["n_freq_only"] else ""
        return None, "RFI database: no other target shows this frequency+drift" + extra
    s = (f"RFI database: seen toward {v['n_targets']} OTHER target(s) in {v['n_sessions']} "
         f"session(s) ({', '.join(v['targets'][:3])}{'...' if v['n_targets'] > 3 else ''}); "
         f"chance {v['lam']:.2g} expected, score {v['score']:.1f}")
    return (s + " - local interference") if v["veto"] else None, \
        s + (" - VETO" if v["veto"] else " - not significant")


# ----------------------------------------------------------------------------- loaders
def rows_turboseti_dat(path, telescope="GBT", source=None, session=None):
    """turboSETI .dat -> schema rows (snr_units turboseti)."""
    hdr = {}
    rows = []
    for ln in open(path):
        if ln.startswith("#"):
            for key, pat in (("source", r"Source:\s*(\S+)"), ("mjd", r"MJD:\s*([\d.]+)"),
                             ("df", r"DELTAF\(Hz\):\s*(-?[\d.]+)"),
                             ("tobs", r"obs_length:\s*([\d.]+)")):
                mm = re.search(pat, ln)
                if mm:
                    hdr[key] = mm.group(1)
            continue
        p = ln.split()
        if len(p) < 4:
            continue
        df = abs(float(hdr["df"]))
        mjd = float(hdr["mjd"])
        rows.append(dict(freq_mhz=float(p[3]), drift_hz_s=float(p[1]), snr=float(p[2]),
                         snr_units="turboseti", target=hdr["source"],
                         session=session or f"{telescope}:{int(mjd)}", mjd=mjd,
                         telescope=telescope,
                         product="fine 0000" if df < 10 else ("mid 0002" if df > 1000 else "other"),
                         res_hz=df, tobs_s=float(hdr["tobs"]), drift_valid=int(df < 10),
                         source=source or Path(path).name))
    return rows


_NUMS = re.compile(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?")


def rows_price2020(path, telescope="GBT", tobs=300.0, res=2.79):
    """Price et al. 2020 (AJ 159, 86) event CSVs (lband2019_events.csv,
    pks10cm2019_events.csv): one row per constituent hit."""
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            fid = r["FileID"]
            mm = re.search(r"guppi_(\d{5})_(\d{5})_", fid)
            if mm:
                mjd = int(mm.group(1)) + int(mm.group(2)) / 86400
            else:                        # e.g. ..._guppi_57403_HIP32769_0018 (day only)
                mm = re.search(r"(?:guppi|Parkes)_(\d{5})_", fid) or re.search(r"_(5\d{4})_", fid)
                if not mm:
                    continue
                mjd = int(mm.group(1)) + 0.5
            fr = [float(x) for x in _NUMS.findall(r["Freqs"])]
            dr = [float(x) for x in _NUMS.findall(r["DriftRates"])]
            sn = [float(x) for x in _NUMS.findall(r["SNR"])]
            for i, f0 in enumerate(fr):
                rows.append(dict(freq_mhz=f0, drift_hz_s=dr[min(i, len(dr) - 1)],
                                 snr=sn[min(i, len(sn) - 1)], snr_units="turboseti",
                                 target=r["Source"], session=f"{telescope}:{int(mjd)}",
                                 mjd=mjd, telescope=telescope, product="fine 0000", res_hz=res,
                                 tobs_s=tobs, drift_valid=1,
                                 source=f"Price2020 {Path(path).name}"))
    return rows


def enriquez_db(chunk_dir, res=2.79, tobs=300.0):
    """Enriquez et al. 2017 (ApJ 849, 104) all 28.9 M hits, from the npz chunks
    written by tools/rfidb_reduce_enriquez.py. Strings are encoded chunk by chunk
    (a 29 M-row string array would need ~2 GB)."""
    fs = sorted(Path(chunk_dir).glob("chunk*.npz"))
    freq, drift, snr, mjd, tcode = [], [], [], [], []
    tvocab = {}
    for f in fs:
        z = np.load(f)
        u, inv = np.unique(z["source"], return_inverse=True)
        m = np.array([tvocab.setdefault(str(t), len(tvocab)) for t in u], np.int32)
        tcode.append(m[inv])
        freq.append(z["freq"]); drift.append(z["drift"]); snr.append(z["snr"]); mjd.append(z["mjd"])
    freq, mjd = np.concatenate(freq), np.concatenate(mjd)
    n = len(freq)
    days = np.floor(mjd).astype(np.int64)
    du, sinv = np.unique(days, return_inverse=True)
    one = lambda s: (np.zeros(n, np.int32), [s])       # noqa: E731
    cols = dict(freq_mhz=freq, drift_hz_s=np.concatenate(drift), snr=np.concatenate(snr), mjd=mjd,
                res_hz=np.full(n, res, np.float32), tobs_s=np.full(n, tobs, np.float32),
                drift_valid=np.ones(n, np.int8), snr_units=one("turboseti"),
                target=(np.concatenate(tcode), sorted(tvocab, key=tvocab.get)),
                session=(sinv.astype(np.int32), [f"GBT:{d}" for d in du]),
                telescope=one("GBT"), product=one("fine 0000"),
                source=one("Enriquez2017 All_hits_turbo_seti.csv"))
    return RFIDB(cols)


# Measured peak resident memory while loading seed + Enriquez npz and concatenating
# (PC, 2026-11-07, see RFI_DB_RESULT.md) divided by the npz's uncompressed size.
PEAK_PER_NPZ_BYTE = 1.25


def mem_available():
    """MemAvailable in bytes (NOT MemFree, which ignores reclaimable cache), or
    None where /proc/meminfo does not exist."""
    try:
        for ln in open("/proc/meminfo"):
            if ln.startswith("MemAvailable:"):
                return int(ln.split()[1]) * 1024
    except OSError:
        pass
    return None


def ram_needed(paths):
    """Estimated peak bytes to load these DB files (npz: uncompressed size x factor)."""
    import zipfile
    need = 0
    for p in paths:
        p = Path(p)
        if p.suffix == ".npz":
            with zipfile.ZipFile(p) as z:
                need += PEAK_PER_NPZ_BYTE * sum(i.file_size for i in z.infolist())
        else:
            need += 30 * p.stat().st_size      # gz CSV: rows -> python dicts, small files only
    return int(need)


def default_db(path=None):
    """DB to use. No path: the checked-in SEED ONLY (weak: it has few hits besides
    GJ699's own, which are left out for GJ699). "full": seed + the Enriquez 2017
    npz (strong; needs the 10.8 GB stream reduced on a big machine and several GB
    of RAM). Otherwise a comma-separated list of files. Big files are refused when
    MemAvailable is below the estimated peak (set SETI_COOKBOOK_RFIDB_FORCE=1 to
    load anyway)."""
    import os
    if not path:
        paths = [SEED]
    elif str(path) == "full":
        if not ENRIQUEZ_NPZ.exists():
            raise SystemExit(f"rfi_db: --db full needs {ENRIQUEZ_NPZ} "
                             "(build it with tools/rfidb_reduce_enriquez.py + rfi_db.py build)")
        paths = [SEED, ENRIQUEZ_NPZ]
    else:
        paths = [Path(p) for p in str(path).split(",")]
    need, have = ram_needed(paths), mem_available()
    if have is not None and need > 0.8 * have and not os.environ.get("SETI_COOKBOOK_RFIDB_FORCE"):
        raise SystemExit(f"rfi_db: loading {', '.join(Path(p).name for p in paths)} needs ~{need / 1e9:.1f} GB "
                         f"but MemAvailable is {have / 1e9:.1f} GB; refusing (SETI_COOKBOOK_RFIDB_FORCE=1 overrides)")
    dbs = [RFIDB.load(p) for p in paths]
    return dbs[0] if len(dbs) == 1 else RFIDB.concat(dbs)


def describe(db):
    """One line saying which DB is in use and how strong it is."""
    srcs = {db.vocab["source"][c].split(" ")[0] for c in np.unique(db.c["source"])}
    if "Enriquez2017" in srcs:
        return f"RFI database: FULL ({db.n:,} hits, {len(db.fam_vocab):,} target families, incl. Enriquez 2017)"
    return (f"RFI database: NO Enriquez 2017 ({db.n:,} hits, {len(db.fam_vocab):,} target families) - weak: "
            "the seed alone explains 0.9 % of GJ699 filter-2 rows vs 40.7 % for the full DB ('full')")


# ----------------------------------------------------------------------------- allocations
def allocations():
    """recipe_api.RFI_BANDS merged unchanged, with review notes (nothing deleted)."""
    import recipe_api as R
    notes = {
        (1435.0, 1535.0): "QUESTIONABLE upper edge: US AMT is 1435-1525 MHz (FCC 47 CFR 2.106, "
                          "footnote US78); 1525-1535 is MSS downlink, already in the Inmarsat entry",
        (1610.0, 1626.5): "QUESTIONABLE label: 1610-1626.5 is the MSS uplink band (Globalstar "
                          "1610-1618.725 up); Iridium is TDD 1616-1626.5 (up AND down). Range fine.",
        (1215.0, 1400.0): "GBT L-band receiver notch filter covers ~1200-1341 MHz; above 1341 "
                          "the radar label is a broad assumption",
        (2180.0, 2290.0): "overlaps the next entry (2200-2300); harmless",
    }
    return [(lo, hi, why, notes.get((lo, hi), "")) for lo, hi, why in R.RFI_BANDS]


def band_fraction(lo=1100.0, hi=1900.0):
    """Fraction of [lo, hi] MHz covered by RFI_BANDS (union)."""
    iv = sorted((max(a, lo), min(b, hi)) for a, b, *_ in allocations() if b > lo and a < hi)
    tot, cur = 0.0, None
    for a, b in iv:
        if cur is None or a > cur[1]:
            if cur:
                tot += cur[1] - cur[0]
            cur = [a, b]
        else:
            cur[1] = max(cur[1], b)
    if cur:
        tot += cur[1] - cur[0]
    return tot / (hi - lo)


# ----------------------------------------------------------------------------- commands
def cmd_info(db):
    print(describe(db))
    print(f"{db.n:,} hits")
    for k in ("source", "telescope", "snr_units", "product"):
        codes, cnt = np.unique(db.c[k], return_counts=True)
        tops = sorted(zip(cnt, codes), reverse=True)[:12]
        print(f"  by {k}: " + "; ".join(f"{db.vocab[k][c]} {n:,}" for n, c in tops))
    print(f"  targets {len(db.vocab['target']):,} (families {len(db.fam_vocab):,}), "
          f"sessions {len(db.vocab['session']):,}")
    print(f"  freq {db.c['freq_mhz'].min():.3f} - {db.c['freq_mhz'].max():.3f} MHz; "
          f"drift-unresolved hits {int((db.c['drift_valid'] == 0).sum()):,}")
    print(f"RFI_BANDS cover {100 * band_fraction():.1f} % of 1.1-1.9 GHz")


def cmd_query(db, a):
    v = db.query(a.freq, a.drift, mjd=a.mjd, target=a.target, tobs_s=a.tobs, res_hz=a.res)
    _, line = explain_line(v)
    print(line)
    print(f"  n_hits {v['n_hits']}  n_targets {v['n_targets']}  n_sessions {v['n_sessions']}  "
          f"lambda {v['lam']:.3g}  p {v['p']:.3g}  veto {v['veto']}")
    if a.session:
        r = db.relatives(a.freq, a.drift, a.session, a.tobs, a.res, mjd=a.mjd)
        print(f"  intermod relatives in {a.session}: {r['n_k']} (chance {r['lam']:.2g}, "
              f"score {r['score']:.1f}) " + " ".join(f"k={k:+d}:{f:.6f}" for _, k, f, *_ in r["found"][:10]))
    return 0


def read_events(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def group_events(ev):
    """Group find_event rows (one per ON-scan hit) into events: rows from different
    ON tables are one event when the later one lies on the earlier one's drift
    track (|f_j - f_i - d*(t_j - t_i)| <= 2 channels + f*10 m/s/c + 0.02 Hz/s*dt)
    with compatible drift. Uses Hit_ID when the table has it. Returns a group id per row."""
    if all(e.get("Hit_ID") for e in ev):
        u = {}
        return [u.setdefault(e["Hit_ID"], len(u)) for e in ev]
    n = len(ev)
    par = list(range(n))

    def root(i):
        while par[i] != i:
            par[i] = par[par[i]]
            i = par[i]
        return i
    f = np.array([float(e["Freq"]) for e in ev])
    d = np.array([float(e["DriftRate"]) for e in ev])
    t = np.array([float(e.get("delta_t") or 0) for e in ev])
    tab = [e.get("status", "") for e in ev]
    r = np.array([abs(float(e.get("DELTAF", 2.8) or 2.8)) for e in ev])
    o = np.argsort(f)
    for ii, i in enumerate(o):
        for j in o[ii + 1:]:
            if f[j] - f[i] > 4 * 4.0 * 1000 * 1e-6:        # beyond any 4 Hz/s x 1000 s track
                break
            if tab[i] == tab[j]:
                continue
            dt = t[j] - t[i]
            tol = (2 * max(r[i], r[j]) + f[i] * 1e6 * V_SLACK / C_M_S + DRIFT_ABS * abs(dt)) * 1e-6
            dtol = DRIFT_ABS + DRIFT_REL * max(abs(d[i]), abs(d[j])) + 2 * max(r[i], r[j]) / 300
            if abs(f[j] - f[i] - d[i] * dt * 1e-6) <= tol and abs(d[i] - d[j]) <= dtol:
                par[root(j)] = root(i)
    roots = {}
    return [roots.setdefault(root(i), len(roots)) for i in range(n)]


def cmd_explain(db, a):
    """Leave-the-target-out explanation of a find_event_pipeline CSV."""
    import recipe_api as R
    print(describe(db))
    ev = read_events(a.events)
    excl = set(a.exclude_session or [])
    n_v = n_band = n_v_noband = 0
    vetoed = []
    for e in ev:
        f, d = float(e["Freq"]), float(e["DriftRate"])
        v = db.query(f, d, mjd=float(e["MJD"]), target=e["Source"],
                     tobs_s=float(e.get("obs_length", 300) or 300),
                     res_hz=abs(float(e.get("DELTAF", 2.8) or 2.8)), exclude_sessions=excl)
        band = R.rfi_reason(f) is not None
        n_v += v["veto"]
        vetoed.append(v["veto"])
        n_band += band
        n_v_noband += v["veto"] and not band
        if a.verbose:
            print(f"{f:12.6f} {d:+.4f} veto={int(v['veto'])} n_t={v['n_targets']} "
                  f"n_s={v['n_sessions']} lam={v['lam']:.2g} band={int(band)} {v['targets'][:3]}")
    nb = len(ev) - n_band
    print(f"{a.events}: {len(ev)} rows; excluded sessions {sorted(excl) or 'none'}")
    print(f"  explained by OTHER targets (veto): {n_v}/{len(ev)} = {100 * n_v / max(1, len(ev)):.1f} %")
    print(f"  in an RFI_BANDS allocation: {n_band}; outside any band: {nb}, of which "
          f"DB-vetoed {n_v_noband} ({100 * n_v_noband / max(1, nb):.1f} %)")
    g = np.array(group_events(ev))
    ng = int(g.max()) + 1 if len(g) else 0
    vv = np.array(vetoed, bool)
    any_v = sum(bool(vv[g == k].any()) for k in range(ng))
    all_v = sum(bool(vv[g == k].all()) for k in range(ng))
    print(f"  grouped into {ng} events (rows on one drift track across ON tables): "
          f"{any_v} = {100 * any_v / max(1, ng):.1f} % have a vetoed row, "
          f"{all_v} = {100 * all_v / max(1, ng):.1f} % have every row vetoed")
    return 0


def cmd_coverage(db, a):
    """Fraction of (frequency, drift) space the veto removes, by Monte-Carlo probes."""
    rng = np.random.default_rng(a.seed)
    print(describe(db))
    lo, hi = a.lo, a.hi
    print(f"veto coverage, {a.n} random probes per drift in {lo:g}-{hi:g} MHz, "
          f"probe target 'PROBE' (no target excluded), T=300 s, r=2.8 Hz")
    print(f"  RFI_BANDS remove {100 * band_fraction(lo, hi):.1f} % of the band (frequency only)")
    for d in a.drifts:
        f = rng.uniform(lo, hi, a.n)
        dd = d if d != 0 else 0.0
        nv = nf = 0
        for x in f:
            v = db.query(x, dd, target="PROBE")
            nv += v["veto"]
            nf += v["n_targets"] > 0
        print(f"  drift {d:+5.2f} Hz/s: vetoed {100 * nv / a.n:6.3f} %   "
              f"(any other-target match, before chance test {100 * nf / a.n:6.3f} %)")
    # frequency-only comparison: ignore drift (what a freq-only DB would remove)
    global DRIFT_ABS
    keep = DRIFT_ABS
    DRIFT_ABS = 1e3
    try:
        f = rng.uniform(lo, hi, a.n)
        nv = sum(db.query(x, 0.5, target="PROBE")["veto"] for x in f)
    finally:
        DRIFT_ABS = keep
    print(f"  FREQUENCY-ONLY keying (drift ignored), probe drift +0.5: vetoed {100 * nv / a.n:6.3f} %")
    return 0


BLC1_HALF_KHZ = 5.0      # tools/rfidb_blc1_fetch.py fetches +-5 kHz per window
BLC1_DRIFT_MAX = 0.5     # ... and its de-drift search cuts drift_max*T at each edge
BLC1_KS = (-135, -63, 40)


def cmd_blc1(db, a):
    """End-to-end control on the public blc1 data (rows from tools/rfidb_blc1_fetch.py).

    The Parkes data were fetched ONLY in +-5 kHz windows at blc1 and at the three
    lookalikes the paper names, so the intermod chance rate must be measured inside
    those windows (relatives(neigh_hz = covered half-width)), not spread over the
    band. If the half-clock CONTROL windows are in the DB too (fetch --windows ctl),
    the empirical rate there is printed as well."""
    blc1 = 982.0026
    sel = [i for i in range(db.n) if family(db.str_of("target", i)) == "proxcen"
           and abs(db.c["freq_mhz"][i] - blc1) < 0.0005 and db.c["drift_hz_s"][i] > 0.005]
    if not sel:
        print("no blc1 rows in this DB (run tools/rfidb_blc1_fetch.py and rebuild)")
        return 1
    print(f"blc1 control: {len(sel)} ProxCen hits within 500 Hz of 982.0026 MHz, drift > 0")
    print("  intermod chance rate measured INSIDE the fetched windows (the data are a selection)")
    ok = present = 0
    ctl_pos = ctl_hit = rel_pos = rel_hit = 0
    for i in sel:
        f, d = float(db.c["freq_mhz"][i]), float(db.c["drift_hz_s"][i])
        T, res, mjd = float(db.c["tobs_s"][i]), float(db.c["res_hz"][i]), float(db.c["mjd"][i])
        ses = db.str_of("session", i)
        v = db.query(f, d, target="ProxCen", tobs_s=T, res_hz=res)
        v2 = db.query(f, d, target="ProxCen", tobs_s=T, res_hz=res, exclude_sessions={ses})
        half = BLC1_HALF_KHZ * 1e3 - BLC1_DRIFT_MAX * T - res      # covered half-width, Hz
        r = db.relatives(f, d, ses, T, res, mjd=mjd, neigh_hz=half)
        r_naive = db.relatives(f, d, ses, T, res, mjd=mjd)          # band-wide density (wrong here)
        kk = sorted({k for _, k, *_ in r["found"]})
        present += set(BLC1_KS) <= set(kk)
        rel_pos += len(BLC1_KS)
        nc = 0
        for off in (-0.5, 0.5):
            rc = db.relatives(f, d, ses, T, res, mjd=mjd, neigh_hz=half, ks=[k for k in BLC1_KS],
                              offset=off)
            nc += len({k for _, k, *_ in rc["found"]})
        ctl_pos += 2 * len(BLC1_KS)
        ctl_hit += nc
        nr = len(set(kk) & set(BLC1_KS))
        rel_hit += nr
        from scipy.stats import fisher_exact
        pf = fisher_exact([[nr, len(BLC1_KS) - nr], [nc, 2 * len(BLC1_KS) - nc]], alternative="greater")[1]
        flag = v["veto"] or r["score"] >= SCORE_VETO
        ok += flag
        print(f"  {mjd:.4f} {f:.6f} MHz {d:+.4f} Hz/s SNR {db.c['snr'][i]:.0f}: "
              f"other targets {v['n_targets']} {v['targets'][:4]} score {v['score']:.1f}; "
              f"other-session-only {v2['n_targets']} score {v2['score']:.1f}; "
              f"relatives k={kk} chance {r['lam']:.2g} score {r['score']:.1f} "
              f"(band-wide density would say {r_naive['lam']:.2g}, score {r_naive['score']:.1f}); "
              f"relatives {nr}/3 vs half-clock controls {nc}/6 (Fisher p {pf:.2g}) "
              f"-> {'FLAGGED' if flag else 'not flagged'}")
    print(f"blc1: relatives present at all of k={list(BLC1_KS)} in {present}/{len(sel)} scans; "
          f"flagged at score >= {SCORE_VETO:g} in {ok}/{len(sel)} scans")
    c = 2.000004
    centres = [blc1 + (k + o) * c for k in BLC1_KS for o in (-0.5, 0.5)]
    fr = db.c["freq_mhz"]
    n_ctl_rows = sum(int(np.searchsorted(fr, x + BLC1_HALF_KHZ * 1e-3) - np.searchsorted(fr, x - BLC1_HALF_KHZ * 1e-3))
                     for x in centres)
    if n_ctl_rows:
        from scipy.stats import fisher_exact
        pp = fisher_exact([[rel_hit, rel_pos - rel_hit], [ctl_hit, ctl_pos - ctl_hit]], alternative="greater")[1]
        print(f"  empirical: relative positions with a matching hit {rel_hit}/{rel_pos}, half-clock CONTROL "
              f"positions {ctl_hit}/{ctl_pos} ({n_ctl_rows} DB hits inside the control windows); pooled "
              f"Fisher p {pp:.2g} - NOT independent tests (the scans share the same persistent lines)")
    else:
        print("  no hits inside the half-clock control windows: fetch them with "
              "tools/rfidb_blc1_fetch.py fetch CTL.csv --windows ctl")
    return 0


# ----------------------------------------------------------------------------- selftest
def _synthetic(seed=0):
    rng = np.random.default_rng(seed)
    rows = []

    def add(f, d, t, s, mjd, snr=30.0, res=2.8, T=300.0, valid=1, src="synthetic"):
        rows.append(dict(freq_mhz=f, drift_hz_s=d, snr=snr, snr_units="true", target=t,
                         session=s, mjd=mjd, telescope="GBT", product="fine", res_hz=res,
                         tobs_s=T, drift_valid=valid, source=src))
    # background: 20,000 random hits, random targets/sessions, random drifts
    for k in range(20000):
        add(rng.uniform(1100, 1900), rng.uniform(-4, 4), f"HIP{rng.integers(1, 400)}",
            f"GBT:{57400 + rng.integers(0, 30)}", 57400.5)
    # recurring zero-drift line at 1500.000100 MHz on 5 other targets, 3 sessions
    for j, (t, s) in enumerate([("HIPa", "GBT:57401"), ("HIPb", "GBT:57401"), ("HIPc", "GBT:57405"),
                                ("HIPd", "GBT:57409"), ("HIPe_OFF", "GBT:57409")]):
        add(1500.000100 + rng.normal(0, 3e-6), 0.0, t, s, 57401 + j)
    # recurring drifting RFI at 1600.0 MHz, -0.5 Hz/s, 4 targets
    for j, t in enumerate(["HIPf", "HIPg", "HIPh", "HIPi"]):
        add(1600.0 + rng.normal(0, 5e-6), -0.5 + rng.normal(0, 0.01), t, f"GBT:{57410 + j}", 57410 + j)
    # a broad persistent band: 3000 zero-drift hits densely over 1650.0-1650.1 MHz, 12 targets
    for k in range(3000):
        add(rng.uniform(1650.0, 1650.1), rng.normal(0, 0.005), f"BAND{k % 12}",
            f"GBT:{57420 + k % 4}", 57420.5)
    # the query star's OWN hits (other date + its OFF scan) at the unique frequency
    add(1420.300000, -0.12, "STARX", "GBT:57450", 57450.5)
    add(1420.300000, -0.12, "STARX_OFF", "GBT:57455", 57455.5)
    # mid-res, drift unresolved, at another frequency: must never veto
    for t in ("HIPm1", "HIPm2", "HIPm3"):
        add(1700.000000, 9.8, t, "GBT:57430", 57430.5, res=2861.0, valid=0)
    # blc1 known-answer family (published frequencies; Sheikh et al. 2021)
    c = 2.000004
    for k in (-135, -63, 40):
        add(982.0026 + k * c + rng.normal(0, 2e-5), 0.03, "PROXCEN_S", "Parkes:58602", 58602.6,
            res=3.81, T=1800.0)
    add(982.0026, 0.03, "PROXCEN_S", "Parkes:58602", 58602.6, res=3.81, T=1800.0)
    # relatives distractors that must NOT count: same clock family in ANOTHER session,
    # and in the same session but a scan that does not overlap the query scan in time
    add(982.0026 - 100 * c, 0.03, "PROXCEN_S", "Parkes:58598", 58598.6, res=3.81, T=1800.0)
    add(982.0026 + 77 * c, 0.03, "HIPq", "Parkes:58602", 58602.0, res=3.81, T=1800.0)
    # ONE other-target hit, alone in an empty neighbourhood (track-sign + rate-floor tests):
    # +0.5 Hz/s for 300 s, track [2000.000000, 2000.000150] MHz
    add(2000.0, 0.5, "HIPlone", "GBT:57440", 57440.5)
    add(2010.0, 0.0, "HIPlone0", "GBT:57440", 57440.5)     # zero drift: narrowest match, highest score
    for k in range(400):     # unrelated same-session hits for a chance rate
        add(rng.uniform(704, 1088), rng.uniform(-0.5, 0.5), "1421-490_S", "Parkes:58602", 58602.6,
            res=3.81, T=1800.0)
    return RFIDB.from_rows(rows)


def selftest():
    db = _synthetic()
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    v = db.query(1500.000100, 0.0, target="STARX")
    chk(f"recurring zero-drift line (5 other targets) is VETOED [n_t={v['n_targets']}, score {v['score']:.1f}]",
        v["veto"] and v["n_targets"] >= 4)
    v = db.query(1600.0, -0.5, target="STARX")
    chk(f"recurring drifting line (-0.5 Hz/s, 4 targets) is VETOED [n_t={v['n_targets']}]", v["veto"])
    v = db.query(1420.300000, -0.12, target="STARX")
    chk(f"unique signal SURVIVES although the star's own other-date + OFF hits match "
        f"[n_t={v['n_targets']}]", not v["veto"] and v["n_targets"] == 0)
    v_own = db.query(1420.300000, -0.12, target="HIPz")
    chk("...and the same hits DO count when the query is a different star (no hidden exclusion)",
        v_own["n_targets"] == 1 and v_own["n_hits"] == 2)
    # a +2 Hz/s signal whose track CROSSES the 1500.000100 line mid-scan
    v = db.query(1500.000100 - 300e-6, 2.0, target="STARX")
    chk(f"drifting signal (+2 Hz/s) crossing the persistent line is NOT vetoed [n_t={v['n_targets']}]",
        not v["veto"])
    v = db.query(1650.05, 1.0, target="STARX")
    chk(f"1 Hz/s signal inside the broad zero-drift band is NOT vetoed [n_t={v['n_targets']}]",
        not v["veto"])
    v = db.query(1650.05, 0.0, target="STARX")
    chk(f"zero-drift signal in the band IS matched but the chance test holds it back "
        f"[n_t={v['n_targets']}, lambda {v['lam']:.1f}, score {v['score']:.1f}]",
        v["n_targets"] >= 1 and not v["veto"])
    v = db.query(1500.000100 + 30e-6, 0.0, target="STARX")
    chk("Doppler/oscillator slack: same line 30 Hz away still matches", v["veto"])
    v = db.query(1500.000100 + 2e-3, 0.0, target="STARX")
    chk(f"2 kHz away does NOT match [n_t={v['n_targets']}]", v["n_targets"] == 0)
    v = db.query(1700.0, 9.8, target="STARX", res_hz=2861.0)
    chk(f"drift-unresolved (mid-res) hits never veto [n_freq_only={v['n_freq_only']}]",
        not v["veto"] and v["n_freq_only"] == 3)
    r = db.relatives(982.0026, 0.03, "Parkes:58602", 1800.0, 3.81, mjd=58602.6)
    ks = sorted(k for _, k, *_ in r["found"])
    chk(f"blc1 known-answer: intermod relatives at k={ks} x 2.000004 MHz, score {r['score']:.1f}",
        ks == [-135, -63, 40] and r["score"] >= SCORE_VETO)
    chk("relatives: same clock family in ANOTHER session (k=-100) and in a NON-overlapping scan "
        "of the same session (k=+77) do not count", -100 not in ks and 77 not in ks)
    r0 = db.relatives(982.0026, 0.03, "Parkes:58602", 1800.0, 3.81)
    chk("...without mjd the non-overlapping same-session scan does count (k=+77), the other session never",
        77 in [k for _, k, *_ in r0["found"]] and -100 not in [k for _, k, *_ in r0["found"]])
    rn = db.relatives(982.0026, 0.03, "Parkes:58602", 1800.0, 3.81, mjd=58602.6, neigh_hz=5e3)
    chk(f"relatives chance rate is measured locally: a +-5 kHz window holding only the match costs "
        f"2 tol/10 kHz each [lambda {rn['lam']:.3f}, score {rn['score']:.1f} < 3]",
        0.2 < rn["lam"] < 0.35 and rn["score"] < SCORE_VETO)
    rc = db.relatives(982.0026, 0.03, "Parkes:58602", 1800.0, 3.81, mjd=58602.6, offset=0.5,
                      ks=[-135, -63, 40])
    chk(f"half-clock control positions find nothing [n={rc['n_k']}]", rc["n_k"] == 0)
    r = db.relatives(982.0026 + 0.37, 0.03, "Parkes:58602", 1800.0, 3.81)
    chk(f"control frequency 370 kHz off blc1 has no significant relatives [n={r['n_k']}, "
        f"score {r['score']:.1f}]", r["score"] < SCORE_VETO)
    # drift-track SIGN and the reported-frequency convention (turboSETI reports f at t=0):
    # DB hit track [F, F+150 Hz]; slack at 2 GHz ~ 72 Hz < the 100 Hz gaps below
    F = 2000.0
    for off, want, why in ((100, 1, "query starting inside the DB track"),
                           (-100, 1, "query whose track ends inside the DB track"),
                           (250, 0, "query starting 100 Hz after the DB track ends"),
                           (-250, 0, "query ending 100 Hz before the DB track starts")):
        v = db.query(F + off * 1e-6, 0.5, target="STARX")
        chk(f"track sign: +0.5 Hz/s {why} ({off:+d} Hz) -> n_t {v['n_targets']} (want {want})",
            v["n_targets"] == want)
    v = db.query(2010.0, 0.0, target="STARX")
    chk(f"ONE other-target hit in an empty neighbourhood does NOT veto (rate floor, lambda "
        f"{v['lam']:.1e} > 0, score {v['score']:.1f}, background {v['n_background']})",
        not v["veto"] and v["lam"] > 0 and v["score"] < 10)
    # loader reference: turboSETI .dat columns (Drift_Rate = col 1, Uncorrected_Frequency = col 3)
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        dp = Path(td) / "x.dat"
        dp.write_text("# Source:GJ699\n# MJD: 57457.555625\tRA: 17h57m47s\tDEC: 4d44m21s\n"
                      "# DELTAT:  17.986224\tDELTAF(Hz):  -2.835503\tmax_drift_rate:   4.0\t"
                      "obs_length: 287.779586\n# --------------------------\n"
                      "# Top_Hit_#\tDrift_Rate\tSNR\tUncorrected_Frequency\tCorrected_Frequency\tIndex\n"
                      "001\t -0.251305\t 30.5\t 1420.123456\t 1420.123999\t 123\t 1420.1\t 1420.2\n")
        rr = rows_turboseti_dat(dp)
        chk("turboSETI .dat loader: Uncorrected (topocentric) frequency, drift, res, T, session",
            len(rr) == 1 and rr[0]["freq_mhz"] == 1420.123456 and rr[0]["drift_hz_s"] == -0.251305
            and abs(rr[0]["res_hz"] - 2.835503) < 1e-9 and abs(rr[0]["tobs_s"] - 287.779586) < 1e-9
            and rr[0]["session"] == "GBT:57457" and rr[0]["drift_valid"] == 1)
    # coverage: drift-keyed veto removes far less than a frequency-only one would
    rng = np.random.default_rng(5)
    probes = rng.uniform(1100, 1900, 400)
    keyed = np.mean([db.query(x, 0.5, target="P")["n_targets"] > 0 for x in probes])
    global DRIFT_ABS
    keep, DRIFT_ABS = DRIFT_ABS, 1e3
    try:
        fonly = np.mean([db.query(x, 0.5, target="P")["n_targets"] > 0 for x in probes])
    finally:
        DRIFT_ABS = keep
    chk(f"drift keying matches fewer random probes than frequency-only ({100 * keyed:.1f} % vs "
        f"{100 * fonly:.1f} %)", keyed < fonly)
    vb = np.mean([db.query(x, 0.5, target="P")["veto"] for x in probes])
    chk(f"random probes vetoed (false-veto proxy) <= 2 % [{100 * vb:.1f} %]", vb <= 0.02)
    # round trip + seed file
    with tempfile.TemporaryDirectory() as td:
        for ext in (".csv.gz", ".npz"):
            p = Path(td) / ("t" + ext)
            db.save(p)
            d2 = RFIDB.load(p)
            chk(f"save/load round trip {ext}", d2.n == db.n and np.allclose(d2.c["freq_mhz"], db.c["freq_mhz"])
                and d2.vocab["target"] == db.vocab["target"])
    if SEED.exists():
        s = RFIDB.load(SEED)
        srcs = {x.split(" ")[0] for x in s.vocab["source"]}
        chk(f"seed file loads: {s.n:,} hits, {len(s.fam_vocab)} target families", s.n > 1000)
        chk("seed has our own GJ699 hits and Price 2020 events",
            any("GJ699" in x for x in s.vocab["source"]) and "Price2020" in srcs)
    else:
        chk("seed file present", False)
    chk(f"RFI_BANDS merged unchanged ({len(allocations())} entries, "
        f"{100 * band_fraction():.0f} % of 1.1-1.9 GHz)", len(allocations()) >= 14)
    chk("family(): GJ699_OFF/ProxCen_S/1421-490_ND_S", family("GJ699_OFF") == "gj699"
        and family("ProxCen_S") == "proxcen" and family("1421-490_ND_S") == "1421-490")
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd")
    sub.add_parser("selftest")
    p = sub.add_parser("info"); p.add_argument("db", nargs="?")
    p = sub.add_parser("query")
    p.add_argument("freq", type=float); p.add_argument("drift", type=float)
    p.add_argument("--mjd", type=float); p.add_argument("--target")
    p.add_argument("--tobs", type=float, default=300.0); p.add_argument("--res", type=float, default=2.8)
    p.add_argument("--session", help="also look for intermod relatives in this session")
    p.add_argument("--db")
    p = sub.add_parser("build"); p.add_argument("out"); p.add_argument("sources", nargs="*")
    p.add_argument("--enriquez-dir")
    p = sub.add_parser("explain"); p.add_argument("events"); p.add_argument("--db")
    p.add_argument("--exclude-session", action="append"); p.add_argument("-v", "--verbose", action="store_true")
    p = sub.add_parser("coverage"); p.add_argument("--db")
    p.add_argument("--n", type=int, default=2000); p.add_argument("--seed", type=int, default=1)
    p.add_argument("--lo", type=float, default=1100.0); p.add_argument("--hi", type=float, default=1900.0)
    p.add_argument("--drifts", type=lambda s: [float(x) for x in s.split(",")],
                   default=[0.0, -0.12, 0.5, -1.0, 2.0])
    p = sub.add_parser("blc1"); p.add_argument("--db")
    a = ap.parse_args(argv)
    if a.cmd in (None, "selftest"):
        return selftest()
    if a.cmd == "build":
        return cmd_build(a)
    db = default_db(getattr(a, "db", None))
    return {"info": lambda: cmd_info(db), "query": lambda: cmd_query(db, a),
            "explain": lambda: cmd_explain(db, a), "coverage": lambda: cmd_coverage(db, a),
            "blc1": lambda: cmd_blc1(db, a)}[a.cmd]()


def cmd_build(a):
    """Sources: *.dat (turboSETI), *events.csv (Price 2020), rows CSV(.gz) in schema,
    *.npz (saved DB). --enriquez-dir adds the reduced Enriquez 2017 chunks."""
    dbs = []
    for s in a.sources:
        p = Path(s)
        if p.suffix == ".dat":
            dbs.append(RFIDB.from_rows(rows_turboseti_dat(p, source=f"seti-cookbook GJ699 {p.name}")))
        elif p.name.endswith("_events.csv"):
            tel = "Parkes" if "pks" in p.name else "GBT"
            dbs.append(RFIDB.from_rows(rows_price2020(p, telescope=tel,
                                                      res=3.0 if tel == "Parkes" else 2.79)))
        else:
            d = RFIDB.load(p)
            dbs.append(d)
        print(f"  {p.name}: {dbs[-1].n:,} hits")
    if a.enriquez_dir:
        dbs.append(enriquez_db(a.enriquez_dir))
        print(f"  Enriquez 2017: {dbs[-1].n:,} hits")
    db = RFIDB.concat(dbs)
    db.save(a.out)
    print(f"wrote {a.out}: {db.n:,} hits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
