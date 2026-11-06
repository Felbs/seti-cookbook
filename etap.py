#!/usr/bin/env python3
"""etap.py - the UCLA injection-recovery protocol, for any narrowband search.

eta_P ("signal recovery rate") is the fraction of planted signals a search
gets back when the plants are drawn the way an unknown transmitter would be:
start frequency uniform over the band, drift uniform over +-8.86 Hz/s, and a
strength set so the search ITSELF should report SNR ~20. The UCLA SETI group
(Margot et al. 2023; Choza et al. 2024) measured it for their own pipeline and
for a Breakthrough-Listen-like turboSETI configuration. This file
measures it for ours, on real Green Bank data, and turns it into an upper
bound on the fraction of stars with a transmitter (reach.py --etap).

What is different from UCLA, and labelled so everywhere:
  * UCLA planted into raw VOLTAGES and channelised; we plant into the already
    channelised fine-resolution spectra (setigen, Doppler-smeared, inject.py).
    Our eta_P therefore leaves out losses from scalloping and PFB leakage
    between fine channels that voltage injection includes.
  * Strength is reported twice: in the search's own SNR units (what UCLA sets
    to ~20) AND in TRUE radiometer SNR (inject.radiometer_noise), which is what
    reach.py and the radiometer equation mean by SNR (EVIDENCE.md rule 4).

Search backends (plug-ins, all return hits as dict(ch, drift, snr) with ch =
the hit's channel at the start of the scan, in the data's own channel order):
  turboseti  - turboSETI 2.3.2 exactly as R9 ran it (max_drift 4, min_drift 0,
               snr 10, n_coarse_chan 1) on a one-coarse-channel HDF5 slice.
  reference  - a small brute-force de-Doppler in numpy (this file): shift and
               sum every drift on a grid, optional box sum over the per-spectrum
               smear, one hit per suppression window. Lets the selftest run
               without turboSETI.
  callable   - any function fn(data, meta) -> hits, e.g. a future fastsearch;
               use CallableBackend(fn, ...).

Matching, reported both ways:
  UCLA       - one-to-one; a hit is credited if its start frequency is within
               +-6 Hz and its drift within +-0.05 Hz/s of the plant's.
  path       - inject.match (one-to-one along the whole track) and, stricter,
               inject.explain_misses' "found" (path AND drift_ok).

    python etap.py selftest                       # synthetic, no downloads, no turboSETI
    python etap.py run FILE.npz --seeds 1,2,3 --json out.json
    python etap.py h5 FULL.h5 --coarse 10,20,30 --json out.json   # read coarse channels
    python etap.py summary out1.json out2.json ... [--json etap_summary.json]
"""
import argparse
import json
import math
import os
import shutil
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import inject  # noqa: E402

UCLA_DRIFT = 8.86            # Hz/s, Margot 2023 / Choza 2024 drift range
UCLA_TOL_HZ = 6.0            # start-frequency tolerance
UCLA_TOL_DRIFT = 0.05        # Hz/s
GBT_NOTCH = (1200.0, 1340.0)  # MHz, the L-band receiver's radar notch
DRIFT_SUBSETS = (8.86, 4.0, 2.0)
UCLA_BL_LIKE = 0.057         # quoted BL-like control, +-8.86 Hz/s (not re-derived here)


# ------------------------------------------------------------------ statistics
def binom_ci(k, n, cl=0.95):
    """Exact (Clopper-Pearson) two-sided interval for k successes in n trials."""
    from scipy.stats import beta
    if n == 0:
        return (0.0, 1.0)
    a = 1 - cl
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return (lo, hi)


# ------------------------------------------------------------------ backends
class Backend:
    """A search. Subclasses implement search(); snr_scale() says how many TRUE
    radiometer SNR units one unit of THIS search's SNR is for an ideal
    (zero-drift, channel-centred) tone, so plants can be set to 'search SNR 20'."""
    name = "backend"
    max_drift = 4.0          # Hz/s the search is asked for
    snr = 10.0               # its threshold, in its own units

    def search(self, data, meta):
        raise NotImplementedError

    def snr_scale(self, data, meta):
        return 1.0

    def suppression_ch(self, meta):
        """Half-width (channels) inside which this search keeps only the
        strongest hit (turboSETI: 2*ceil(max_drift / drift resolution))."""
        T = meta["nt"] * meta["tsamp"]
        res = abs(meta["foff_hz"]) / T
        return 2 * int(math.ceil(self.max_drift / res))


class TurboSETIBackend(Backend):
    """turboSETI 2.3.2 with the R9 settings, on a one-coarse-channel slice."""
    name = "turboseti"

    def __init__(self, max_drift=4.0, snr=10.0, min_drift=0.0, workdir=None):
        self.max_drift, self.snr, self.min_drift = max_drift, snr, min_drift
        self.workdir = workdir

    def snr_scale(self, data, meta):
        # turboSETI divides by the 5-95 % trimmed std of the time-summed coarse
        # channel (KNOWN_ISSUES #1); the radiometer noise excludes the bandpass.
        return float(inject.turboseti_noise(data) / inject.radiometer_noise(data))

    def search(self, data, meta):
        import recovery
        from turbo_seti.find_doppler.find_doppler import FindDoppler
        out = Path(tempfile.mkdtemp(prefix="etap_ts_", dir=self.workdir))
        try:
            fn = out / "slice.h5"
            recovery.write_h5(fn, data, dict(fch1=meta["fch1_mhz"], foff=meta["foff_hz"] / 1e6,
                                             tsamp=meta["tsamp"], tstart=meta.get("tstart", 57457.0),
                                             ra_h=meta.get("ra_h", 0.0), dec_deg=meta.get("dec_deg", 0.0)))
            cwd = os.getcwd()
            os.chdir(out)
            try:
                import contextlib
                import io
                with contextlib.redirect_stdout(io.StringIO()):     # turboSETI prints headers
                    FindDoppler(fn.name, max_drift=self.max_drift, min_drift=self.min_drift,
                                snr=self.snr, out_dir=".", n_coarse_chan=1, append_output=False,
                                log_level_int=40).search()
            finally:
                os.chdir(cwd)
            return inject.read_dat(fn.with_suffix(".dat"), meta["fch1_mhz"], meta["foff_hz"] / 1e6)
        finally:
            shutil.rmtree(out, ignore_errors=True)


def dedoppler_brute(data, foff_hz, tsamp, max_drift, snr, drift_step=None, boxcar=True,
                    radius=None):
    """Reference de-Doppler search, numpy only.

    Bandpass out (each channel over the 65-channel running median of the mean
    spectrum), then for every drift on a grid (default step: one channel over
    the scan) sum the track through all spectra. With boxcar=True each
    spectrum contributes a box of 1+round(|r|) channels, r = channels moved per
    spectrum, so a fast tone's smeared power is collected (turboSETI sums one
    channel). SNR = (sum - median) / robust std of that drift's row, i.e. TRUE
    radiometer units for a zero-drift channel-centred tone. One hit per
    +-radius channels (default turboSETI's own window). Returns hits."""
    from scipy.ndimage import median_filter, maximum_filter1d
    nt, nc = data.shape
    x = np.asarray(data, dtype=np.float64)
    bp = median_filter(x.mean(0), size=65)
    bp[bp <= 0] = np.median(bp[bp > 0]) if np.any(bp > 0) else 1.0
    x = x / bp - 1.0
    c = np.zeros((nt, nc + 1))
    np.cumsum(x, axis=1, out=c[:, 1:])
    del x
    step = drift_step or abs(foff_hz) / ((nt - 1) * tsamp)
    nd = int(math.floor(max_drift / step))
    drifts = np.arange(-nd, nd + 1) * step
    best = np.full(nc, -np.inf)
    bestk = np.zeros(nc, dtype=np.int32)
    sub = slice(None, None, 7)
    i = np.arange(nt)
    for k, d in enumerate(drifts):
        r = d * tsamp / foff_hz                     # channels per spectrum, signed
        w = 1 + int(round(abs(r))) if boxcar else 1
        off = np.round(r * (i + 0.5) - (w - 1) / 2.0).astype(int)
        lo, hi = max(0, -off.min()), nc - (off.max() + w)
        if hi <= lo:
            continue
        row = np.zeros(hi - lo)
        for t in range(nt):
            a = lo + off[t]
            row += c[t, a + w:a + w + (hi - lo)] - c[t, a:a + (hi - lo)]
        med = np.median(row[sub])
        sd = inject.robust_std(row[sub])
        if sd <= 0:
            continue
        s = (row - med) / sd
        seg = best[lo:hi]
        m = s > seg
        seg[m] = s[m]
        bestk[lo:hi][m] = k
    if radius is None:
        T = nt * tsamp
        radius = 2 * int(math.ceil(max_drift / (abs(foff_hz) / T)))
    peak = (best >= snr) & (best >= maximum_filter1d(best, size=2 * radius + 1, mode="constant",
                                                     cval=-np.inf))
    return [dict(ch=float(j), drift=float(drifts[bestk[j]]), snr=float(best[j]))
            for j in np.nonzero(peak)[0]]


class ReferenceBackend(Backend):
    """Brute-force numpy de-Doppler (dedoppler_brute)."""
    name = "reference"

    def __init__(self, max_drift=4.0, snr=10.0, boxcar=True, drift_step=None):
        self.max_drift, self.snr, self.boxcar, self.drift_step = max_drift, snr, boxcar, drift_step

    def search(self, data, meta):
        return dedoppler_brute(data, meta["foff_hz"], meta["tsamp"], self.max_drift, self.snr,
                               drift_step=self.drift_step, boxcar=self.boxcar,
                               radius=self.suppression_ch(meta))


class CallableBackend(Backend):
    """Wrap any search function fn(data, meta) -> [dict(ch, drift, snr), ...].
    snr_scale: TRUE SNR per unit of fn's SNR; None = measure it (calibrate())."""

    def __init__(self, fn, name="callable", max_drift=4.0, snr=10.0, snr_scale=None):
        self.fn, self.name, self.max_drift, self.snr = fn, name, max_drift, snr
        self._scale = snr_scale

    def search(self, data, meta):
        return list(self.fn(data, meta))

    def snr_scale(self, data, meta):
        if self._scale is None:
            self._scale = calibrate_scale(self, data, meta)
        return self._scale


def calibrate_scale(backend, data, meta, n=12, true_snr=100.0):
    """Measure TRUE SNR per unit of the backend's SNR: plant n channel-centred
    zero-drift tones at true_snr into a copy and take the median ratio."""
    nc = data.shape[1]
    chs = np.linspace(0.1 * nc, 0.9 * nc, n).round()
    plants = [dict(ch=float(c), drift=0.0, snr=true_snr) for c in chs]
    d = data.copy()
    inject.plant(d, meta["foff_hz"], meta["tsamp"], plants)
    hits = backend.search(d, meta)
    r = [true_snr / h["snr"] for p in plants for h in hits if abs(h["ch"] - p["ch"]) <= 2]
    return float(np.median(r)) if r else float("nan")


# ------------------------------------------------------------------ the protocol
def allowed_ranges_mhz(lo_mhz, hi_mhz, exclude=(GBT_NOTCH,)):
    """[lo, hi] minus excluded MHz ranges -> list of (lo, hi)."""
    out = [(lo_mhz, hi_mhz)]
    for a, b in exclude:
        nxt = []
        for x, y in out:
            if b <= x or a >= y:
                nxt.append((x, y))
                continue
            if a > x:
                nxt.append((x, a))
            if b < y:
                nxt.append((b, y))
        out = nxt
    return out


def design_plants(meta, per_mhz, drift_max, gap_ch, rng, exclude=(GBT_NOTCH,), edge_ch=64,
                  max_tries=200000):
    """UCLA-style plants: start frequency uniform over the slice (minus excluded
    ranges), drift uniform in +-drift_max, n = per_mhz * usable MHz; tracks kept
    gap_ch apart (beyond the search's suppression window) by rejection.
    Returns (plants, n_wanted)."""
    nt, nc, foff_hz, ts = meta["nt"], meta["nchan"], meta["foff_hz"], meta["tsamp"]
    f_a = meta["fch1_mhz"]
    f_b = f_a + (nc - 1) * foff_hz / 1e6
    ranges = allowed_ranges_mhz(min(f_a, f_b), max(f_a, f_b), exclude)
    widths = np.array([b - a for a, b in ranges])
    usable = float(widths.sum())
    n_want = int(round(per_mhz * usable))
    if n_want == 0:
        return [], 0
    T = nt * ts
    taken = []                # (lo_ch, hi_ch) of placed tracks
    plants = []
    tries = 0
    while len(plants) < n_want and tries < max_tries:
        tries += 1
        r = rng.choice(len(ranges), p=widths / usable)
        f = rng.uniform(*ranges[r])
        dr = rng.uniform(-drift_max, drift_max)
        ch0 = (f - f_a) * 1e6 / foff_hz
        ch1 = ch0 + dr * T / foff_hz
        lo, hi = min(ch0, ch1), max(ch0, ch1)
        if lo < edge_ch or hi > nc - 1 - edge_ch:
            continue
        if any(lo < b + gap_ch and hi > a - gap_ch for a, b in taken):
            continue
        taken.append((lo, hi))
        plants.append(dict(ch=float(ch0), drift=float(dr), f_mhz=float(f)))
    return plants, n_want


def match_ucla(plants, hits, foff_hz, tol_hz=UCLA_TOL_HZ, tol_drift=UCLA_TOL_DRIFT):
    """One-to-one: start frequency within tol_hz AND drift within tol_drift."""
    from scipy.optimize import linear_sum_assignment
    out = [None] * len(plants)
    if not plants or not hits:
        return out
    pc = np.array([p["ch"] for p in plants])
    pd = np.array([p["drift"] for p in plants])
    hc = np.array([h["ch"] for h in hits])
    hd = np.array([h["drift"] for h in hits])
    dfz = np.abs(pc[:, None] - hc[None, :]) * abs(foff_hz)
    ddr = np.abs(pd[:, None] - hd[None, :])
    ok = (dfz <= tol_hz) & (ddr <= tol_drift)
    if not ok.any():
        return out
    rows = np.nonzero(ok.any(1))[0]
    cols = np.nonzero(ok.any(0))[0]
    cost = np.where(ok, dfz / tol_hz + ddr / tol_drift, 1e9)[np.ix_(rows, cols)]
    r, c = linear_sum_assignment(cost)
    for i, j in zip(r, c):
        if cost[i, j] < 1e9:
            out[int(rows[i])] = int(cols[j])
    return out


def in_rfi(f_mhz):
    import recipe_api as R
    return R.rfi_reason(f_mhz) is not None


def trial(data, meta, backend, seed, per_mhz=75.0, target_snr=20.0, drift_max=UCLA_DRIFT,
          exclude=(GBT_NOTCH,), null_hits=None, verbose=True):
    """One injection-recovery pass on one slice. Returns (records, info)."""
    rng = np.random.default_rng(seed)
    nt, nc = data.shape
    foff_hz, ts = meta["foff_hz"], meta["tsamp"]
    scale = backend.snr_scale(data, meta)
    true_snr = target_snr * scale
    gap = backend.suppression_ch(meta) + 8
    plants, n_want = design_plants(meta, per_mhz, drift_max, gap, rng, exclude)
    for p in plants:
        p["snr"] = true_snr
    if not plants:
        return [], dict(n_want=n_want, n=0)
    t0 = time.time()
    d = data.copy()
    inject.plant(d, foff_hz, ts, plants, seed=seed)
    t_plant = time.time() - t0
    t0 = time.time()
    hits = backend.search(d, meta)
    t_search = time.time() - t0
    del d
    m_u = match_ucla(plants, hits, foff_hz)
    m_p = inject.match(plants, hits, foff_hz, ts, nt)
    why = inject.explain_misses(plants, hits, m_p, foff_hz, ts, nt, backend.max_drift)
    own = null_hits or []
    recs = []
    T = (nt - 1) * ts
    for i, p in enumerate(plants):
        f_mid = p["f_mhz"] + p["drift"] * T / 2 / 1e6
        hu = hits[m_u[i]] if m_u[i] is not None else None
        hp = hits[m_p[i]] if m_p[i] is not None else None
        near_own = any(abs(o["ch"] - p["ch"]) <= abs(p["drift"]) * T / abs(foff_hz) + gap for o in own)
        recs.append(dict(
            f_mhz=round(p["f_mhz"], 7), drift=round(p["drift"], 5), true_snr=round(true_snr, 2),
            search_snr=target_snr, rfi=bool(in_rfi(f_mid)), near_own=bool(near_own),
            ucla=hu is not None, path=hp is not None, path_drift=why[i]["status"] == "found",
            status=why[i]["status"],
            hit_snr=None if hp is None else round(hp["snr"], 2),
            hit_dhz=None if hp is None else round((hp["ch"] - p["ch"]) * foff_hz, 2),
            hit_ddrift=None if hp is None else round(hp["drift"] - p["drift"], 4)))
    info = dict(n_want=n_want, n=len(plants), hits=len(hits), scale=scale, true_snr=true_snr,
                t_plant_s=round(t_plant, 1), t_search_s=round(t_search, 1),
                max_abs_hit_drift=max([abs(h["drift"]) for h in hits], default=0.0))
    if verbose:
        k = sum(r["ucla"] for r in recs)
        print(f"  seed {seed}: {len(plants)}/{n_want} plants at search SNR {target_snr:g} "
              f"(= true SNR {true_snr:.1f}, scale {scale:.2f}); {len(hits)} hits; "
              f"UCLA-matched {k} ({k / len(plants):.3f}); plant {t_plant:.0f}s search {t_search:.0f}s",
              flush=True)
    return recs, info


def summarize(recs, cl=0.95, quiet=False):
    """eta_P with exact binomial CIs: drift subsets x (all / outside RFI) x match method."""
    out = {}
    for dmax in DRIFT_SUBSETS:
        for where in ("all", "outside_rfi", "outside_rfi_and_own_lines"):
            sel = [r for r in recs if abs(r["drift"]) <= dmax
                   and (where == "all" or not r["rfi"])
                   and (where != "outside_rfi_and_own_lines" or not r["near_own"])]
            n = len(sel)
            row = dict(n=n)
            for meth in ("ucla", "path", "path_drift"):
                k = sum(1 for r in sel if r[meth])
                lo, hi = binom_ci(k, n, cl)
                row[meth] = dict(k=k, eta=(k / n if n else float("nan")), lo=lo, hi=hi)
            out[f"{dmax:g}|{where}"] = row
    # by |drift| bin, UCLA match
    bins = [0, 0.15, 0.5, 1, 2, 3, 3.78, 4, 6, 8.86]
    by = []
    for a, b in zip(bins[:-1], bins[1:]):
        sel = [r for r in recs if a <= abs(r["drift"]) < b]
        k = sum(r["ucla"] for r in sel)
        kp = sum(r["path_drift"] for r in sel)
        by.append(dict(lo=a, hi=b, n=len(sel), k_ucla=k, k_path_drift=kp))
    out["by_drift"] = by
    st = {}
    for r in recs:
        st[r["status"]] = st.get(r["status"], 0) + 1
    out["statuses"] = st
    if not quiet:
        print_summary(out, cl)
    return out


def print_summary(s, cl=0.95):
    pct = int(round(cl * 100))
    print(f"\neta_P (fraction recovered), exact {pct}% binomial intervals")
    print(f"{'|drift|<=':>9s} {'where':26s} {'n':>6s}  {'UCLA match':>22s}  {'path match':>22s}  {'path+drift':>22s}")
    for dmax in DRIFT_SUBSETS:
        for where in ("all", "outside_rfi", "outside_rfi_and_own_lines"):
            row = s[f"{dmax:g}|{where}"]
            cells = "  ".join(f"{row[m]['eta']:.3f} [{row[m]['lo']:.3f},{row[m]['hi']:.3f}]".rjust(22)
                              for m in ("ucla", "path", "path_drift"))
            print(f"{dmax:9.2f} {where:26s} {row['n']:6d}  {cells}")
    print("\nby |drift| (Hz/s):   n   UCLA   path+drift")
    for b in s["by_drift"]:
        if b["n"]:
            print(f"  {b['lo']:5.2f}-{b['hi']:5.2f}  {b['n']:5d}  {b['k_ucla'] / b['n']:.3f}  "
                  f"{b['k_path_drift'] / b['n']:.3f}")
    print("miss/find reasons (path match):", s["statuses"])


def band_fractions(lo_mhz, hi_mhz, exclude=(GBT_NOTCH,)):
    """D for reach.f_tx_upper: fraction of [lo, hi] searched (minus notch), and
    fraction outside notch AND outside recipe_api.RFI_BANDS."""
    import recipe_api as R
    tot = hi_mhz - lo_mhz
    usable = allowed_ranges_mhz(lo_mhz, hi_mhz, exclude)
    clean = allowed_ranges_mhz(lo_mhz, hi_mhz, list(exclude) + [(a, b) for a, b, _ in R.RFI_BANDS])
    return (sum(b - a for a, b in usable) / tot, sum(b - a for a, b in clean) / tot)


# ------------------------------------------------------------------ data access
def meta_from_npz(path):
    import recovery
    data, m = recovery.load(path)
    meta = dict(fch1_mhz=float(m["fch1"]), foff_hz=float(m["foff"]) * 1e6, tsamp=float(m["tsamp"]),
                tstart=float(m["tstart"]), ra_h=float(m["ra_h"]), dec_deg=float(m["dec_deg"]),
                nt=data.shape[0], nchan=data.shape[1], source=Path(path).name)
    return data, meta


def read_coarse(h5path, k, nfine=1048576):
    """Coarse channel k (0-based, nfine channels) of a full BL fine-res HDF5."""
    import h5py
    import hdf5plugin  # noqa: F401
    with h5py.File(h5path, "r") as f:
        ds = f["data"]
        a = dict(ds.attrs)
        data = ds[:, 0, k * nfine:(k + 1) * nfine].astype(np.float32)
    meta = dict(fch1_mhz=float(a["fch1"]) + float(a["foff"]) * k * nfine, foff_hz=float(a["foff"]) * 1e6,
                tsamp=float(a["tsamp"]), tstart=float(a["tstart"]), ra_h=float(a["src_raj"]),
                dec_deg=float(a["src_dej"]), nt=data.shape[0], nchan=data.shape[1],
                source=f"{Path(h5path).name}#coarse{k}", coarse=k)
    return data, meta


def make_backend(name, max_drift, snr):
    if os.environ.get("SETI_COOKBOOK_GPU", "0") in ("1", "true", "yes"):
        print("note: SETI_COOKBOOK_GPU=1 set, but etap backends run on CPU (GPU optional, not needed)")
    if name == "turboseti":
        return TurboSETIBackend(max_drift=max_drift, snr=snr)
    if name == "reference":
        return ReferenceBackend(max_drift=max_drift, snr=snr)
    if name == "reference-1ch":
        return ReferenceBackend(max_drift=max_drift, snr=snr, boxcar=False)
    raise SystemExit(f"unknown backend {name}")


def run_slices(slices, a):
    """slices: iterable of (data, meta). Runs null probe + seeds on each."""
    be = make_backend(a.backend, a.max_drift, a.snr)
    allrecs, infos = [], []
    for data, meta in slices:
        print(f"{meta['source']}: {meta['nt']} x {meta['nchan']} ch from {meta['fch1_mhz']:.6f} MHz, "
              f"{abs(meta['foff_hz']):.4f} Hz x {meta['tsamp']:.2f} s; backend {be.name} "
              f"max_drift {be.max_drift:g} snr {be.snr:g}", flush=True)
        t0 = time.time()
        own = be.search(data, meta)
        print(f"  null probe: {len(own)} hits of the data's own ({time.time() - t0:.0f}s)", flush=True)
        for seed0 in a.seed_list:
            # a different draw in every slice: the same seed everywhere would plant the SAME
            # drift list in every coarse channel (correlated samples, CI too narrow)
            seed = seed0 if (a.fixed_seed or "coarse" not in meta) else seed0 + 7919 * (meta["coarse"] + 1)
            recs, info = trial(data, meta, be, seed, per_mhz=a.per_mhz, target_snr=a.target_snr,
                               drift_max=a.drift_max, null_hits=own,
                               exclude=() if a.no_notch else (GBT_NOTCH,))
            for r in recs:
                r["source"], r["seed"] = meta["source"], seed
            info.update(source=meta["source"], seed=seed, null_hits=len(own))
            allrecs += recs
            infos.append(info)
            if a.json:   # write as we go: a long run that dies keeps what it measured
                Path(a.json).write_text(json.dumps(dict(protocol=protocol_dict(a, be), infos=infos,
                                                        records=allrecs)))
        del data
    s = summarize(allrecs)
    if a.json:
        Path(a.json).write_text(json.dumps(dict(protocol=protocol_dict(a, be), infos=infos,
                                                records=allrecs, summary=s)))
        print(f"wrote {a.json}")
    return 0


def protocol_dict(a, be):
    return dict(backend=be.name, max_drift=be.max_drift, snr=be.snr, per_mhz=a.per_mhz,
                target_search_snr=a.target_snr, drift_max=a.drift_max,
                tol_hz=UCLA_TOL_HZ, tol_drift=UCLA_TOL_DRIFT, notch_excluded=not a.no_notch,
                injection="post-channelisation (setigen, Doppler-smeared), NOT voltage")


def cmd_summary(paths, out_json=None, band=(1126.4, 1876.5), n_stars=1, cl=0.95):
    recs, infos, prot = [], [], None
    for p in paths:
        z = json.loads(Path(p).read_text())
        recs += z["records"]
        infos += z["infos"]
        prot = prot or z["protocol"]
    print(f"{len(paths)} file(s), {len(infos)} trials, {len(recs)} plants; protocol {prot}")
    scales = [i["scale"] for i in infos]
    print(f"search-SNR -> true-SNR scale: median {np.median(scales):.2f} "
          f"(range {min(scales):.2f}-{max(scales):.2f}); true SNR of plants "
          f"{np.median([i['true_snr'] for i in infos]):.1f} median")
    print(f"largest |drift| reported by the search: {max(i['max_abs_hit_drift'] for i in infos):.4f} Hz/s")
    s = summarize(recs)
    f_rfi = np.mean([r["rfi"] for r in recs])
    print(f"plants inside RFI_BANDS: {f_rfi:.3f}")
    d_all, d_out = band_fractions(*band)
    import reach
    key = f"{UCLA_DRIFT:g}"
    e_all = s[f"{key}|all"]["ucla"]["eta"]
    e_out = s[f"{key}|outside_rfi"]["ucla"]["eta"]
    res = dict(n_plants=len(recs), protocol=prot, summary=s, band_mhz=list(band), D_all=d_all,
               D_outside_rfi=d_out, eta_all=e_all, eta_outside_rfi=e_out,
               eta_all_ci=[s[f"{key}|all"]["ucla"]["lo"], s[f"{key}|all"]["ucla"]["hi"]],
               eta_outside_rfi_ci=[s[f"{key}|outside_rfi"]["ucla"]["lo"], s[f"{key}|outside_rfi"]["ucla"]["hi"]])
    print(f"\nD (band {band[0]:g}-{band[1]:g} MHz minus the 1200-1340 notch) = {d_all:.3f}; "
          f"outside RFI_BANDS too = {d_out:.3f}")
    reach.print_ftx(res, n_stars=n_stars, cl=cl)
    if out_json:
        Path(out_json).write_text(json.dumps(res, indent=1))
        print(f"wrote {out_json}")
    return 0


# ------------------------------------------------------------------ selftest
def selftest():
    checks = []

    def chk(name, ok, extra=""):
        checks.append((name, bool(ok)))
        print(f"  {'PASS' if ok else 'FAIL'}  {name}{('  ' + extra) if extra else ''}", flush=True)

    # 1. binomial CI against hand / textbook values
    lo, hi = binom_ci(0, 10)
    chk("CI k=0 n=10 upper = 1-0.025^(1/10) = 0.3085", abs(hi - (1 - 0.025 ** 0.1)) < 1e-4 and lo == 0)
    lo, hi = binom_ci(5, 10)
    chk("CI k=5 n=10 = [0.1871, 0.8129] (Clopper-Pearson)", abs(lo - 0.1871) < 1e-3 and abs(hi - 0.8129) < 1e-3)
    lo, hi = binom_ci(57, 1000)
    chk("CI k=57 n=1000 brackets 0.057 and is ~+-1.5 %", lo < 0.057 < hi and 0.04 < lo and hi < 0.075,
        f"[{lo:.4f},{hi:.4f}]")

    # 2. f_TX formula against hand-computed values
    import reach
    v1 = reach.f_tx_upper(1, 0.5, 0.8, 0.95)      # 0.95 / 0.4
    v2 = reach.f_tx_upper(100, 0.5, 0.8, 0.95)    # (1 - 0.05^0.01) / 0.4 = 0.029513 / 0.4
    chk("f_TX n=1 eta=0.5 D=0.8 CL=95% = 2.375 (no constraint, >1)", abs(v1 - 2.375) < 1e-9)
    chk("f_TX n=100 eta=0.5 D=0.8 CL=95% = 0.07378", abs(v2 - 0.073782) < 2e-6, f"{v2:.6f}")

    # 3. band fractions: notch only, 1100-1900
    d_all, d_out = band_fractions(1100, 1900)
    chk("D: 1100-1900 MHz minus 140 MHz notch = 0.825", abs(d_all - 0.825) < 1e-9)
    chk("D outside RFI < D all", d_out < d_all)

    # 4. the reference backend on synthetic noise
    rng = np.random.default_rng(7)
    nt, nc, foff, ts = 16, 131072, -2.835503418452676, 17.986224128
    noise = rng.chisquare(64, (nt, nc)).astype(np.float32)
    bp = (1 + 0.3 * np.cos(np.linspace(0, 2 * np.pi, nc))).astype(np.float32)   # a bandpass
    noise *= bp
    meta = dict(fch1_mhz=1500.0, foff_hz=foff, tsamp=ts, nt=nt, nchan=nc, source="synthetic")
    be = ReferenceBackend(max_drift=2.2, snr=10.0, drift_step=0.02)
    t0 = time.time()
    own = be.search(noise, meta)
    chk("control: pure noise gives <= 1 false hit at SNR 10", len(own) <= 1,
        f"{len(own)} hits, {time.time() - t0:.1f}s")
    sc = calibrate_scale(be, noise, meta, n=8)
    chk("reference backend SNR is TRUE radiometer SNR (scale 0.85-1.15)", 0.85 < sc < 1.15, f"{sc:.3f}")
    be_c = CallableBackend(lambda d, m: dedoppler_brute(d, m["foff_hz"], m["tsamp"], 2.2, 10.0,
                                                        drift_step=0.02),
                           name="callable(reference)", max_drift=2.2)
    sc_c = be_c.snr_scale(noise, meta)
    chk("generic callable backend calibrates itself the same way", abs(sc_c - sc) < 0.05, f"{sc_c:.3f}")

    hi_r, _ = trial(noise, meta, be, seed=1, per_mhz=100, target_snr=60, drift_max=2.0,
                    exclude=(), verbose=False)
    eta_hi = np.mean([r["ucla"] for r in hi_r])
    chk("eta_P ~1 at search SNR 60 (|drift| <= 2)", eta_hi >= 0.9, f"eta {eta_hi:.3f} of {len(hi_r)}")
    lo_r, _ = trial(noise, meta, be, seed=2, per_mhz=100, target_snr=3, drift_max=2.0,
                    exclude=(), verbose=False)
    eta_lo = np.mean([r["ucla"] for r in lo_r])
    chk("eta_P ~0 at search SNR 3", eta_lo <= 0.06, f"eta {eta_lo:.3f} of {len(lo_r)}")

    # 5. matching controls: a hit with the opposite drift, or 10 Hz away, must NOT match
    plants = [dict(ch=1000.0, drift=1.0), dict(ch=5000.0, drift=-0.5), dict(ch=9000.0, drift=0.3)]
    hits = [dict(ch=1000.5, drift=1.02, snr=20), dict(ch=5000.0, drift=+0.5, snr=20),
            dict(ch=9000.0 + 10 / abs(foff), drift=0.3, snr=20)]
    m = match_ucla(plants, hits, foff)
    chk("UCLA match: right hit credited", m[0] == 0)
    chk("UCLA match: opposite drift NOT credited", m[1] is None)
    chk("UCLA match: 10 Hz off NOT credited", m[2] is None)
    m2 = match_ucla([dict(ch=1000.0, drift=1.0), dict(ch=1000.3, drift=1.01)], [hits[0]], foff)
    chk("UCLA match: one hit never credited to two plants", sum(x is not None for x in m2) == 1)

    # 6. plant design: density, drift range, spacing beyond suppression window, notch respected
    gap = be.suppression_ch(meta) + 8
    pl, nw = design_plants(dict(meta, nchan=1048576), 75, UCLA_DRIFT, gap, np.random.default_rng(3))
    T = nt * ts
    iv = sorted((min(p["ch"], p["ch"] + p["drift"] * T / foff), max(p["ch"], p["ch"] + p["drift"] * T / foff))
                for p in pl)
    sep = min(b[0] - a[1] for a, b in zip(iv[:-1], iv[1:]))
    chk("design: 75 plants/MHz placed, all kept apart beyond the window",
        len(pl) == nw and sep >= gap, f"{len(pl)}/{nw}, min gap {sep:.0f} ch >= {gap}")
    chk("design: drifts span +-8.86 Hz/s", max(abs(p["drift"]) for p in pl) <= UCLA_DRIFT and
        max(abs(p["drift"]) for p in pl) > 8.0)
    pn, nwn = design_plants(dict(meta, fch1_mhz=1338.0, nchan=1048576), 75, 2.0, gap, np.random.default_rng(4))
    chk("design: nothing planted inside the 1200-1340 MHz notch",
        all(not (GBT_NOTCH[0] <= p["f_mhz"] <= GBT_NOTCH[1]) for p in pn) and nwn < 0.75 * 75 * 2.97)

    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL", f"({sum(c[1] for c in checks)}/{len(checks)})")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    for nm in ("run", "h5"):
        p = sub.add_parser(nm)
        p.add_argument("path")
        if nm == "h5":
            p.add_argument("--coarse", required=True, help="comma list of coarse-channel indices")
        p.add_argument("--backend", default="turboseti", choices=["turboseti", "reference", "reference-1ch"])
        p.add_argument("--max-drift", type=float, default=4.0)
        p.add_argument("--snr", type=float, default=10.0)
        p.add_argument("--seeds", default="1")
        p.add_argument("--per-mhz", type=float, default=75.0)
        p.add_argument("--target-snr", type=float, default=20.0, help="in the SEARCH's own SNR units")
        p.add_argument("--drift-max", type=float, default=UCLA_DRIFT)
        p.add_argument("--no-notch", action="store_true")
        p.add_argument("--fixed-seed", action="store_true",
                       help="use each seed as given in every slice (cross-machine reproduction)")
        p.add_argument("--json")
    p = sub.add_parser("summary")
    p.add_argument("paths", nargs="+")
    p.add_argument("--json")
    p.add_argument("--band", default="1126.4,1876.5", help="MHz lo,hi of the searched file")
    p.add_argument("--n-stars", type=int, default=1)
    p.add_argument("--cl", type=float, default=0.95)
    a = ap.parse_args()
    if a.cmd == "selftest":
        return selftest()
    if a.cmd == "summary":
        lo, hi = (float(v) for v in a.band.split(","))
        return cmd_summary(a.paths, a.json, (lo, hi), a.n_stars, a.cl)
    a.seed_list = [int(s) for s in a.seeds.split(",")]
    if a.cmd == "run":
        return run_slices([meta_from_npz(a.path)], a)
    ks = [int(k) for k in a.coarse.split(",")]
    return run_slices((read_coarse(a.path, k) for k in ks), a)


if __name__ == "__main__":
    sys.exit(main())
