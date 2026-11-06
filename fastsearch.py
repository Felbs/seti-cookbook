#!/usr/bin/env python3
"""fastsearch.py - single-read narrowband drift search, in TRUE radiometer SNR.

turboSETI re-reads every coarse channel many times (KNOWN_ISSUES.md #4), never
searches drift exactly 0 (#2), stops short of the drift asked for (#3) and
quotes an SNR ~4x smaller than the radiometer SNR (#1). This search:

  1. reads each HDF5 file ONCE, one coarse channel at a time (h5py hyperslab
     data[:, 0, c0:c1]), never the whole file, with a MemAvailable guard;
  2. flattens the bandpass per coarse channel with a 1024-channel BLOCK MEDIAN
     of the time-mean spectrum across FREQUENCY (never a per-channel time
     median, which erases steady carriers - law), divides each spectrum by it
     and by its own gain, and blanks the DC spike at the coarse-channel centre;
  3. de-Dopplers by brute force over EVERY drift k * step, k = -K..K, including
     drift exactly 0; step = |foff| / ((nt-1) * tsamp) (one channel over the
     scan, 0.0105 Hz/s on GBT fine-res), K = ceil(max_drift / step) so the
     whole asked range is covered; each spectrum is shifted by round(k t/(nt-1)).
     Where the smear is >= 16 channels the rows are computed on data summed
     over s = 2^j <= w/8 channels with drift step s*step (drift_plan; --widths
     full turns this off);
  4. sums each drift track over a BOXCAR of w = max(1, round(|drift| tsamp/df))
     channels, the width a drifting tone is smeared over inside one spectrum,
     so fast drifters keep their power (measured gain: see FASTSEARCH_RESULT.md);
  5. quotes SNR in true radiometer units: the zero-drift row's noise is the
     channel-to-channel noise of the bandpass-removed time-summed spectrum,
     exactly what inject.radiometer_noise measures; boxcar rows are divided by
     the measured noise of the same boxcar on that spectrum;
  6. keeps a hit only if it is the strongest whose TRACK comes within
     (half-widths + 2) channels of it at any time in the scan (track-crossing
     suppression; turboSETI instead drops everything within a fixed
     +-max_drift*T/df channels of a stronger hit);
  7. writes a turboSETI-compatible .dat (find_event_pipeline reads it) and a
     CSV with detection SNR, power SNR, drift, freq at t=0, coarse channel and
     the local noise around each hit.

numpy by default; numba is used if installed (same results, faster); cupy only
with SETI_COOKBOOK_GPU=1 (never required, not verified on real hardware here).

    python fastsearch.py selftest
    python fastsearch.py FILE.h5 --max-drift 4 --snr 8 --out-dir out/ [--workers 3]
    python fastsearch.py FILE.h5 --coarse 150-160        # a subset of coarse channels
"""
import argparse
import csv
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

BLOCK = 1024          # bandpass block-median width (channels)
DC_HALF = 1           # blank DC channel +-1
SUPPRESS_TOL = 2.0    # extra channels in the track-crossing suppression
MAX_HITS_PER_COARSE = 3000
SCRUNCH_AT = 8      # scrunch by s = 2^j only while s <= w / SCRUNCH_AT (w = smear in channels)
SCRUNCH_PHASES = 1  # scrunch phases (offsets of s/phases channels)
SCRUNCH_EXTRA = 0   # extra scrunched channels in the boxcar of a scrunched row
GBT_COARSE_MHZ = 2.9296875


# ---------------------------------------------------------------- utilities
def robust_std(x):
    x = np.asarray(x, dtype=np.float64)
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))


def mem_available_bytes():
    try:
        for ln in open("/proc/meminfo"):
            if ln.startswith("MemAvailable:"):
                return int(ln.split()[1]) * 1024
    except OSError:
        pass
    return None


def wait_for_memory(need, reserve=1 << 30, timeout=600):
    """Block until MemAvailable >= need + reserve (law: MemAvailable, not MemFree)."""
    t0 = time.time()
    while True:
        av = mem_available_bytes()
        if av is None or av >= need + reserve:
            return av
        if time.time() - t0 > timeout:
            raise MemoryError(f"MemAvailable {av/2**30:.2f} GB < need {need/2**30:.2f} GB + reserve")
        time.sleep(5)


def backend():
    """'cupy' only on explicit opt-in, else 'numba' if importable, else 'numpy'."""
    if os.environ.get("SETI_COOKBOOK_GPU") == "1":
        try:
            import cupy  # noqa: F401
            return "cupy"
        except Exception:
            print("SETI_COOKBOOK_GPU=1 but cupy is not importable - using the CPU", file=sys.stderr)
    if os.environ.get("FASTSEARCH_NO_NUMBA") != "1":
        try:
            import numba  # noqa: F401
            return "numba"
        except Exception:
            pass
    return "numpy"


# ---------------------------------------------------------------- per-slab steps
def normalise(raw, dc_idx=()):
    """raw (nt, nf) -> x (nt, nf) float32 with x = raw / (bandpass * gain_t) - 1,
    DC channels blanked to 0. Also returns the time-summed spectrum S of x."""
    raw = np.asarray(raw, dtype=np.float32)
    nt, nf = raw.shape
    m = raw.mean(0, dtype=np.float64)
    for c in dc_idx:                                  # keep the spike out of the bandpass fit
        lo, hi = max(0, c - DC_HALF), min(nf, c + DC_HALF + 1)
        m[lo:hi] = np.nan
    nb = max(1, nf // BLOCK)
    edges = np.linspace(0, nf, nb + 1).astype(int)
    cen = 0.5 * (edges[:-1] + edges[1:])
    med = np.array([np.nanmedian(m[a:b]) for a, b in zip(edges[:-1], edges[1:])])
    good = np.isfinite(med) & (med > 0)
    bp = np.interp(np.arange(nf), cen[good], med[good]).astype(np.float32)
    x = raw / bp
    gain = np.median(x[:, ::7], axis=1).astype(np.float32)   # per-spectrum gain
    x /= gain[:, None]
    x -= 1.0
    for c in dc_idx:
        x[:, max(0, c - DC_HALF):min(nf, c + DC_HALF + 1)] = 0.0
    return x, x.sum(0, dtype=np.float64)


def boxcar(a, w, sign):
    """sum of a[c .. c+w-1] (sign>0) or a[c-w+1 .. c] (sign<0); zero outside."""
    if w <= 1:
        return a
    cs = np.concatenate(([0.0], np.cumsum(a, dtype=np.float64)))
    n = a.size
    c = np.arange(n)
    if sign >= 0:
        hi = np.minimum(c + w, n)
        return (cs[hi] - cs[c]).astype(np.float32)
    lo = np.maximum(c - w + 1, 0)
    return (cs[c + 1] - cs[lo]).astype(np.float32)


def drift_plan(K, nt, mode="auto"):
    """Which drifts are searched on which data. Returns a list of
    (s, ks, ws): scrunch factor s (channels summed in frequency first), drift
    indices ks (in full-resolution units, drift = k * step), and boxcar widths
    ws in scrunched channels.

    mode "auto":  boxcar w = round(|k|/(nt-1)) channels (the per-spectrum smear);
                  once w >= 16 the search runs on data pre-summed over s = 2^j <= w/8
                  channels, drift step s (one scrunched channel over the scan).
                  Measured cost vs "full": 0 % at <= 1 Hz/s, 6-15 % at 1.5-3.5 Hz/s
                  on GBT fine-res, for ~3x less CPU.
    mode "full":  same boxcar, every drift at full resolution (exact, ~5x slower).
    mode "1":     no boxcar, full resolution (what a one-channel-per-spectrum
                  tree sums; kept to MEASURE the boxcar gain)."""
    groups = {}
    for k in range(-K, K + 1):
        u = abs(k) / (nt - 1)
        w = 1 if mode == "1" else max(1, int(round(u)))
        sc = 1
        if mode == "auto":
            while sc * SCRUNCH_AT <= w:
                sc *= 2
        if sc > 1 and k % sc and abs(k) != K:
            continue
        g = groups.setdefault(sc, ([], []))
        g[0].append(k)
        g[1].append((max(1, int(round(u / sc))) + (SCRUNCH_EXTRA if sc > 1 else 0)) if mode != "1" else 1)
    return [(sc, np.array(ks, dtype=np.int64), np.array(ws, dtype=np.int64))
            for sc, (ks, ws) in sorted(groups.items())]


def noise_for_widths(S, ws):
    """Robust noise of the time-summed spectrum S after a boxcar of each width,
    with the local trend removed (65-channel running median, as
    inject.radiometer_noise). Width 1 on unscrunched data is the TRUE radiometer
    noise in these units; drift rows sum the same independent spectra, so a
    boxcar on S has the same noise as the boxcar on any drift row."""
    from scipy.ndimage import median_filter
    sub = (S - median_filter(S, size=65, mode="nearest")).astype(np.float32)
    return {int(w): robust_std(boxcar(sub, int(w), +1)[::3]) for w in sorted(set(int(v) for v in ws))}


def _shift(k, t, nt, sc):
    return int(np.floor(k * t / ((nt - 1) * sc) + 0.5))


def _dd_numpy(x, ks, sc, widths, sig, xp=np):
    nt, nf = x.shape
    P = int(np.max(np.abs(ks))) // sc + 2
    xpad = xp.zeros((nt, nf + 2 * P), dtype=xp.float32)
    xpad[:, P:P + nf] = xp.asarray(x)
    best = xp.full(nf, -np.inf, dtype=xp.float32)
    bestk = xp.zeros(nf, dtype=xp.int32)
    D = xp.empty(nf, dtype=xp.float32)
    c = xp.arange(nf)
    for i, k in enumerate(ks):
        k = int(k)
        D[:] = 0
        for t in range(nt):
            o = P + _shift(k, t, nt, sc)
            D += xpad[t, o:o + nf]
        w = int(widths[i])
        if w > 1:
            cs = xp.concatenate((xp.zeros(1, dtype=xp.float64), xp.cumsum(D, dtype=xp.float64)))
            if k >= 0:
                B = cs[xp.minimum(c + w, nf)] - cs[c]
            else:
                B = cs[c + 1] - cs[xp.maximum(c - w + 1, 0)]
            B = (B / float(sig[i])).astype(xp.float32)
        else:
            B = D / float(sig[i])
        upd = B > best
        best = xp.where(upd, B, best)
        bestk = xp.where(upd, k, bestk)
    if xp is not np:
        best, bestk = xp.asnumpy(best), xp.asnumpy(bestk)
    return best, bestk


_NUMBA_KERNEL = None


def _numba_kernel():
    global _NUMBA_KERNEL
    if _NUMBA_KERNEL is None:
        from numba import njit

        @njit(cache=True)
        def kern(xpad, P, nf, ks, sc, widths, inv_sig, best, bestk, D):
            nt = xpad.shape[0]
            for i in range(ks.size):
                k = ks[i]
                for c in range(nf):
                    D[c] = 0.0
                for t in range(nt):
                    o = P + int(np.floor(k * t / ((nt - 1) * sc) + 0.5))
                    row = xpad[t]
                    for c in range(nf):
                        D[c] += row[o + c]
                w = widths[i]
                g = inv_sig[i]
                if w <= 1:
                    for c in range(nf):
                        v = D[c] * g
                        if v > best[c]:
                            best[c] = v
                            bestk[c] = k
                elif k >= 0:
                    acc = 0.0
                    for c in range(min(w, nf)):
                        acc += D[c]
                    for c in range(nf):          # window D[c .. c+w-1]
                        v = acc * g
                        if v > best[c]:
                            best[c] = v
                            bestk[c] = k
                        acc -= D[c]
                        if c + w < nf:
                            acc += D[c + w]
                else:
                    acc = 0.0
                    for c in range(nf):          # window D[c-w+1 .. c]
                        acc += D[c]
                        if c - w >= 0:
                            acc -= D[c - w]
                        v = acc * g
                        if v > best[c]:
                            best[c] = v
                            bestk[c] = k
        _NUMBA_KERNEL = kern
    return _NUMBA_KERNEL


def _dd_numba(x, ks, sc, widths, sig):
    nt, nf = x.shape
    P = int(np.max(np.abs(ks))) // sc + 2
    xpad = np.zeros((nt, nf + 2 * P), dtype=np.float32)
    xpad[:, P:P + nf] = x
    best = np.full(nf, -np.inf, dtype=np.float32)
    bestk = np.zeros(nf, dtype=np.int32)
    D = np.empty(nf, dtype=np.float32)
    _numba_kernel()(xpad, P, nf, np.asarray(ks, dtype=np.int64), sc, np.asarray(widths, dtype=np.int64),
                    (1.0 / np.asarray(sig, dtype=np.float64)), best, bestk, D)
    return best, bestk


def dedoppler(x, ks, sc, widths, sig, how=None):
    how = how or backend()
    if how == "numba":
        return _dd_numba(x, ks, sc, widths, sig)
    if how == "cupy":
        import cupy as cp
        return _dd_numpy(x, ks, sc, widths, np.asarray(sig), xp=cp)
    return _dd_numpy(x, ks, sc, widths, np.asarray(sig))


def select_hits(best, bestk, bestw, thr, nt, local_win=3, cap=MAX_HITS_PER_COARSE):
    """Local maxima above thr, then greedy track-crossing suppression: a weaker
    candidate is dropped if its track comes within (half-width_a + half-width_b
    + SUPPRESS_TOL) channels of an already kept stronger hit's track at the first
    or last spectrum, or crosses it in between. Returns list of (chan, k, snr) and
    the number of channels above thr."""
    above = np.flatnonzero(best >= thr)
    if above.size == 0:
        return [], 0
    from scipy.ndimage import maximum_filter1d
    lm = maximum_filter1d(best, size=2 * local_win + 1, mode="constant", cval=-np.inf)
    cand = above[best[above] >= lm[above]]
    n_above = int(above.size)
    order = cand[np.argsort(-best[cand], kind="stable")]
    if order.size > 50 * cap:
        order = order[:50 * cap]
    c = order.astype(np.float64)
    u = bestk[order] / (nt - 1)                       # channels per spectrum
    hw = np.maximum(1, bestw[order]) / 2.0
    # track centre in spectrum t: c + u (t + 0.5); compare at t = 0 and t = nt - 1
    a0 = c + 0.5 * u
    a1 = c + (nt - 0.5) * u
    alive = np.ones(order.size, bool)
    kept = []
    for i in range(order.size):
        if not alive[i]:
            continue
        kept.append((int(order[i]), int(bestk[order[i]]), float(best[order[i]])))
        if len(kept) >= cap:
            break
        d0 = a0 - a0[i]
        d1 = a1 - a1[i]
        cross = np.sign(d0) != np.sign(d1)
        dmin = np.where(cross, 0.0, np.minimum(np.abs(d0), np.abs(d1)))
        alive &= ~(dmin <= hw + hw[i] + SUPPRESS_TOL)
    return kept, n_above


def local_noise(S, chans, half=4096):
    """Robust channel-to-channel noise of the time-summed spectrum around each hit
    (+-half channels, the hit's own +-64 channels excluded)."""
    from scipy.ndimage import median_filter
    sub = S - median_filter(S, size=65, mode="nearest")
    out = []
    for c in chans:
        lo, hi = max(0, c - half), min(S.size, c + half)
        seg = np.concatenate((sub[lo:max(lo, c - 64)], sub[min(hi, c + 65):hi]))
        out.append(robust_std(seg) if seg.size > 100 else np.nan)
    return out


def search_slab(raw, foff_hz, tsamp, max_drift=4.0, thr=8.0, dc_idx=(), widths="auto", how=None):
    """Search one coarse channel (or any contiguous block) of shape (nt, nf).
    Returns (hits, info). hit: dict(ch = start channel at t=0, k, drift Hz/s,
    snr = detection SNR (true radiometer units), snr_power = boxcar power over the
    zero-drift radiometer noise (estimates the tone's true SNR), w = boxcar
    width in channels, s = scrunch, local_noise, noise_ratio)."""
    nt, nf = raw.shape
    step = abs(foff_hz) / ((nt - 1) * tsamp)
    K = int(np.ceil(max_drift / step - 1e-9))
    x, S = normalise(raw, dc_idx)
    plan = drift_plan(K, nt, widths)
    best = np.full(nf, -np.inf, dtype=np.float32)
    bestk = np.zeros(nf, dtype=np.int32)
    bestw = np.ones(nf, dtype=np.int32)          # boxcar width, full-resolution channels
    bests = np.ones(nf, dtype=np.int32)          # scrunch factor
    bestp = np.zeros(nf, dtype=np.int32)         # scrunch phase
    sig1 = None
    noise = {}
    n_drifts = 0
    for sc, ks, ws in plan:
        for ph in ([0] if sc == 1 else list(range(0, sc, max(1, sc // SCRUNCH_PHASES)))):
            if sc == 1:
                xs, Ss = x, S
                m = nf
            else:
                m = ((nf - ph) // sc) * sc
                xs = x[:, ph:ph + m].reshape(nt, m // sc, sc).sum(2, dtype=np.float32)
                Ss = xs.sum(0, dtype=np.float64)
            if (sc, ph) == (sc, 0) or sc == 1:
                nz = noise_for_widths(Ss, ws)
                if sc == 1:
                    sig1 = nz.get(1) or noise_for_widths(S, [1])[1]
                for w, v in nz.items():
                    noise[(sc, w)] = v
                sig = np.array([nz[int(w)] for w in ws])
            b, bk = dedoppler(xs, ks, sc, ws, sig, how)
            n_drifts += len(ks) / sc
            if sc > 1:
                b = np.repeat(b, sc)
                bk = np.repeat(bk, sc)
            wlut = np.ones(2 * K + 1, dtype=np.int32)
            wlut[ks + K] = ws * sc
            sl = slice(ph, ph + m)
            upd = b > best[sl]
            best[sl] = np.where(upd, b, best[sl])
            bestk[sl] = np.where(upd, bk, bestk[sl])
            bestw[sl] = np.where(upd, wlut[bk + K], bestw[sl])
            bests[sl] = np.where(upd, sc, bests[sl])
            bestp[sl] = np.where(upd, ph, bestp[sl])
            del xs
    if sig1 is None:
        sig1 = noise_for_widths(S, [1])[1]
    kept, n_above = select_hits(best, bestk, bestw, thr, nt)
    # a scrunched hit's start is known to its s-channel block: report the block centre
    rep = [((c - bestp[c]) // bests[c]) * bests[c] + bestp[c] + (bests[c] - 1) / 2.0 for c, _, _ in kept]
    ln = local_noise(S, [int(r) for r in rep])
    hits = []
    for (c, k, snr), ch, lnz in zip(kept, rep, ln):
        sc = int(bests[c])
        w = int(bestw[c])
        drift = (k / (nt - 1)) * foff_hz / tsamp       # Hz/s; foff<0 -> k>0 is a negative drift
        hits.append(dict(ch=float(ch), k=k, drift=float(drift), snr=float(snr),
                         snr_power=float(snr * noise[(sc, w // sc)] / sig1), w=w, s=sc,
                         local_noise=float(lnz), noise_ratio=float(lnz / sig1)))
    info = dict(K=K, step=step, drift_searched=K * step, sigma1=sig1, n_above=n_above,
                n_drifts_equiv=round(n_drifts, 1), plan=[(sc, len(ks)) for sc, ks, _ in plan])
    return hits, info


# ---------------------------------------------------------------- file driver
def header_of(f):
    a = dict(f["data"].attrs)
    h = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in a.items()}
    return h


def coarse_width(h, n_coarse=None):
    """Fine channels per coarse channel, the way blimpy decides (calc_n_coarse_chan)."""
    nchans = int(h["nchans"])
    if n_coarse:
        n = int(n_coarse)
    elif nchans >= 2 ** 20 and nchans % 2 ** 20 == 0:
        n = nchans // 2 ** 20
    elif int(h.get("telescope_id", 0)) == 6:
        n = int(round(abs(float(h["foff"])) * nchans / GBT_COARSE_MHZ))
    else:
        n = max(1, nchans // 2 ** 20)
    if nchans % n:
        raise SystemExit(f"nchans {nchans} not divisible by {n} coarse channels; pass --n-coarse")
    return nchans // n, n


def _deg_str(v, hours):
    from astropy.coordinates import Angle
    import astropy.units as u
    a = Angle(float(v), unit=u.hourangle if hours else u.deg)
    return a.to_string(unit=u.hourangle if hours else u.deg, sep=("h", "m", "s") if hours else ("d", "m", "s"),
                       precision=3)


def write_dat(path, rows, h, file_id, max_drift_searched, obs_length):
    with open(path, "w") as fo:
        fo.write("# -------------------------- o --------------------------\n")
        fo.write(f"# File ID: {file_id} \n")
        fo.write("# -------------------------- o --------------------------\n")
        fo.write(f"# Source:{h.get('source_name', 'unknown')}\n")
        fo.write(f"# MJD: {float(h['tstart']):18.12f}\tRA: {_deg_str(h['src_raj'], True)}\t"
                 f"DEC: {_deg_str(h['src_dej'], False)}\n")
        fo.write(f"# DELTAT: {float(h['tsamp']):10.6f}\tDELTAF(Hz): {float(h['foff'])*1e6:10.6f}\t"
                 f"max_drift_rate: {max_drift_searched:10.6f}\tobs_length: {obs_length:10.6f}\n")
        fo.write("# --------------------------\n")
        fo.write("# Top_Hit_# \tDrift_Rate \tSNR \tUncorrected_Frequency \tCorrected_Frequency \tIndex \t"
                 "freq_start \tfreq_end \tSEFD \tSEFD_freq \tCoarse_Channel_Number \tFull_number_of_hits \t\n")
        fo.write("# --------------------------\n")
        for i, r in enumerate(rows, 1):
            fo.write(f"{i:06d}\t{r['drift']:10.6f}\t{r['snr']:10.6f}\t{r['freq']:14.6f}\t{r['freq']:14.6f}\t"
                     f"{r['index']:d}\t{r['freq']:14.6f}\t{r['freq_end']:14.6f}\t0.0\t{0.0:14.6f}\t"
                     f"{r['coarse']:d}\t{r['n_above']:d}\t\n")


CSV_COLS = ["coarse_chan", "index", "global_chan", "freq_mhz_t0", "freq_mhz_end", "drift_hz_s", "snr",
            "snr_power", "boxcar_w", "local_noise", "noise_ratio", "n_above_in_coarse"]


def _worker(args):
    path, coarse_list, nfpc, opts = args
    import h5py
    try:
        import hdf5plugin  # noqa: F401
    except ImportError:
        pass
    out = []
    with h5py.File(path, "r") as f:
        h = header_of(f)
        ds = f["data"]
        nt = ds.shape[0]
        foff_hz = float(h["foff"]) * 1e6
        for cc in coarse_list:
            c0 = cc * nfpc
            wait_for_memory(nt * nfpc * 4 * 6)
            t0 = time.time()
            raw = ds[:, 0, c0:c0 + nfpc]
            t_read = time.time() - t0
            hits, info = search_slab(raw, foff_hz, float(h["tsamp"]), opts["max_drift"], opts["snr"],
                                     dc_idx=(nfpc // 2,), widths=opts["widths"])
            del raw
            out.append((cc, hits, info, t_read, time.time() - t0))
            if opts.get("verbose"):
                print(f"  coarse {cc:4d}: {len(hits):5d} hits, read {t_read:5.1f} s, total "
                      f"{time.time()-t0:5.1f} s", flush=True)
    return out


def run_file(path, out_dir, max_drift=4.0, snr=8.0, workers=1, coarse=None, n_coarse=None,
             widths="auto", verbose=True):
    import h5py
    try:
        import hdf5plugin  # noqa: F401
    except ImportError:
        pass
    path = Path(path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "r") as f:
        h = header_of(f)
        nt = f["data"].shape[0]
    nfpc, ncoarse = coarse_width(h, n_coarse)
    todo = list(range(ncoarse)) if coarse is None else [c for c in coarse if 0 <= c < ncoarse]
    per_slab = nt * nfpc * 4 * 6
    av = mem_available_bytes()
    if av is not None:
        workers = max(1, min(workers, int((av - (1 << 30)) // per_slab)))
    t0 = time.time()
    opts = dict(max_drift=max_drift, snr=snr, widths=widths, verbose=verbose)
    print(f"{path.name}: {nt} spectra x {int(h['nchans'])} ch = {ncoarse} coarse x {nfpc}; "
          f"{len(todo)} to search with {workers} worker(s), backend {backend()}", flush=True)
    # contiguous runs per worker (each byte read once)
    chunks = [todo[i::1] for i in range(1)]
    if workers > 1:
        n = len(todo)
        bounds = np.linspace(0, n, workers * 4 + 1).astype(int)
        chunks = [todo[a:b] for a, b in zip(bounds[:-1], bounds[1:]) if b > a]
    jobs = [(str(path), ch, nfpc, opts) for ch in chunks]
    results = []
    if workers > 1:
        from multiprocessing import get_context
        with get_context("spawn").Pool(workers) as pool:
            for r in pool.imap_unordered(_worker, jobs):
                results.extend(r)
    else:
        for j in jobs:
            results.extend(_worker(j))
    results.sort(key=lambda r: r[0])
    foff = float(h["foff"])
    fch1 = float(h["fch1"])
    tsamp = float(h["tsamp"])
    obs_length = nt * tsamp
    rows = []
    K = step = None
    t_read = sum(r[3] for r in results)
    for cc, hits, info, _, _ in results:
        K, step = info["K"], info["step"]
        for hh in hits:
            g = cc * nfpc + int(hh["ch"])
            fr = fch1 + foff * g
            rows.append(dict(coarse=cc, index=int(hh["ch"]), global_chan=g, freq=fr,
                             freq_end=fr + hh["drift"] * obs_length / 1e6, drift=hh["drift"],
                             snr=hh["snr"], snr_power=hh["snr_power"], w=hh["w"],
                             local_noise=hh["local_noise"], noise_ratio=hh["noise_ratio"],
                             n_above=info["n_above"]))
    rows.sort(key=lambda r: -r["freq"] if foff < 0 else r["freq"])
    searched = K * step if K is not None else 0.0
    stem = path.name[:-3] if path.name.endswith(".h5") else path.stem
    write_dat(out_dir / f"{stem}.dat", rows, h, path.name, searched, obs_length)
    with open(out_dir / f"{stem}.fastsearch.csv", "w", newline="") as fo:
        wr = csv.writer(fo)
        wr.writerow(CSV_COLS)
        for r in rows:
            wr.writerow([r["coarse"], r["index"], r["global_chan"], f"{r['freq']:.6f}", f"{r['freq_end']:.6f}",
                         f"{r['drift']:.6f}", f"{r['snr']:.3f}", f"{r['snr_power']:.3f}", r["w"],
                         f"{r['local_noise']:.6g}", f"{r['noise_ratio']:.3f}", r["n_above"]])
    wall = time.time() - t0
    print(f"done: {len(rows)} hits in {len(results)} coarse channels, wall {wall/60:.2f} min "
          f"(read {t_read/60:.2f} worker-min); drift searched +-{searched:.4f} Hz/s in steps of "
          f"{step:.6f} (incl. 0); threshold SNR {snr} (true radiometer units)", flush=True)
    return rows, dict(wall_s=wall, n_coarse=len(results), drift_searched=searched, step=step)


def parse_coarse(s):
    if not s:
        return None
    out = []
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


# ---------------------------------------------------------------- selftest
def _synthetic(rng, nt=16, nf=131072, dc=True):
    """chi-square noise (~100 dof like GBT fine-res) x a smooth scalloped bandpass
    (0.6 at the coarse-channel edges, as measured on GJ699 2016-03-10, plus a
    ripple), per-spectrum gain drift, optional DC spike at the centre."""
    f = np.linspace(-1, 1, nf)
    bp = (1.0 - 0.4 * f ** 2) * (1 + 0.05 * np.sin(9 * f))
    gain = 1 + 0.03 * np.sin(np.arange(nt))
    d = (rng.chisquare(100, (nt, nf)) / 100).astype(np.float32)
    d *= (bp[None, :] * gain[:, None]).astype(np.float32) * 1e9
    if dc:
        d[:, nf // 2] *= 25.0
    return d


def selftest():
    import inject
    foff, ts = -2.835503418452676, 17.986224128
    nt, nf = 16, 131072
    checks = []
    how = backend()
    print(f"backend: {how}")

    # 1. control: pure noise (+ bandpass + gain) and NO DC spike -> false-alarm count
    rng = np.random.default_rng(7)
    noise = _synthetic(rng, nt, nf, dc=False)
    t0 = time.time()
    h0, info = search_slab(noise, foff, ts, 4.0, 8.0, dc_idx=(nf // 2,))
    dt = time.time() - t0
    print(f"  noise only: {len(h0)} hits >= 8, max {max([h['snr'] for h in h0], default=0):.2f}; "
          f"{info['K']*2+1} drifts, +-{info['drift_searched']:.4f} Hz/s step {info['step']:.6f}; {dt:.1f} s")
    checks.append(("pure noise: 0 hits at SNR 8", len(h0) == 0))
    checks.append(("drift range covers +-4 Hz/s and includes 0", info["drift_searched"] >= 4.0))
    h6, _ = search_slab(noise, foff, ts, 4.0, 6.5, dc_idx=(nf // 2,))
    print(f"  noise only at SNR 6.5: {len(h6)} hits (Gaussian tail, ~763 correlated drifts per channel)")
    checks.append(("pure noise at SNR 6.5: few hits (<= 10 per 131k channels)", len(h6) <= 10))
    # radiometer-unit noise: on a FLAT bandpass, sigma1 must equal inject.radiometer_noise
    flat = (np.random.default_rng(12).chisquare(100, (nt, nf)) / 100).astype(np.float32)
    _, inf2 = search_slab(flat, foff, ts, 0.5, 8.0)
    rel = inject.radiometer_noise(flat) / np.median(flat.sum(0))
    print(f"  flat noise: sigma1 {inf2['sigma1']/nt:.5f} (units of the mean) vs inject.radiometer_noise "
          f"{rel:.5f}; radiometer eq. 1/sqrt(16*50) = {1/np.sqrt(16*50):.5f}")
    checks.append(("noise unit matches inject.radiometer_noise within 5 %",
                   abs(inf2["sigma1"] / nt / rel - 1) < 0.05))

    # 2. DC spike control: spike present, no plants -> no hit at the DC channel
    rng = np.random.default_rng(8)
    spk = _synthetic(rng, nt, nf, dc=True)
    hd, _ = search_slab(spk, foff, ts, 4.0, 8.0, dc_idx=(nf // 2,))
    near_dc = [h for h in hd if abs(h["ch"] - nf // 2) < 30]
    print(f"  DC spike x25: {len(hd)} hits, {len(near_dc)} near the DC channel")
    checks.append(("DC spike does not become a hit", len(near_dc) == 0))

    # 3. smeared plants: found, correct drift, SNR in true units
    rng = np.random.default_rng(9)
    d = _synthetic(rng, nt, nf, dc=True)
    base = d.copy()
    plants = [dict(ch=9000.3, drift=0.0, snr=40), dict(ch=21000.0, drift=0.1, snr=40),
              dict(ch=33000.6, drift=-0.3, snr=40), dict(ch=45000.2, drift=1.0, snr=80),
              dict(ch=57000.0, drift=-1.0, snr=80), dict(ch=80000.5, drift=2.0, snr=120),
              dict(ch=95000.0, drift=-3.5, snr=180), dict(ch=110000.4, drift=3.5, snr=180)]
    inject.plant(d, foff, ts, plants)
    hp, info = search_slab(d, foff, ts, 4.0, 8.0, dc_idx=(nf // 2,))
    m = inject.match(plants, hp, foff, ts, nt)
    why = inject.explain_misses(plants, hp, m, foff, ts, nt, 4.0)
    T = (nt - 1) * ts
    for p, j, wy in zip(plants, m, why):
        h = hp[j] if j is not None else None
        print(f"  plant drift {p['drift']:+5.2f} true SNR {p['snr']:4d}: {wy['status']:16s}"
              + (f" drift {h['drift']:+.4f} snr {h['snr']:6.1f} power-snr {h['snr_power']:6.1f} w {h['w']}"
                 if h else ""))
    checks.append(("every plant found with the correct drift",
                   all(w["status"] == "found" for w in why)))
    extra = len(hp) - sum(j is not None for j in m)
    near_dc = [h for h in hp if abs(h["ch"] - nf // 2) < 30]
    print(f"  extra (unplanted) hits: {extra}")
    checks.append(("no extra hits beside the plants", extra == 0))
    checks.append(("DC spike not a hit when plants are present", len(near_dc) == 0))
    # SNR calibration on the zero-drift plant: zero-drift channel-centred tone -> ~planted SNR
    rng2 = np.random.default_rng(10)
    d2 = _synthetic(rng2, nt, nf, dc=False)
    pz = [dict(ch=float(c), drift=0.0, snr=50) for c in range(5000, 125000, 6000)]
    inject.plant(d2, foff, ts, pz)
    hz, _ = search_slab(d2, foff, ts, 4.0, 8.0)
    mz = inject.match(pz, hz, foff, ts, nt)
    got = [hz[j]["snr"] for j in mz if j is not None]
    ratio = np.median(got) / 50 if got else 0
    print(f"  zero-drift true SNR 50 x{len(pz)}: reported median {np.median(got):.1f} (ratio {ratio:.3f}; "
          f"bandpass 0.55-1.3 x makes per-plant true SNR vary)")
    checks.append(("zero-drift SNR in true units (median within 25 % of planted)", 0.75 < ratio < 1.25))
    # boxcar gain at 3.5 Hz/s
    hb, _ = search_slab(d, foff, ts, 4.0, 5.0, dc_idx=(nf // 2,), widths="1")
    mb = inject.match(plants, hb, foff, ts, nt)
    for i in (6, 7):
        s_box = hp[m[i]]["snr"] if m[i] is not None else 0
        s_one = hb[mb[i]]["snr"] if mb[i] is not None else 0
        print(f"  drift {plants[i]['drift']:+.1f}: SNR with boxcar {s_box:.1f}, without {s_one:.1f}")
    s_box = np.mean([hp[m[i]]["snr"] for i in (6, 7) if m[i] is not None] or [0])
    s_one = np.mean([hb[mb[i]]["snr"] for i in (6, 7) if mb[i] is not None] or [1e-9])
    checks.append(("boxcar raises SNR of a 3.5 Hz/s smeared tone by > 1.5x", s_box > 1.5 * s_one))

    # 4. numpy and numba agree (if numba is here), full-res and scrunched rows
    if how == "numba":
        x, S = normalise(d[:, :16384])
        for sc, ks, ws in drift_plan(381, nt, "auto")[::2]:
            xs = x if sc == 1 else x[:, :16384].reshape(nt, -1, sc).sum(2)
            sig = np.ones(len(ks))
            b1, k1 = _dd_numpy(xs, ks[::7], sc, ws[::7], sig)
            b2, k2 = _dd_numba(xs, ks[::7], sc, ws[::7], sig)
            checks.append((f"numba kernel == numpy kernel (scrunch {sc})",
                           np.allclose(b1, b2, atol=1e-3) and np.mean(k1 == k2) > 0.999))

    # 5. file path: tiny 2-coarse-channel HDF5 -> .dat (inject.read_dat parses it) + CSV
    import tempfile
    import h5py
    tmp = Path(tempfile.mkdtemp(prefix="fastsearch_selftest_"))
    try:
        nfc = 32768
        rng = np.random.default_rng(11)
        two = np.concatenate([_synthetic(rng, nt, nfc), _synthetic(rng, nt, nfc)], axis=1)
        pl = [dict(ch=40000.0, drift=-0.5, snr=60)]
        inject.plant(two, foff, ts, pl)
        fn = tmp / "tiny.h5"
        with h5py.File(fn, "w") as f:
            ds = f.create_dataset("data", data=two[:, None, :], chunks=(1, 1, 8192))
            for k, v in dict(fch1=1500.0, foff=foff / 1e6, nchans=2 * nfc, tsamp=ts, tstart=57457.5,
                             src_raj=17.96, src_dej=4.74, source_name="SELFTEST", telescope_id=6,
                             nifs=1, nbits=32).items():
                ds.attrs[k] = v
        rows, run = run_file(fn, tmp, 4.0, 8.0, workers=1, n_coarse=2, verbose=False)
        dat = list(tmp.glob("*.dat"))[0]
        back = inject.read_dat(dat, 1500.0, foff / 1e6)
        ok = (len(back) == 1 and abs(back[0]["ch"] - 40000) < 3 and abs(back[0]["drift"] + 0.5) < 0.05)
        print(f"  tiny h5: {len(back)} hit(s) via .dat: "
              + ", ".join(f"ch {b['ch']:.1f} drift {b['drift']:+.4f} snr {b['snr']:.1f}" for b in back))
        checks.append(("file driver: .dat round-trips the one planted tone", ok))
        try:
            from turbo_seti.find_event.find_event import read_dat as ts_read
            df = ts_read(str(dat))
            checks.append(("turboSETI's own read_dat parses the .dat",
                           len(df) == len(back) and abs(float(df["DriftRate"].iloc[0]) + 0.5) < 0.05))
        except ImportError:
            print("  (turbo_seti not installed: skipped its read_dat check)")
    finally:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)

    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    if sys.argv[1:2] == ["selftest"]:
        return selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("h5")
    ap.add_argument("--max-drift", type=float, default=4.0)
    ap.add_argument("--snr", type=float, default=8.0, help="threshold, TRUE radiometer SNR")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--coarse", help="e.g. 0-9,150 (default: all)")
    ap.add_argument("--n-coarse", type=int)
    ap.add_argument("--widths", default="auto", choices=["auto", "full", "1"],
                    help="auto = boxcar + scrunch fast drifts; full = boxcar at full resolution; 1 = no boxcar")
    a = ap.parse_args()
    run_file(a.h5, a.out_dir, a.max_drift, a.snr, a.workers, parse_coarse(a.coarse), a.n_coarse, a.widths)
    return 0


if __name__ == "__main__":
    sys.exit(main())
