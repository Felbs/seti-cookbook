#!/usr/bin/env python3
"""periodicity.py - search each channel of a dynamic spectrum for periodic pulse trains.

A transmitter that pulses (a beacon, a rotating beam) puts a pulse train into
one frequency channel: on for a fraction d (the duty cycle) of every period P.
This tool folds every channel's time series at every trial period of a
PRE-REGISTERED range and scores a boxcar matched filter over pulse phase and
width, then compares the scans of an ON/OFF cadence.

Method (ideas from BLIPSS, Suresh et al. 2023, AJ 165, 255, MIT; reimplemented
here in numpy, no code copied):
  1. normalise: each spectrum is divided by its median over its own coarse
     channel (removes common gain changes; also removes signals that cover the
     whole coarse channel), then each channel's time series has its median
     subtracted and is divided by its robust spread (1.4826 MAD); the DC channel
     of each coarse channel is blanked;
  2. fold (numpy brute force, = what an FFA computes): for every trial period p
     on a grid whose step keeps the phase error over the scan below half a bin
     (step 0.5 p^2 / (N nb)), the samples are binned into nb = min(bins_max,
     floor(p / tsamp)) >= bins_min phase bins;
  3. score: every circular window of w = 1 .. floor(ducy_max nb) bins and every
     start phase: snr = sum(window) / sqrt(n (1 - n/N)), n = samples in the window,
     N = all samples (the series is mean-subtracted first). On white noise this has
     mean 0 and spread 1 for every trial (selftest control);
  4. per channel: keep periods >= cutoff, merge neighbours within the period
     tolerance (friends-of-friends), keep the best per cluster, label harmonics
     (P ~ p0/n, 'H') and sub-harmonics (P ~ n p0, 'S') of a stronger period ('F');
     channels next to each other with the same fundamental form one group
     (span > BROAD_CH channels is labelled BROAD: not narrowband);
  5. cadence: every ON fundamental >= on_cutoff is MEASURED in every scan at its
     channel (+-1) and period (+- tolerance) - a forced fold, not "was a candidate
     logged" (KNOWN_ISSUES #9). Detection code, one digit per scan in time order:
     1 if the forced snr >= on_cutoff (ON scans) or >= off_cutoff (OFF scans,
     asymmetric as in BLIPSS: an OFF vetoes more easily). '101010' = ON-only.
  6. look-elsewhere: trials = channels x sum over periods of (phases x widths) are
     counted and printed with the Gaussian expectation; the real null is measured
     on unplanted data, time-shuffled data and with the OFF scans playing ON.

Products. A period P needs >= bins_min samples per period and >= fpmin periods
per scan: tsamp * bins_min <= P <= N tsamp / fpmin. GBT fine-resolution (.0000:
16 x 18 s) reaches nothing (162 s > 96 s) and is refused. Mid-resolution (.0002:
1.07 s x 2.86 kHz, ~280 spectra per 5-min scan) reaches ~9.5-98 s, so the
pre-registered 10-70 s (BLIPSS's range) fits. High-time-resolution (.0001:
~349 us x 366 kHz) reaches ms periods but is ~1.3 GB per scan for 8 channels
per coarse channel (not used here).

Backends: numpy (default, CI). riptide-ffa (pip install riptide-ffa; MIT) is used
only with --backend riptide, for a cross-check; it is never required.

    python periodicity.py selftest
    python periodicity.py search FILE.0002.h5 --out-dir out/ [--coarse 150-165] [--workers 3]
    python periodicity.py cadence DIR --out-dir out/ [--coarse ...]       # six .0002 files
    python periodicity.py inject FILE.0002.h5 --coarse 150-165 --snrs 5,6,7,8,10,12,15,20 --json r.json
    python periodicity.py reach                                          # reachable periods per product
"""
import argparse
import functools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Pre-registered 2026-11-07, before any real-data search with this tool.
PREREG = dict(
    date="2026-11-07",
    product="GBT mid-resolution .0002 (2016: 1.058 s x 2.861 kHz, 277 spectra per scan)",
    period_min_s=10.0, period_max_s=70.0,      # BLIPSS's range
    fpmin=3, bins_min=9, bins_max=10, ducy_max=0.5,
    oversample=0.5,                            # period step = oversample * p^2 / (N nb)
    on_cutoff=8.0, off_cutoff=6.0,             # asymmetric presence bars (BLIPSS: 7 / 6 with riptide's S/N)
    cand_floor=6.0,                            # per-scan candidates are logged from here
    common_mode=True,
)
BROAD_CH = 4            # a group wider than this many channels is labelled BROAD
HARM_MAX = 8            # harmonic orders 2..HARM_MAX


# ------------------------------------------------------------------ utilities
def mem_available_bytes():
    try:
        for ln in open("/proc/meminfo"):
            if ln.startswith("MemAvailable:"):
                return int(ln.split()[1]) * 1024
    except OSError:
        pass
    return None


def reachable(nt, tsamp, fpmin=PREREG["fpmin"], bins_min=PREREG["bins_min"]):
    """(P_min, P_max) in s that a scan of nt samples of tsamp can host."""
    return bins_min * tsamp, nt * tsamp / fpmin


def period_tol(P, T, nb=10):
    """Period tolerance (s) for 'same period': 1.5x the resolution of one scan of
    length T (a one-bin phase slip, P^2 / (T nb)), at least 0.2 s (BLIPSS's epsilon)."""
    return max(0.2, 1.5 * P * P / (T * nb))


def period_grid(nt, tsamp, pmin_s, pmax_s, fpmin=PREREG["fpmin"], bins_min=PREREG["bins_min"],
                bins_max=PREREG["bins_max"], oversample=PREREG["oversample"]):
    """Trial periods in SAMPLES and their bin counts. Raises if the range is unreachable."""
    lo_s, hi_s = reachable(nt, tsamp, fpmin, bins_min)
    a, b = max(pmin_s, lo_s), min(pmax_s, hi_s)
    if a > b:
        raise SystemExit(f"period range {pmin_s:g}-{pmax_s:g} s unreachable with {nt} samples of {tsamp:.4g} s "
                         f"(reachable {lo_s:.4g}-{hi_s:.4g} s: needs >= {bins_min} samples per period and >= "
                         f"{fpmin} periods per scan). Use a product with finer time resolution (.0002 mid-res).")
    out = []
    p = a / tsamp
    while p <= b / tsamp + 1e-9:
        nb = int(min(bins_max, math.floor(p)))
        out.append((p, nb))
        p += oversample * p * p / (nt * nb)
    return out


def trials_per_channel(grid, ducy_max=PREREG["ducy_max"]):
    return int(sum(nb * max(1, int(ducy_max * nb)) for _, nb in grid))


# ------------------------------------------------------------------ data
def read_header(path):
    import h5py
    try:
        import hdf5plugin  # noqa: F401
    except ImportError:
        pass
    with h5py.File(path, "r") as f:
        d = f["data"]
        h = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in d.attrs.items()}
        h["nt"] = d.shape[0]
        h["nchans"] = d.shape[-1]
    return h


def read_block(path, c0, c1):
    import h5py
    try:
        import hdf5plugin  # noqa: F401
    except ImportError:
        pass
    with h5py.File(path, "r") as f:
        d = f["data"]
        x = d[:, 0, c0:c1] if d.ndim == 3 else d[:, c0:c1]
    return np.asarray(x, dtype=np.float32)


def nfpc_of(h):
    import coarse_align
    return coarse_align.nfine_per_coarse(float(h["foff"]), int(h["nchans"]), int(h.get("telescope_id", 6)))


def normalise(raw, nfpc, common=True):
    """raw (nt, nc), nc a multiple of nfpc (block starts on a coarse boundary) ->
    Z (nc, nt) float32, mean-subtracted, unit robust spread; DC channels and dead
    channels are 0. Also returns sigma_raw (nc,): robust spread of each raw channel
    after the common-mode division, in raw units (for planting at measured noise)."""
    nt, nc = raw.shape
    x = raw.astype(np.float32)
    g = None
    if common and nc % nfpc == 0:
        xr = x.reshape(nt, nc // nfpc, nfpc)
        g = np.median(xr, axis=2, keepdims=True)
        g = np.where(g > 0, g, 1.0)
        x = (xr / g).reshape(nt, nc)
    med = np.median(x, axis=0)
    mad = 1.4826 * np.median(np.abs(x - med), axis=0)
    ok = mad > 0
    z = np.zeros_like(x)
    z[:, ok] = (x[:, ok] - med[ok]) / mad[ok]
    if nfpc > 2:
        z[:, nfpc // 2::nfpc] = 0.0                      # DC channel of every coarse channel
    z -= z.mean(axis=0)
    gm = np.ones(nc, np.float32) if g is None else np.repeat(np.median(g[:, :, 0], axis=0), nfpc)[:nc]
    return np.ascontiguousarray(z.T), np.where(ok, mad * gm, 0.0).astype(np.float32)


# ------------------------------------------------------------------ fold + score
def _fold_one(Z, p, nb, wmax, t):
    """Best snr over phase and width at one trial period, all channels. Returns
    (best (nc,), best start bin, best width)."""
    nt = Z.shape[1]
    b = np.floor(((t / p) % 1.0) * nb).astype(np.int64)
    oh = np.zeros((nt, nb), np.float32)
    oh[np.arange(nt), b] = 1.0
    B = Z @ oh                                              # (nc, nb) bin sums
    n = oh.sum(0)
    B2 = np.concatenate([B, B[:, :wmax]], axis=1)
    cs = np.concatenate([np.zeros((B.shape[0], 1), np.float32), np.cumsum(B2, axis=1, dtype=np.float32)], axis=1)
    cn = np.concatenate([[0.0], np.cumsum(np.concatenate([n, n[:wmax]]))])
    J = np.tile(np.arange(nb), wmax)
    W = np.repeat(np.arange(1, wmax + 1), nb)
    N = cn[J + W] - cn[J]
    good = (N > 0) & (N < nt)
    inv = np.zeros_like(N)
    inv[good] = 1.0 / np.sqrt(N[good] * (1.0 - N[good] / nt))
    S = (cs[:, J + W] - cs[:, J]) * inv.astype(np.float32)[None, :]
    k = np.argmax(S, axis=1)
    return S[np.arange(S.shape[0]), k], J[k], W[k]


def fold_search(Z, tsamp, grid, ducy_max=PREREG["ducy_max"], backend="numpy"):
    """Z (nc, nt) -> best (nc, nP) snr, start bin (nc, nP), width in bins (nc, nP)."""
    if backend == "riptide":
        return _riptide_search(Z, tsamp, grid, ducy_max)
    nc, nt = Z.shape
    t = np.arange(nt) + 0.5
    best = np.empty((nc, len(grid)), np.float32)
    bj = np.empty((nc, len(grid)), np.int16)
    bw = np.empty((nc, len(grid)), np.int16)
    for i, (p, nb) in enumerate(grid):
        s, j, w = _fold_one(Z, p, nb, max(1, int(ducy_max * nb)), t)
        best[:, i], bj[:, i], bw[:, i] = s, j, w
    return best, bj, bw


def _riptide_search(Z, tsamp, grid, ducy_max):
    """Optional cross-check: riptide's ffa_search per channel, its S/N mapped onto our
    period grid (max of riptide's periodogram within each of our grid cells). riptide
    defines S/N its own way: compare the two only through injections."""
    from riptide import TimeSeries, ffa_search
    ps = np.array([p for p, _ in grid]) * tsamp
    edges = np.concatenate([[ps[0] - 1e-9], 0.5 * (ps[1:] + ps[:-1]), [ps[-1] + 1e-9]])
    nc = Z.shape[0]
    best = np.full((nc, len(grid)), -np.inf, np.float32)
    for c in range(nc):
        if not Z[c].any():
            continue
        ts = TimeSeries.from_numpy_array(Z[c].astype(np.float32), tsamp=tsamp)
        _, pg = ffa_search(ts, period_min=ps[0], period_max=ps[-1] * 1.0001, fpmin=PREREG["fpmin"],
                           bins_min=PREREG["bins_min"], bins_max=PREREG["bins_max"], ducy_max=ducy_max,
                           deredden=False, already_normalised=False)
        s = pg.snrs.max(axis=1)
        idx = np.clip(np.searchsorted(edges, pg.periods) - 1, 0, len(grid) - 1)
        np.maximum.at(best[c], idx, s.astype(np.float32))
    z = np.zeros(best.shape, np.int16)
    return best, z, z + 1


# ------------------------------------------------------------------ candidates
def label_harmonics(periods, snrs, tol_fn):
    """Greedy, strongest first: each unlabelled period becomes a fundamental 'F';
    periods ~ p0/n are its harmonics 'H', ~ n p0 its sub-harmonics 'S' (n = 2..HARM_MAX).
    Returns labels and the index of the fundamental each one belongs to."""
    order = np.argsort(-np.asarray(snrs))
    lab = np.array([""] * len(periods), dtype=object)
    of = np.full(len(periods), -1)
    for i in order:
        if lab[i]:
            continue
        lab[i], of[i] = "F", i
        p0 = periods[i]
        for j in order:
            if lab[j]:
                continue
            for n in range(2, HARM_MAX + 1):
                if abs(periods[j] - p0 / n) <= tol_fn(p0 / n):
                    lab[j], of[j] = "H", i
                    break
                if abs(periods[j] - n * p0) <= tol_fn(n * p0):
                    lab[j], of[j] = "S", i
                    break
    return lab, of


def channel_candidates(best, bj, bw, grid, tsamp, nt, chan0, thr):
    """Per channel: periods >= thr, friends-of-friends merged within period_tol,
    best per cluster, harmonics labelled. Returns a list of dict rows."""
    T = nt * tsamp
    tol_fn = lambda P: period_tol(P, T)   # noqa: E731
    ps = np.array([p for p, _ in grid]) * tsamp
    nbs = np.array([nb for _, nb in grid])
    rows = []
    cs, ks = np.nonzero(best >= thr)
    if cs.size == 0:
        return rows
    for c in np.unique(cs):
        k = np.sort(ks[cs == c])
        P, S = ps[k], best[c, k]
        # friends-of-friends along period
        clusters, cur = [], [0]
        for i in range(1, len(k)):
            if P[i] - P[cur[-1]] <= tol_fn(P[i]):
                cur.append(i)
            else:
                clusters.append(cur)
                cur = [i]
        clusters.append(cur)
        rep = [cl[int(np.argmax(S[cl]))] for cl in clusters]
        Pr, Sr = P[rep], S[rep]
        lab, of = label_harmonics(Pr, Sr, tol_fn)
        for i, r in enumerate(rep):
            kk = k[r]
            nb = int(nbs[kk])
            rows.append(dict(chan=int(chan0 + c), period_s=float(Pr[i]), snr=float(Sr[i]),
                             duty=float(bw[c, kk]) / nb, phase=(float(bj[c, kk]) + 0.5 * float(bw[c, kk])) / nb,
                             nbins=nb, label=str(lab[i]),
                             fundamental_s=float(Pr[of[i]]) if lab[i] != "F" else float(Pr[i])))
    return rows


def group_channels(rows, T):
    """Fundamentals in neighbouring channels (|dc| <= 1) with the same period form one
    group; returns one row per group (the strongest) with chan_lo/chan_hi/span/broad."""
    f = sorted([r for r in rows if r["label"] == "F"], key=lambda r: (r["chan"], r["period_s"]))
    used = [False] * len(f)
    out = []
    for i, r in enumerate(f):
        if used[i]:
            continue
        grp = [i]
        used[i] = True
        changed = True
        while changed:
            changed = False
            for j in range(len(f)):
                if used[j]:
                    continue
                if any(abs(f[j]["chan"] - f[g]["chan"]) <= 1 and
                       abs(f[j]["period_s"] - f[g]["period_s"]) <= period_tol(f[g]["period_s"], T) for g in grp):
                    grp.append(j)
                    used[j] = True
                    changed = True
        g = max(grp, key=lambda q: f[q]["snr"])
        lo, hi = min(f[q]["chan"] for q in grp), max(f[q]["chan"] for q in grp)
        out.append(dict(f[g], chan_lo=lo, chan_hi=hi, span=hi - lo + 1, broad=int(hi - lo + 1 > BROAD_CH)))
    return out


# ------------------------------------------------------------------ one file
def coarse_blocks(h, coarse, per_block=16):
    nfpc = nfpc_of(h)
    ncoarse = int(h["nchans"]) // nfpc
    todo = list(range(ncoarse)) if coarse is None else [c for c in coarse if 0 <= c < ncoarse]
    blocks = []
    i = 0
    while i < len(todo):          # runs of consecutive coarse channels, at most per_block long
        j = i
        while j + 1 < len(todo) and todo[j + 1] == todo[j] + 1 and j + 1 - i < per_block:
            j += 1
        blocks.append((todo[i] * nfpc, (todo[j] + 1) * nfpc))
        i = j + 1
    return blocks, nfpc


def _search_block(args):
    path, c0, c1, opts = args
    h = read_header(path)
    nfpc = nfpc_of(h)
    raw = read_block(path, c0, c1)
    if opts.get("shuffle"):
        rng = np.random.default_rng(opts.get("seed", 0) + c0)
        raw = raw[rng.permutation(raw.shape[0])]
    Z, _ = normalise(raw, nfpc, opts["common"])
    grid = period_grid(Z.shape[1], float(h["tsamp"]), opts["pmin"], opts["pmax"])
    best, bj, bw = fold_search(Z, float(h["tsamp"]), grid, backend=opts["backend"])
    bm = best.max(axis=1)
    counts = {str(v): int((bm >= v).sum()) for v in (5, 6, 7, 8, 10)}
    rows = channel_candidates(best, bj, bw, grid, float(h["tsamp"]), Z.shape[1], c0, opts["floor"])
    return c0, c1, rows, counts, int((Z.any(axis=1)).sum())


def search_file(path, out_dir=None, coarse=None, workers=1, pmin=PREREG["period_min_s"],
                pmax=PREREG["period_max_s"], backend="numpy", common=True, shuffle=False, quiet=False):
    """Search one file. Returns (rows, summary). rows: per-channel candidates >= cand_floor
    (all labels); summary: trials, counts of channels above 5/6/7/8/10, timing."""
    path = Path(path)
    h = read_header(path)
    nt, tsamp = int(h["nt"]), float(h["tsamp"])
    grid = period_grid(nt, tsamp, pmin, pmax)
    blocks, nfpc = coarse_blocks(h, coarse)
    opts = dict(pmin=pmin, pmax=pmax, backend=backend, common=common, floor=PREREG["cand_floor"], shuffle=shuffle)
    jobs = [(str(path), a, b, opts) for a, b in blocks]
    av = mem_available_bytes()
    per_job = nt * 16 * nfpc * 4 * 8
    if av is not None:
        workers = max(1, min(workers, int((av - (1 << 30)) // max(per_job, 1))))
    t0 = time.time()
    res = []
    if workers > 1 and len(jobs) > 1:
        from multiprocessing import get_context
        with get_context("spawn").Pool(workers) as pool:
            for r in pool.imap_unordered(_search_block, jobs):
                res.append(r)
    else:
        res = [_search_block(j) for j in jobs]
    res.sort()
    rows = [r for x in res for r in x[2]]
    nch = sum(x[4] for x in res)
    counts = {k: sum(x[3][k] for x in res) for k in res[0][3]} if res else {}
    tpc = trials_per_channel(grid)
    from scipy.stats import norm
    T = nt * tsamp
    for r in rows:
        r["freq_mhz"] = float(h["fch1"]) + float(h["foff"]) * r["chan"]
    groups = group_channels(rows, T)
    summ = dict(file=path.name, nt=nt, tsamp=tsamp, T_s=T, nchan_searched=nch, n_periods=len(grid),
                period_range_s=[grid[0][0] * tsamp, grid[-1][0] * tsamp], trials_per_channel=tpc,
                trials_total=tpc * nch, channels_with_max_snr_ge=counts, backend=backend, common_mode=common,
                shuffled=shuffle,
                gaussian_expected_trials_ge={str(v): float(tpc * nch * norm.sf(v)) for v in (6, 8)},
                n_candidates=len(rows), n_fundamental_groups=len(groups),
                n_groups_ge_on_cutoff=sum(g["snr"] >= PREREG["on_cutoff"] for g in groups),
                wall_s=round(time.time() - t0, 1), prereg=PREREG)
    if out_dir:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        stem = path.name[:-3] if path.name.endswith(".h5") else path.stem
        if shuffle:
            stem += ".shuffled"
        _write_csv(out / f"{stem}.periodic.csv", groups)
        (out / f"{stem}.periodic.json").write_text(json.dumps(summ, indent=1))
    if not quiet:
        print(f"{path.name}: {nch} channels x {len(grid)} periods ({grid[0][0]*tsamp:.2f}-{grid[-1][0]*tsamp:.2f} s), "
              f"{tpc} trials/channel = {tpc*nch:.3g} trials; channels with max snr >= 6/8: {counts.get('6')}/"
              f"{counts.get('8')} (Gaussian expectation for independent trials: <= "
              f"{summ['gaussian_expected_trials_ge']['6']:.2g}/{summ['gaussian_expected_trials_ge']['8']:.2g}); "
              f"{len(groups)} fundamental groups, {summ['n_groups_ge_on_cutoff']} >= {PREREG['on_cutoff']}; "
              f"{summ['wall_s']:.0f} s", flush=True)
    return groups, summ


def _write_csv(path, rows):
    import csv
    cols = ["chan", "chan_lo", "chan_hi", "span", "broad", "freq_mhz", "period_s", "snr", "duty", "phase", "nbins",
            "label", "fundamental_s"]
    with open(path, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _read_csv(path):
    import csv
    out = []
    for r in csv.DictReader(open(path)):
        out.append({k: (float(v) if k not in ("label",) else v) for k, v in r.items()})
        for k in ("chan", "chan_lo", "chan_hi", "span", "broad", "nbins"):
            out[-1][k] = int(out[-1][k])
    return out


# ------------------------------------------------------------------ forced measurement + cadence
@functools.lru_cache(maxsize=96)
def _coarse_Z(path, cc, nfpc, common):
    """Normalised coarse channel cc of one scan (cached: events are measured channel-sorted)."""
    h = read_header(path)
    Z, _ = normalise(read_block(path, cc * nfpc, (cc + 1) * nfpc), nfpc, common)
    return Z, float(h["tsamp"])


def forced_snr(path, chan, period_s, tol_s, nfpc, common=True, half=1):
    """Best snr at channels chan-half..chan+half and periods within +-tol (fine grid),
    in one scan. The block read is the whole coarse channel (so the common-mode and
    the robust spread are the same as in the search)."""
    cc = chan // nfpc
    Z, tsamp = _coarse_Z(str(path), cc, nfpc, common)
    lo, hi = max(0, chan - half - cc * nfpc), min(nfpc, chan + half + 1 - cc * nfpc)
    Zs = Z[lo:hi]
    nt = Z.shape[1]
    ps = np.linspace(max(period_s - tol_s, 1.0 * tsamp), period_s + tol_s, 9) / tsamp
    t = np.arange(nt) + 0.5
    best = -np.inf
    for p in ps:
        nb = int(min(PREREG["bins_max"], max(PREREG["bins_min"], math.floor(p))))
        if p < nb:
            continue
        s, _, _ = _fold_one(Zs, p, nb, max(1, int(PREREG["ducy_max"] * nb)), t)
        best = max(best, float(s.max()))
    return best


def cadence_codes(files, kinds, cand_per_file, nfpc, T, as_on="ON", common=True):
    """For every group >= on_cutoff found in a scan whose kind == as_on, measure it in all
    scans and build the detection code. as_on="OFF" is the OFF-as-ON null: the OFF scans
    play the ON role (their candidates seed, their bar is on_cutoff, the ON scans' bar
    is off_cutoff)."""
    out = []
    for i, (cands, k) in enumerate(zip(cand_per_file, kinds)):
        if k != as_on:
            continue
        for g in sorted(cands, key=lambda r: r["chan"]):
            if g["snr"] < PREREG["on_cutoff"] or g["label"] != "F":
                continue
            tol = period_tol(g["period_s"], T)
            snrs, code = [], ""
            for j, (f, kj) in enumerate(zip(files, kinds)):
                s = g["snr"] if j == i else forced_snr(f, g["chan"], g["period_s"], tol, nfpc, common)
                bar = PREREG["on_cutoff"] if kj == as_on else PREREG["off_cutoff"]
                snrs.append(round(s, 2))
                code += "1" if s >= bar else "0"
            role = "".join("1" if kj == as_on else "0" for kj in kinds)
            verdict = ("ON-ONLY" if code == role else "IN-OFF" if any(c == "1" and r == "0" for c, r in zip(code, role))
                       else "PARTIAL")
            out.append(dict(g, found_in=i, code=code, verdict=verdict, snr_per_scan=snrs,
                            as_on=as_on))
    return out


def dedupe_events(ev, T):
    """The same signal seeded from two ON scans gives two rows: keep the strongest."""
    ev = sorted(ev, key=lambda e: -e["snr"])
    keep = []
    for e in ev:
        if any(abs(e["chan"] - k["chan"]) <= 1 and abs(e["period_s"] - k["period_s"]) <= period_tol(e["period_s"], T)
               for k in keep):
            continue
        keep.append(e)
    return keep


def write_cadence_outputs(out, res):
    """periodic_cadence.json keeps counts and every ON-ONLY / PARTIAL event in full;
    every event (IN-OFF included) goes one line each into periodic_cadence_events.csv."""
    import csv
    cols = ["as_on", "verdict", "code", "found_in", "chan", "chan_lo", "chan_hi", "span", "broad", "freq_mhz",
            "period_s", "snr", "duty", "phase", "snr_per_scan"]
    with open(Path(out) / "periodic_cadence_events.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for e in res["events"] + res["off_as_on_events"]:
            w.writerow(dict(e, snr_per_scan=" ".join(f"{v:.2f}" for v in e["snr_per_scan"]),
                            freq_mhz=f"{e['freq_mhz']:.6f}", period_s=f"{e['period_s']:.3f}", snr=f"{e['snr']:.2f}"))
    slim = dict(res, events=[e for e in res["events"] if e["verdict"] != "IN-OFF"],
                off_as_on_events=[e for e in res["off_as_on_events"] if e["verdict"] != "IN-OFF"],
                note="IN-OFF events are only in periodic_cadence_events.csv")
    (Path(out) / "periodic_cadence.json").write_text(json.dumps(slim, indent=1))


def run_cadence(d, out_dir, coarse=None, workers=1, backend="numpy", common=True, null=True):
    import fastsearch_cadence as FC
    scans, problems = FC.order_cadence(d)
    if problems:
        raise SystemExit("NOT A VALID CADENCE: " + "; ".join(problems))
    files = [s["path"] for s in scans]
    kinds = [s["kind"] for s in scans]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cands, summs = [], []
    for f in files:
        stem = f.name[:-3]
        csvp, js = out / f"{stem}.periodic.csv", out / f"{stem}.periodic.json"
        if csvp.exists() and js.exists() and json.loads(js.read_text()).get("backend") == backend:
            print(f"have {csvp.name} - skipped")
            cands.append(_read_csv(csvp))
            summs.append(json.loads(js.read_text()))
            continue
        g, s = search_file(f, out, coarse, workers, backend=backend, common=common)
        cands.append(g)
        summs.append(s)
    h = read_header(files[0])
    nfpc = nfpc_of(h)
    T = int(h["nt"]) * float(h["tsamp"])
    t0 = time.time()
    ev = dedupe_events(cadence_codes(files, kinds, cands, nfpc, T, "ON", common), T)
    evn = dedupe_events(cadence_codes(files, kinds, cands, nfpc, T, "OFF", common), T) if null else []
    res = dict(files=[f.name for f in files], kinds=kinds, prereg=PREREG, trials_per_scan=[s["trials_total"] for s in summs],
               events=ev, off_as_on_events=evn,
               counts={v: sum(e["verdict"] == v for e in ev) for v in ("ON-ONLY", "IN-OFF", "PARTIAL")},
               off_as_on_counts={v: sum(e["verdict"] == v for e in evn) for v in ("ON-ONLY", "IN-OFF", "PARTIAL")},
               forced_s=round(time.time() - t0, 1))
    write_cadence_outputs(out, res)
    print(f"cadence {' '.join(kinds)}: {len(ev)} ON events >= {PREREG['on_cutoff']}: {res['counts']}; "
          f"OFF-as-ON null: {len(evn)} events: {res['off_as_on_counts']}")
    for e in ev:
        print(f"  {e['verdict']:8s} code {e['code']} chan {e['chan']} ({e['freq_mhz']:.6f} MHz, span {e['span']}"
              f"{' BROAD' if e['broad'] else ''}) P {e['period_s']:.2f} s duty {e['duty']:.2f} snr per scan "
              f"{e['snr_per_scan']}")
    return res


# ------------------------------------------------------------------ injection
def plant_trains(raw, sigma_raw, plants, tsamp):
    """Add boxcar pulse trains to raw (nt, nc) in place. plant: dict(chan (in raw),
    period_s, duty, phase (pulse centre, 0-1), snr = INTEGRATED matched-filter SNR).
    Amplitude per on-pulse sample A = snr * sigma / sqrt(n_on), sigma = the channel's
    measured robust spread; BLIPSS's 'pulse SNR' (per-sample peak) is A / sigma."""
    nt = raw.shape[0]
    t = np.arange(nt) + 0.5
    for p in plants:
        ph = (t * tsamp / p["period_s"]) % 1.0
        d = ((ph - p["phase"] + 0.5) % 1.0) - 0.5
        on = np.abs(d) < 0.5 * p["duty"]
        n_on = int(on.sum())
        if n_on == 0:
            p["n_on"], p["pulse_snr"] = 0, 0.0
            continue
        A = p["snr"] * float(sigma_raw[p["chan"]]) / math.sqrt(n_on * (1 - n_on / nt))
        raw[on, p["chan"]] += A
        p["n_on"], p["pulse_snr"] = n_on, A / float(sigma_raw[p["chan"]])
    return raw


def score_plants(plants, rows, T, c0=0):
    """found (fundamental at the planted period within tolerance, +-1 channel, snr >=
    on_cutoff), found-harmonic (a harmonic/sub-harmonic relation instead), or missed."""
    st = []
    for p in plants:
        best = "missed"
        for r in rows:
            if abs(r["chan"] - (p["chan"] + c0)) > 1 or r["snr"] < PREREG["on_cutoff"]:
                continue
            P = r["period_s"]
            if abs(P - p["period_s"]) <= period_tol(p["period_s"], T):
                best = "found"
                break
            if any(abs(P - p["period_s"] / n) <= period_tol(p["period_s"] / n, T) or
                   abs(P - p["period_s"] * n) <= period_tol(p["period_s"] * n, T) for n in range(2, HARM_MAX + 1)):
                best = "found-harmonic"
        st.append(best)
    return st


def inject_run(path, coarse, snrs, per_level, duty=0.1, seed=1, backend="numpy", pmin=PREREG["period_min_s"],
               pmax=PREREG["period_max_s"], common=True, quiet=False):
    """Plant pulse trains into one block of REAL data at its measured per-channel noise
    and search it with the production settings. One plant per channel; channels with a
    candidate in the unplanted search (>= cand_floor) and the DC channels are avoided."""
    h = read_header(path)
    nfpc = nfpc_of(h)
    blocks, _ = coarse_blocks(h, coarse, per_block=len(coarse))
    c0, c1 = blocks[0]
    raw0 = read_block(path, c0, c1)
    tsamp, nt = float(h["tsamp"]), raw0.shape[0]
    T = nt * tsamp
    Z, sig = normalise(raw0, nfpc, common)
    grid = period_grid(nt, tsamp, pmin, pmax)
    best, bj, bw = fold_search(Z, tsamp, grid, backend=backend)
    rows0 = channel_candidates(best, bj, bw, grid, tsamp, nt, 0, PREREG["cand_floor"])
    busy = {r["chan"] + o for r in rows0 for o in (-2, -1, 0, 1, 2)}
    rng = np.random.default_rng(seed)
    free = [c for c in range(raw0.shape[1]) if c % nfpc not in (nfpc // 2 - 1, nfpc // 2, nfpc // 2 + 1)
            and c % nfpc > 8 and c % nfpc < nfpc - 8 and c not in busy and sig[c] > 0]
    rng.shuffle(free)
    need = len(snrs) * per_level
    if len(free) < 3 * need:
        raise SystemExit(f"only {len(free)} free channels for {need} plants (spaced 3 apart)")
    chans = sorted(free[:3 * need])[::3][:need]
    rng.shuffle(chans)
    plants = []
    for i, c in enumerate(chans):
        plants.append(dict(chan=int(c), snr=float(snrs[i % len(snrs)]), duty=duty,
                           period_s=float(np.exp(rng.uniform(np.log(pmin), np.log(pmax)))),
                           phase=float(rng.uniform(0, 1))))
    raw = plant_trains(raw0.copy(), sig, plants, tsamp)
    Z1, _ = normalise(raw, nfpc, common)
    t0 = time.time()
    b1, j1, w1 = fold_search(Z1, tsamp, grid, backend=backend)
    rows = channel_candidates(b1, j1, w1, grid, tsamp, nt, 0, PREREG["cand_floor"])
    st = score_plants(plants, rows, T)
    table = {}
    for p, s in zip(plants, st):
        c = table.setdefault(f"{p['snr']:g}", dict(n=0, found=0, harmonic=0, pulse_snr=[]))
        c["n"] += 1
        c["found"] += s == "found"
        c["harmonic"] += s == "found-harmonic"
        c["pulse_snr"].append(p["pulse_snr"])
    for c in table.values():
        c["pulse_snr_median"] = float(np.median(c.pop("pulse_snr")))
    planted = {p["chan"] for p in plants}
    stray = [r for r in rows if r["snr"] >= PREREG["on_cutoff"] and r["label"] == "F"
             and all(abs(r["chan"] - c) > 1 for c in planted)]
    res = dict(file=Path(path).name, coarse=list(coarse), seed=seed, duty=duty, backend=backend,
               period_range_s=[pmin, pmax], n_plants=len(plants), table=table,
               null_candidates_unplanted=sum(r["snr"] >= PREREG["on_cutoff"] and r["label"] == "F" for r in rows0),
               stray_after_planting=len(stray), search_s=round(time.time() - t0, 1), plants=plants, status=st)
    if not quiet:
        print(f"{Path(path).name} coarse {coarse[0]}-{coarse[-1]}: {len(plants)} plants, duty {duty}, periods "
              f"{pmin:g}-{pmax:g} s log-uniform, backend {backend}; unplanted fundamentals >= {PREREG['on_cutoff']}: "
              f"{res['null_candidates_unplanted']}")
        print("  integrated SNR : found / harmonic-only / n   (BLIPSS-style per-sample pulse SNR)")
        for k in sorted(table, key=float):
            c = table[k]
            print(f"  {float(k):6.1f} : {c['found']:3d} / {c['harmonic']:3d} / {c['n']:3d}   "
                  f"({c['pulse_snr_median']:.2f})")
    return res


# ------------------------------------------------------------------ selftest
def _synthetic_h5(path, rng, nt=277, nc=2048, tsamp=1.058013184, kind="ON", tstart=57457.5, plants=(),
                  name_target="SYNTH"):
    """Mid-res-like file: chi-square noise (~3000 dof per sample as on .0002), a
    bandpass, a common per-spectrum gain wobble, coarse channels of 1024."""
    import h5py
    f = np.linspace(-1, 1, 1024)
    bp = np.tile((1.0 - 0.4 * f ** 2), nc // 1024)
    gain = 1 + 0.02 * np.sin(np.arange(nt) / 7.0) + 0.01 * rng.standard_normal(nt)
    x = (rng.chisquare(3000, (nt, nc)) / 3000).astype(np.float32) * bp[None, :] * gain[:, None]
    x *= 1e7
    x[:, 512::1024] *= 30.0
    if plants:
        xn = x / (bp[None, :] * gain[:, None])
        sig = 1.4826 * np.median(np.abs(xn - np.median(xn, 0)), 0) * bp
        plant_trains(x, sig.astype(np.float32), [dict(p) for p in plants], tsamp)
    with h5py.File(path, "w") as fo:
        ds = fo.create_dataset("data", data=x[:, None, :], chunks=(5, 1, 1024))
        for k, v in dict(fch1=1650.0, foff=-0.00286102294921875, nchans=nc, tsamp=tsamp, tstart=tstart,
                         src_raj=17.96, src_dej=4.74, source_name=name_target, telescope_id=6, nifs=1,
                         nbits=32).items():
            ds.attrs[k] = v


def selftest():
    import shutil
    import tempfile
    checks = []
    rng = np.random.default_rng(5)
    tsamp, nt = 1.058013184, 277
    # 0. reachability: fine-res refused, mid-res fits 10-70 s
    try:
        period_grid(16, 17.986, 10, 70)
        ok0 = False
    except SystemExit:
        ok0 = True
    lo, hi = reachable(nt, tsamp)
    checks.append(("fine-res (16 x 18 s) cannot host 10-70 s and is refused", ok0))
    checks.append(("mid-res (277 x 1.058 s) reaches 9.5-97.7 s", abs(lo - 9.52) < 0.05 and abs(hi - 97.7) < 0.2))
    grid = period_grid(nt, tsamp, 10, 70)
    tpc = trials_per_channel(grid)
    print(f"  grid: {len(grid)} periods 10-70 s, {tpc} trials per channel")

    # 1. statistic on white noise: mean ~0, spread ~1 at fixed (period, phase, width)
    Z = rng.standard_normal((4000, nt)).astype(np.float32)
    Z -= Z.mean(1, keepdims=True)
    t = np.arange(nt) + 0.5
    zs = []
    for p, nb in grid[::50]:
        b = np.floor(((t / p) % 1.0) * nb).astype(int)
        on = b < 2
        n = on.sum()
        zs.append(Z[:, on].sum(1) / math.sqrt(n * (1 - n / nt)))
    zs = np.concatenate(zs)
    print(f"  fixed-trial statistic on white noise: mean {zs.mean():+.3f}, std {zs.std():.3f}")
    checks.append(("fixed-trial snr on noise: |mean| < 0.03, std within 3 % of 1",
                   abs(zs.mean()) < 0.03 and abs(zs.std() - 1) < 0.03))
    # max over the grid on noise: count above 6 and 8 in 4000 channels
    best, bj, bw = fold_search(Z, tsamp, grid)
    bm = best.max(1)
    from scipy.stats import norm
    print(f"  white noise, 4000 channels: max snr >= 5/6/8 in {(bm >= 5).sum()}/{(bm >= 6).sum()}/{(bm >= 8).sum()} "
          f"channels (independent-trial bound at 6: {4000 * tpc * norm.sf(6):.2f})")
    checks.append(("white noise: no channel >= 8, <= 3 channels >= 6", (bm >= 8).sum() == 0 and (bm >= 6).sum() <= 3))

    # 2. a planted train is found at its period; harmonic labelling
    Zp = rng.standard_normal((50, nt)).astype(np.float32)
    raw = Zp.T.copy() + 10.0
    pl = [dict(chan=10, period_s=23.0, duty=0.1, phase=0.3, snr=15.0),
          dict(chan=30, period_s=41.0, duty=0.25, phase=0.7, snr=12.0)]
    plant_trains(raw, np.ones(50, np.float32), pl, tsamp)
    Zn, _ = normalise(raw, 10 ** 9, common=False)
    b, j, w = fold_search(Zn, tsamp, grid)
    rows = channel_candidates(b, j, w, grid, tsamp, nt, 0, 6.0)
    st = score_plants(pl, rows, nt * tsamp)
    f10 = [r for r in rows if r["chan"] == 10 and r["label"] == "F"]
    print(f"  plants: {st}; chan 10 fundamental {f10[0]['period_s']:.2f} s snr {f10[0]['snr']:.1f} duty "
          f"{f10[0]['duty']:.2f}; labels in chan 10: {sorted({r['label'] for r in rows if r['chan'] == 10})}"
          if f10 else f"  plants: {st}")
    checks.append(("two planted trains found at their periods", st == ["found", "found"]))
    # snr loss: 40 trains of integrated snr 15, duty 0.1, random period and phase
    Zq = rng.standard_normal((120, nt)).astype(np.float32)
    rq = Zq.T.copy() + 10.0
    pq = [dict(chan=3 * i, period_s=float(np.exp(rng.uniform(np.log(10), np.log(70)))), duty=0.1,
               phase=float(rng.uniform()), snr=15.0) for i in range(40)]
    plant_trains(rq, np.ones(120, np.float32), pq, tsamp)
    Zq, _ = normalise(rq, 10 ** 9, common=False)
    bq, _, _ = fold_search(Zq, tsamp, grid)
    ratio = np.median(bq.max(1)[::3] / 15.0)
    print(f"  40 trains of integrated snr 15, duty 0.1: recovered / planted median {ratio:.2f} "
          f"(phase-bin and period-grid loss)")
    checks.append(("recovered snr: median 0.6-1.1 of the planted integrated snr", 0.6 < ratio < 1.1))
    lab, of = label_harmonics(np.array([20.0, 10.0, 40.0, 33.0]), np.array([10.0, 6.0, 7.0, 8.0]),
                              lambda P: 0.2)
    checks.append(("harmonic labels: 20 F, 10 H, 40 S, 33 F", list(lab) == ["F", "H", "S", "F"]))

    # 3. cadence on six tiny files: ON-only train -> 101010; all-scan train -> IN-OFF;
    #    OFF-as-ON null finds the all-scan one as IN-OFF and nothing ON-ONLY
    tmp = Path(tempfile.mkdtemp(prefix="periodic_selftest_"))
    try:
        starts = [48006, 48334, 48661, 48988, 49314, 49640]
        on_tr = dict(chan=300, period_s=27.0, duty=0.1, phase=0.4, snr=14.0)
        all_tr = dict(chan=1500, period_s=17.0, duty=0.2, phase=0.1, snr=14.0)
        for i, s in enumerate(starts):
            kind = "ON" if i % 2 == 0 else "OFF"
            tgt = "SYNTH" if kind == "ON" else "SYNTH_OFF"
            pls = [all_tr] + ([on_tr] if kind == "ON" else [])
            _synthetic_h5(tmp / f"blcXX_guppi_57457_{s}_{tgt}_{i + 3:04d}.gpuspec.0002.h5", rng, kind=kind,
                          tstart=57457 + s / 86400, plants=pls, name_target=tgt)
        res = run_cadence(tmp, tmp / "out", workers=1)
        ev = {e["chan"]: e for e in res["events"]}
        checks.append(("ON-only train: code 101010, ON-ONLY", 300 in ev and ev[300]["code"] == "101010"
                       and ev[300]["verdict"] == "ON-ONLY"))
        checks.append(("train in every scan: IN-OFF", 1500 in ev and ev[1500]["verdict"] == "IN-OFF"))
        checks.append(("no other ON events", set(ev) <= {300, 1500}))
        checks.append(("OFF-as-ON null: nothing ON-ONLY", res["off_as_on_counts"]["ON-ONLY"] == 0))
        js = json.loads((tmp / "out" / "blcXX_guppi_57457_48006_SYNTH_0003.gpuspec.0002.periodic.json").read_text())
        checks.append(("summary counts trials", js["trials_total"] == js["trials_per_channel"] * js["nchan_searched"]))
        # injection gate on a synthetic file (the real-data gate is in PERIODICITY_RESULT.md)
        r = inject_run(tmp / "blcXX_guppi_57457_48334_SYNTH_OFF_0004.gpuspec.0002.h5", [0, 1], [6, 15], 40,
                       quiet=True)
        f15 = r["table"]["15"]["found"] / r["table"]["15"]["n"]
        f6 = r["table"]["6"]["found"] / r["table"]["6"]["n"]
        print(f"  synthetic injection: found {f6:.2f} at integrated snr 6 (below the bar of 8), {f15:.2f} at 15")
        checks.append(("injection: >= 95 % found at snr 15, <= 20 % at snr 6", f15 >= 0.95 and f6 <= 0.2))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    try:
        import riptide  # noqa: F401
        print("  riptide-ffa importable: cross-check available (--backend riptide)")
    except ImportError:
        print("  riptide-ffa not installed: numpy backend only (that is the default)")
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_reach():
    rows = [("GBT fine-res .0000 (2016)", 16, 17.986), ("GBT mid-res .0002 (2016)", 277, 1.058013),
            ("GBT high-time-res .0001 (2016)", int(288 / 349.5e-6), 349.5e-6)]
    for name, nt, ts in rows:
        lo, hi = reachable(nt, ts)
        print(f"{name:32s} {nt:8d} x {ts:.4g} s: periods {lo:.4g} - {hi:.4g} s"
              + ("  (none)" if lo > hi else ""))


def main():
    if sys.argv[1:2] == ["selftest"]:
        return selftest()
    if sys.argv[1:2] == ["reach"]:
        return cmd_reach()
    import fastsearch as F
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for nm in ("search", "cadence", "inject"):
        p = sub.add_parser(nm)
        p.add_argument("path")
        p.add_argument("--coarse", help="coarse channels, e.g. 150-165 (default all; inject needs it)")
        p.add_argument("--backend", default="numpy", choices=["numpy", "riptide"])
        p.add_argument("--no-common", action="store_true", help="do not divide each spectrum by its coarse median")
        if nm != "inject":
            p.add_argument("--out-dir", required=True)
            p.add_argument("--workers", type=int, default=1)
        if nm == "search":
            p.add_argument("--shuffle", action="store_true", help="time-shuffle each block first (null control)")
        if nm == "inject":
            p.add_argument("--snrs", default="5,6,7,8,9,10,12,15,20")
            p.add_argument("--per-level", type=int, default=40)
            p.add_argument("--duty", type=float, default=0.1)
            p.add_argument("--seed", type=int, default=1)
            p.add_argument("--json")
    a = ap.parse_args()
    coarse = F.parse_coarse(a.coarse)
    if a.cmd == "search":
        search_file(a.path, a.out_dir, coarse, a.workers, backend=a.backend, common=not a.no_common,
                    shuffle=a.shuffle)
    elif a.cmd == "cadence":
        run_cadence(a.path, a.out_dir, coarse, a.workers, a.backend, not a.no_common)
    else:
        if not coarse:
            raise SystemExit("inject needs --coarse")
        r = inject_run(a.path, coarse, [float(v) for v in a.snrs.split(",")], a.per_level, a.duty, a.seed,
                       a.backend, common=not a.no_common)
        if a.json:
            Path(a.json).write_text(json.dumps(r, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
