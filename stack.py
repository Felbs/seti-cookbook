#!/usr/bin/env python3
"""stack.py - add several observing dates of one star together, in the star's own frame.

One scan reaches a limit set by its own noise. A beacon held at one frequency
in its star's rest frame is in EVERY scan of every date, so adding the dates
(incoherently: powers, not phases) should lower the noise of the sum by up to
sqrt(N). This does it carefully and then MEASURES what it gained:

  per scan   flatten the bandpass (fastsearch.normalise), de-drift by the
             PREDICTED Earth drift at mid-scan (earth_term physics via
             xepoch.frame) plus each value of a small grid of intrinsic drifts,
             sum the 16 spectra (integer-channel shift-and-add: no smoothing),
             and divide by the measured noise -> z spectra.
  RFI mask   channels where an OFF scan of the same epoch rises above
             --mask-z (any drift row), +-3 channels, are blanked in that
             epoch's ONs (RFI is fixed in OUR frame); DC channel blanked.
  frame      shift each scan to the STAR frame at its own mid-scan time:
             f_star = K(t) f_obs, K from xepoch.frame (Earth's velocity on the
             star's propagated direction, the star's RV including the secular
             term), linear interpolation onto one common grid.
  coverage   an epoch whose star-frame coverage does not contain the window
             is REFUSED (different receivers, different centre frequency), not
             silently stacked with NaNs.
  stack      nan-mean over the scans of an epoch, then over epochs, each
             renormalised by its measured noise (MAD).

What is measured, never assumed:
  gain       the same star-frame tone planted (setigen, Doppler-smeared, at
             TRUE radiometer SNR per scan) in every epoch: its z in the stack
             over its mean z in one epoch, after the noise floor of the same
             statistic at empty places is subtracted. Ideal is sqrt(N).
  threshold  pre-registered: max(z_floor, 99th percentile of the window
             statistic at >= 300 empty windows of the same stack).
  false alarms  the identical pipeline with the OFF scans playing the ONs
             (OFF-as-ON control): fraction of windows above the threshold.

    python stack.py selftest          # synthetic epochs, offline (~1 min)
    python stack.py gj699             # three ON scans of 2016-03-10 as PSEUDO-epochs (data/hi_*GJ699*.npz)
    python stack.py real --star GJ725 # real dates fetched by `xepoch.py fetch`
"""
import argparse
import json
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

PREREG = dict(
    drifts_hz_s=[round(x, 3) for x in np.arange(-0.05, 0.0501, 0.01)],
    window_half_hz=3000.0,      # star-frame +-3 kHz: RV error +-0.5 km/s = +-2.4 kHz, plus margin
    z_floor=6.0, null_percentile=99.0, n_null_min=300,
    mask_z=6.0, mask_pad=3,
)


def robust(z):
    z = np.asarray(z, dtype=np.float64)
    med = np.nanmedian(z, axis=-1, keepdims=True)
    mad = 1.4826 * np.nanmedian(np.abs(z - med), axis=-1, keepdims=True)
    return (z - med) / np.where(mad > 0, mad, 1.0)


def dedrift(xn, drifts_obs, foff_hz, tsamp):
    """Integer-channel shift-and-add of a flattened waterfall along each drift
    (Hz/s, observed frame), referenced to MID-scan. Returns (n_drift, nf)."""
    nt, nf = xn.shape
    t = (np.arange(nt) + 0.5 - nt / 2) * tsamp
    out = np.zeros((len(drifts_obs), nf), dtype=np.float32)
    for i, D in enumerate(drifts_obs):
        for k in range(nt):
            sh = int(np.rint(D * t[k] / foff_hz))      # the track sits at ch + sh at time k
            if sh >= 0:
                out[i, :nf - sh] += xn[k, sh:]
            else:
                out[i, -sh:] += xn[k, :nf + sh]
    return out


def prep_scan(x, freqs_hz, tsamp, t_mid, star, drifts=None, dc_idx=(), site=None):
    """One scan -> z spectra (n_drift, nf) at mid-scan, its star-frame grid and K."""
    import fastsearch as F
    import xepoch as X
    drifts = PREREG["drifts_hz_s"] if drifts is None else drifts
    xn, _ = F.normalise(x, dc_idx)
    K, a, _, _ = X.frame(star, [t_mid], site)
    K, a = float(K[0]), float(a[0])
    foff = float(freqs_hz[1] - freqs_hz[0])
    fc = float(np.mean(freqs_hz))
    earth = fc * a / C
    z = robust(dedrift(xn, [earth + d / K for d in drifts], foff, tsamp)).astype(np.float32)
    for c in dc_idx:
        z[:, max(0, c - 1):c + 2] = np.nan
    pad = int(np.ceil(abs(earth + max(abs(d) for d in drifts)) * x.shape[0] * tsamp / abs(foff))) + 1
    z[:, :pad] = np.nan
    z[:, -pad:] = np.nan
    return dict(z=z, f_star=freqs_hz * K, K=K, earth_drift=earth, t_mid=t_mid)


def rfi_mask(off_preps, nf, mask_z=None, pad=None):
    """Topocentric channel mask from OFF scans: any drift row above mask_z."""
    mask_z = PREREG["mask_z"] if mask_z is None else mask_z
    pad = PREREG["mask_pad"] if pad is None else pad
    m = np.zeros(nf, bool)
    for p in off_preps:
        m |= np.nanmax(np.nan_to_num(p["z"], nan=0.0), axis=0) > mask_z
    if pad:
        k = np.flatnonzero(m)
        for d in range(-pad, pad + 1):
            m[np.clip(k + d, 0, nf - 1)] = True
    return m


def to_grid(prep, grid, mask=None):
    """Interpolate a scan's z rows onto the common ascending star-frame grid; NaN
    where masked or outside the scan."""
    f = prep["f_star"]
    z = prep["z"]
    if f[0] > f[-1]:
        f, z = f[::-1], z[:, ::-1]
        mask = None if mask is None else mask[::-1]
    bad = np.isnan(z).any(axis=0)
    if mask is not None:
        bad |= mask
    zz = np.nan_to_num(z, nan=0.0)
    out = np.empty((z.shape[0], grid.size), dtype=np.float32)
    for i in range(z.shape[0]):
        out[i] = np.interp(grid, f, zz[i], left=np.nan, right=np.nan)
    b = np.interp(grid, f, bad.astype(np.float32), left=1.0, right=1.0) > 0
    out[:, b] = np.nan
    return out


def coverage(prep):
    f = prep["f_star"]
    return float(min(f[0], f[-1])), float(max(f[0], f[-1]))


def epoch_spectrum(rows_list):
    """nan-mean of an epoch's scans (each already unit noise), renormalised."""
    a = np.stack(rows_list)
    n = np.sum(~np.isnan(a), axis=0)
    s = np.nanmean(a, axis=0) * np.sqrt(np.maximum(n, 1))
    s[n == 0] = np.nan
    return robust(s).astype(np.float32)


def stack_epochs(epochs, window=None, step_hz=None):
    """epochs: list of dict(name, on=[prep...], mask=topocentric bool or None).
    Returns dict(grid, per_epoch (E, n_drift, ng), stack (n_drift, ng), coverage).
    Refuses (ValueError) an epoch whose coverage misses the window."""
    covs = []
    for e in epochs:
        lo = max(coverage(p)[0] for p in e["on"])
        hi = min(coverage(p)[1] for p in e["on"])
        covs.append((e["name"], lo, hi))
    lo = max(c[1] for c in covs)
    hi = min(c[2] for c in covs)
    if window is not None:
        bad = [c for c in covs if not (c[1] <= window[0] and window[1] <= c[2])]
        if bad:
            raise ValueError("epochs without the window in their star-frame coverage: " +
                             ", ".join(f"{n} ({a/1e6:.4f}-{b/1e6:.4f} MHz)" for n, a, b in bad))
    if hi <= lo:
        raise ValueError("the epochs share no star-frame coverage")
    step = step_hz or abs(epochs[0]["on"][0]["f_star"][1] - epochs[0]["on"][0]["f_star"][0])
    grid = np.arange(lo, hi, step)
    per = []
    for e in epochs:
        per.append(epoch_spectrum([to_grid(p, grid, e.get("mask")) for p in e["on"]]))
    per = np.stack(per)
    n = np.sum(~np.isnan(per), axis=0)
    st = np.nanmean(per, axis=0) * np.sqrt(np.maximum(n, 1))
    st[n == 0] = np.nan
    return dict(grid=grid, per_epoch=per, stack=robust(st).astype(np.float32), coverage=(lo, hi),
                epoch_coverage=covs)


def window_stat(z, grid, centre, half):
    """max over drift rows and channels inside centre +- half."""
    a, b = np.searchsorted(grid, centre - half), np.searchsorted(grid, centre + half)
    w = z[:, a:b]
    if w.size == 0 or np.all(np.isnan(w)):
        return np.nan
    return float(np.nanmax(w))


def box_stat(z, grid, f, row, nch=1):
    """Sum of z over the 2*nch+1 bins around f in one drift row: linear in the
    signal (unlike a max), and it keeps power split between neighbouring bins."""
    j = int(np.searchsorted(grid, f))
    w = z[row, max(0, j - nch):j + nch + 1]
    return float(np.sum(w)) if w.size == 2 * nch + 1 and not np.any(np.isnan(w)) else np.nan


def snr_gain(z_stack, z_single, grid, plants_f, empty_f, row):
    """MEASURED gain: SNR of the box statistic at the plants (mean over plants /
    spread at empty places) in the stack, over the same in one epoch."""
    bs = np.array([box_stat(z_stack, grid, f, row) for f in plants_f])
    b1 = np.array([[box_stat(z, grid, f, row) for f in plants_f] for z in z_single])
    ns = np.array([box_stat(z_stack, grid, f, row) for f in empty_f])
    n1 = np.array([[box_stat(z, grid, f, row) for f in empty_f] for z in z_single])
    sd_s = 1.4826 * np.nanmedian(np.abs(ns - np.nanmedian(ns)))
    sd_1 = 1.4826 * np.nanmedian(np.abs(n1 - np.nanmedian(n1)))
    snr_s = (np.nanmean(bs) - np.nanmedian(ns)) / sd_s
    snr_1 = (np.nanmean(b1) - np.nanmedian(n1)) / sd_1
    return float(snr_s / snr_1), float(snr_s), float(snr_1)


def null_centres(grid, avoid, half, n, rng=None):
    """Centres of NON-overlapping windows tiled over the grid (every one that fits,
    or n of them spread evenly if more fit), clear of `avoid` (list of (lo, hi))."""
    lo, hi = grid[0] + half + 2e3, grid[-1] - half - 2e3
    c = lo
    out = []
    while c <= hi:
        if not any(a - half - 1e3 < c < b + half + 1e3 for a, b in avoid):
            out.append(float(c))
        c += 2 * half
    if len(out) > n:
        out = [out[i] for i in np.linspace(0, len(out) - 1, n).astype(int)]
    return out


# ------------------------------------------------------------------ measurement
def measure(load_scans, star, epochs_def, prereg_window, levels, n_plants=40, seed=5, label="", site=None):
    """Full pre-registered run: threshold from ON nulls, the window, gain from
    plants at each level, false alarms from OFF-as-ON.

    load_scans(): list of dict(name, kind, epoch, x (nt, nf) raw, freqs_hz, tsamp, t_mid, dc)
    epochs_def: ordered epoch names."""
    import inject
    import xepoch as X
    rng = np.random.default_rng(seed)
    t0 = time.time()
    scans = load_scans()
    half = PREREG["window_half_hz"]

    def build(scs, plants=None, as_on="ON"):
        preps = {}
        for s in scs:
            x = s["x"]
            if plants and s["kind"] == as_on:
                x = x.astype(np.float32, copy=True)
                pl = []
                for f_s in plants["f_star"]:
                    ch, dr = X.plant_track(star, s["t_mid"], None, s["tsamp"], x.shape[0], s["freqs_hz"], f_s,
                                           site=site)
                    pl.append(dict(ch=ch, drift=dr, snr=plants["snr"]))
                inject.plant(x, s["freqs_hz"][1] - s["freqs_hz"][0], s["tsamp"], pl)
            preps[s["name"]] = prep_scan(x, s["freqs_hz"], s["tsamp"], s["t_mid"], star,
                                         dc_idx=[s["dc"]] if s.get("dc") is not None else (), site=site)
        ep = []
        other = "OFF" if as_on == "ON" else "ON"
        for e in epochs_def:
            on = [preps[s["name"]] for s in scs if s["epoch"] == e and s["kind"] == as_on]
            off = [preps[s["name"]] for s in scs if s["epoch"] == e and s["kind"] == other]
            nf = on[0]["z"].shape[1]
            ep.append(dict(name=e, on=on, mask=rfi_mask(off, nf) if off else None))
        return ep

    ep_on = build(scans)
    centre = prereg_window[0]
    st = stack_epochs(ep_on, window=(centre - half, centre + half))
    grid = st["grid"]
    masked = float(np.mean(np.isnan(st["stack"][0])))
    nulls_c = null_centres(grid, [(centre - half, centre + half)], half, 400, rng)
    null_stack = np.array([window_stat(st["stack"], grid, c, half) for c in nulls_c])
    null_single = np.array([[window_stat(st["per_epoch"][i], grid, c, half) for c in nulls_c]
                            for i in range(len(ep_on))])
    ok = np.isfinite(null_stack)
    thr = max(PREREG["z_floor"], float(np.percentile(null_stack[ok], PREREG["null_percentile"])))
    thr1 = max(PREREG["z_floor"], float(np.nanpercentile(null_single, PREREG["null_percentile"])))
    print(f"{label}: {len(ep_on)} epochs, grid {grid.size} bins of {grid[1]-grid[0]:.3f} Hz, "
          f"{100*masked:.1f} % masked; nulls {int(ok.sum())} windows: 99th pct {np.percentile(null_stack[ok], 99):.2f} "
          f"-> threshold {thr:.2f} (single epoch {thr1:.2f})  ({time.time()-t0:.0f} s)")
    # the pre-registered window itself
    w_stack = window_stat(st["stack"], grid, centre, half)
    w_single = [window_stat(st["per_epoch"][i], grid, centre, half) for i in range(len(ep_on))]
    # OFF-as-ON control (the same pipeline; the ONs supply the masks)
    ep_off = build(scans, as_on="OFF")
    so = stack_epochs(ep_off, window=(centre - half, centre + half))
    null_off = np.array([window_stat(so["stack"], so["grid"], c, half) for c in nulls_c])
    okf = np.isfinite(null_off)
    far = float(np.mean(null_off[okf] >= thr))
    w_off = window_stat(so["stack"], so["grid"], centre, half)
    print(f"  window {centre/1e6:.6f} MHz +-{half/1e3:.1f} kHz: stack max z {w_stack:.2f} (thr {thr:.2f}); "
          f"single epochs {', '.join(f'{v:.2f}' for v in w_single)}; OFF-as-ON {w_off:.2f}")
    print(f"  OFF-as-ON control: {int((null_off[okf] >= thr).sum())}/{int(okf.sum())} windows above threshold "
          f"= false-alarm rate {100*far:.2f} % per window  ({time.time()-t0:.0f} s)")
    # plants: star-frame tones in every epoch's ONs, clear of masks, nulls kept apart
    mask_any = np.isnan(st["stack"]).any(axis=0)
    cands = rng.uniform(grid[0] + 20e3, grid[-1] - 20e3, 20000)
    plants_f = []
    for f in cands:
        j = np.searchsorted(grid, f)
        if mask_any[max(0, j - 40):j + 40].any() or abs(f - centre) < half + 5e3:
            continue
        if any(abs(f - g) < 12e3 for g in plants_f):
            continue
        plants_f.append(float(f))
        if len(plants_f) >= n_plants:
            break
    row0 = PREREG["drifts_hz_s"].index(0.0)
    empty = null_centres(grid, [(f - 6e3, f + 6e3) for f in plants_f] + [(centre - half, centre + half)],
                         60.0, 400)
    out_levels = {}
    for snr in levels:
        ep_p = build(scans, plants=dict(f_star=plants_f, snr=snr))
        sp = stack_epochs(ep_p, window=(centre - half, centre + half))
        gain, snr_s, snr_1 = snr_gain(sp["stack"], list(sp["per_epoch"]), grid, plants_f, empty, row0)
        ws = np.array([window_stat(sp["stack"], grid, f, half) for f in plants_f])
        w1 = np.array([[window_stat(sp["per_epoch"][i], grid, f, half) for f in plants_f] for i in range(len(ep_p))])
        det_stack = float(np.mean(ws >= thr))
        det_single = float(np.mean(w1 >= thr1))
        out_levels[str(snr)] = dict(box_snr_stack=snr_s, box_snr_single_epoch=snr_1, gain=gain,
                                    detect_stack=det_stack, detect_single_epoch=det_single, n=len(plants_f))
        print(f"  true SNR {snr:>4}/scan: 3-bin SNR stack {snr_s:5.2f} vs one epoch {snr_1:5.2f} "
              f"-> measured gain {gain:.2f} (sqrt N = {np.sqrt(len(ep_p)):.2f}); detected stack "
              f"{100*det_stack:.0f} %, one epoch {100*det_single:.0f} %  ({time.time()-t0:.0f} s)")
    s90 = next((float(k) for k, v in out_levels.items() if v["detect_stack"] >= 0.9), None)
    s90_1 = next((float(k) for k, v in out_levels.items() if v["detect_single_epoch"] >= 0.9), None)
    verdict = ("CANDIDATE (above threshold in the stack and not in OFF-as-ON): try hard to kill it"
               if w_stack >= thr and not (w_off >= thr) else "NULL - nothing above the pre-registered threshold")
    print(f"  VERDICT: {verdict}; 90 % detection at true SNR/scan {s90} (stack) vs {s90_1} (one epoch)")
    return dict(label=label, prereg=PREREG, window_centre_mhz=centre / 1e6, epochs=epochs_def,
                epoch_coverage_mhz=[(n, a / 1e6, b / 1e6) for n, a, b in st["epoch_coverage"]],
                masked_fraction=masked, threshold=thr, threshold_single_epoch=thr1, n_null=int(ok.sum()),
                null_p99=float(np.percentile(null_stack[ok], 99)), window_stack=w_stack, window_single=w_single,
                window_off_as_on=w_off, off_as_on_false_alarm=far, off_as_on_n=int(okf.sum()),
                levels=out_levels, snr90_stack=s90, snr90_single_epoch=s90_1, verdict=verdict,
                plants=len(plants_f))


# ------------------------------------------------------------------ data sources
def gj699_scans():
    """The six 2016-03-10 GJ699 slices, clean part of coarse channel 155 only.
    PSEUDO-epochs: P1 = (ON 0003, OFF 0004), P2 = (0005, 0006), P3 = (0007, 0008)."""
    import fastsearch_recovery as FR
    import xepoch as X
    files = sorted(DATA.glob("hi_*GJ699*.npz"))
    if len(files) < 6:
        raise SystemExit("needs the six data/hi_*GJ699*.npz slices (python hi_starframe.py fetch)")
    out = []
    for i, p in enumerate(files):
        d, m = X.load_npz(p)
        seg, lo, dc = FR.coarse_segment(np.asarray(d, np.float32), m)
        nt = seg.shape[0]
        out.append(dict(name=p.name, kind="OFF" if "_OFF_" in p.name else "ON", epoch=f"P{i // 2 + 1}",
                        x=seg, freqs_hz=m["freqs_hz"][lo:lo + seg.shape[1]], tsamp=float(m["tsamp"]),
                        t_mid=float(m["tstart"]) + nt * float(m["tsamp"]) / 2 / 86400, dc=dc))
        del d
    return out


def real_scans(token):
    import xepoch as X
    files = sorted((DATA / f"xepoch_{token}").glob("*.npz"))
    if not files:
        raise SystemExit(f"nothing in data/xepoch_{token}: run `xepoch.py fetch --star {token}` first")
    scans = [X.make_scan(p.name, "?", 0, 1, 1, 1, [], [], []) for p in files]
    out = []
    for p in files:
        d, m = X.load_npz(p)
        nt, nf = d.shape
        out.append(dict(name=p.name, kind=str(m["kind"]), x=np.asarray(d, np.float32), freqs_hz=m["freqs_hz"],
                        tsamp=float(m["tsamp"]), tstart=float(m["tstart"]),
                        t_mid=float(m["tstart"]) + nt * float(m["tsamp"]) / 2 / 86400,
                        dc=nf // 2 if int(m.get("nfpc", nf)) == nf else None))
    for s, sc in zip(out, scans):
        sc["t_mid"], sc["tstart"] = s["t_mid"], s["tstart"]
    X.group_epochs(scans)
    for s, sc in zip(out, scans):
        s["epoch"] = sc["epoch"]
    return out


# ------------------------------------------------------------------ selftest
def selftest():
    import inject
    import xepoch as X
    t0 = time.time()
    rng = np.random.default_rng(3)
    checks = []
    star = X.Star("GJ699", 269.4520769586187, 4.693364966576667, -801.551, 10362.394, 546.9759, -110.11)
    tsamp, nt, foff, nf = 18.253611008, 16, -2.7939677238464355, 1 << 16
    days = [57457.55, 57470.40, 57489.60, 57510.2]
    f_star0 = 1421.05e6
    Ks = [X.frame(star, [t])[0][0] for t in days]
    top = 0.5 * (max(f_star0 / k for k in Ks) + min(f_star0 / k for k in Ks)) + 0.5 * nf * abs(foff)
    freqs = top + foff * np.arange(nf)
    # dedrift: a tone at the predicted Earth drift sums into ONE channel
    x = np.zeros((nt, nf), np.float32)
    ch, dr = X.plant_track(star, days[0], None, tsamp, nt, freqs, f_star0)
    c_mid = round(ch + dr * nt / 2 * tsamp / foff)
    for k in range(nt):
        x[k, int(c_mid + np.rint(dr * (k + 0.5 - nt / 2) * tsamp / foff))] = 1.0
    s = dedrift(x, [dr], foff, tsamp)[0]
    checks.append(("Earth-drift de-drift puts a drifting track into one channel", s.max() >= nt - 1))

    def make(epoch_days, snr, plants_f, rng_):
        eps = []
        for k, t in enumerate(epoch_days):
            on = []
            for j in range(2):
                tm = t + j * 11 / 1440
                y = rng_.chisquare(64, (nt, nf)).astype(np.float32)
                if plants_f:
                    pl = []
                    for f in plants_f:
                        c, d = X.plant_track(star, tm, None, tsamp, nt, freqs, f)
                        pl.append(dict(ch=c, drift=d, snr=snr))
                    inject.plant(y, foff, tsamp, pl)
                on.append(prep_scan(y, freqs, tsamp, tm, star))
            eps.append(dict(name=f"E{k}", on=on, mask=None))
        return eps
    plants_f = [f_star0 - 40e3 + 8e3 * i for i in range(10)]
    eps = make(days, 4.0, plants_f, rng)
    st = stack_epochs(eps)
    row0 = PREREG["drifts_hz_s"].index(0.0)
    nb = null_centres(st["grid"], [(f - 3e3, f + 3e3) for f in plants_f], 60.0, 300)
    gain, snr_s, snr_1 = snr_gain(st["stack"], list(st["per_epoch"]), st["grid"], plants_f, nb, row0)
    checks.append((f"4 dates x 2 scans, tone at true SNR 4/scan: measured gain {gain:.2f} "
                   f"(3-bin SNR {snr_1:.1f} -> {snr_s:.1f}; sqrt 4 = 2)", 1.6 < gain < 2.2))
    # noise only: the stack has unit spread (no gain from nothing)
    e0 = make(days, 0, [], rng)
    s0 = stack_epochs(e0)
    sd = 1.4826 * np.nanmedian(np.abs(s0["stack"][row0] - np.nanmedian(s0["stack"][row0])))
    tail = float(np.nanmean(s0["stack"][row0] > 4.5))
    checks.append((f"noise-only stack: robust spread {sd:.3f}, P(z > 4.5) {tail:.1e}", abs(sd - 1) < 0.05
                   and tail < 1e-4))
    # a topocentrically fixed tone does NOT stack (it lands kHz apart in the star frame)
    e_t = []
    f_t = freqs[nf // 3]
    for k, t in enumerate(days):
        on = []
        for j in range(2):
            tm = t + j * 11 / 1440
            y = rng.chisquare(64, (nt, nf)).astype(np.float32)
            d_e = prep_scan(y[:, :64], freqs[:64], tsamp, tm, star)["earth_drift"]   # mimics a sky drift
            inject.plant(y, foff, tsamp, [dict(ch=float(nf // 3) - d_e * (nt - 1) * tsamp / 2 / foff,
                                               drift=d_e, snr=10.0)])
            on.append(prep_scan(y, freqs, tsamp, tm, star))
        e_t.append(dict(name=f"T{k}", on=on, mask=None))
    stt = stack_epochs(e_t)
    pk1 = max(np.nanmax(stt["per_epoch"][i]) for i in range(4))
    pks = float(np.nanmax(stt["stack"]))
    checks.append((f"topocentric tone (sky-like drift): one epoch peak {pk1:.1f}, stack peak {pks:.1f} "
                   f"(lands kHz apart, does not add up)", pks < 0.8 * pk1))
    # coverage refusal
    far = [dict(name="X", on=[dict(z=np.zeros((1, 10), np.float32), f_star=np.linspace(2.3e9, 2.31e9, 10), K=1.0)],
                mask=None)]
    try:
        stack_epochs(eps[:1] + far, window=(f_star0 - 3e3, f_star0 + 3e3))
        refused = False
    except ValueError:
        refused = True
    checks.append(("an epoch from another band is refused, not NaN-stacked", refused))
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print(f"selftest {'PASS' if ok else 'FAIL'} ({time.time()-t0:.0f} s)")
    return 0 if ok else 1


def main():
    if sys.argv[1:2] in (["selftest"], []):
        return selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    g = sp.add_parser("gj699")
    g.add_argument("--snr", type=float, nargs="+", default=[2, 3, 4, 5, 6, 8, 10, 12])
    g.add_argument("--json", default=None)
    r = sp.add_parser("real")
    r.add_argument("--star", required=True)
    r.add_argument("--snr", type=float, nargs="+", default=[2, 3, 4, 5, 6, 8, 10, 12])
    r.add_argument("--json", default=None)
    a = ap.parse_args()
    import xepoch as X
    if a.cmd == "gj699":
        star = X.Star.from_targets("GJ699")
        res = measure(gj699_scans, star, ["P1", "P2", "P3"], [F_HI], a.snr,
                      label="GJ699 2016-03-10, three ON scans as PSEUDO-epochs (one session)")
    else:
        star = X.Star.from_targets(a.star)
        scans = real_scans(a.star)
        eps = sorted({s["epoch"] for s in scans}, key=lambda e: min(s["t_mid"] for s in scans if s["epoch"] == e))
        res = measure(lambda: scans, star, eps, [F_HI], a.snr, label=f"{a.star}, {len(eps)} real dates")
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1, default=float))
        print(f"wrote {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
