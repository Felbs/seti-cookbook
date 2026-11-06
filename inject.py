#!/usr/bin/env python3
"""inject.py - plant fake alien signals into REAL data, then score the search fairly.

A limit ("nothing brighter than X") only means something if we have shown the
search finds an X-strength signal when one is there. This module does that
the way the field does (Choza et al. 2024, Ma et al. 2023, theoxenia 2026):

  plant()  - setigen (Brzycki et al. 2022, MIT) draws each tone with Doppler
             smearing: a drifting tone spreads its power over every channel it
             crosses during each integration, as a real one does (law L10).
  match()  - ONE-TO-ONE assignment (Hungarian) along the whole track: a hit
             counts for a planted tone only if its line stays within a few
             channels of the tone's path from start to end, and no hit is
             credited to two tones.
  explain_misses() - every miss gets a reason: found at the WRONG DRIFT,
             SUPPRESSED (a stronger hit within turboSETI's neighbour radius
             swallowed it), or BELOW THRESHOLD.

Strength is TRUE radiometer SNR: the planted tone's time-summed power over the
channel-to-channel noise with the bandpass removed (radiometer_noise). That is
NOT turboSETI's SNR, which divides by a spread that includes the bandpass shape
(turboseti_noise) - on real GBT coarse channels ~4x larger.

    python inject.py selftest       # synthetic noise, no download, no turboSETI
"""
import numpy as np


def robust_std(x):
    x = np.asarray(x, dtype=np.float64)
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))


def radiometer_noise(data):
    """Noise of the time-summed spectrum from channel to channel, with the
    bandpass removed (65-channel running median). On real GBT fine-res data
    this matches the radiometer equation to ~5 % (GJ699 2016-03-10: 0.0235 vs
    0.0248 of the mean), so SNR in these units converts straight to Jy/EIRP
    with reach.py."""
    from scipy.ndimage import median_filter
    s = np.asarray(data, dtype=np.float64).sum(0)
    return robust_std(s - median_filter(s, size=65))


def turboseti_noise(data):
    """What turboSETI 2.3.2 divides by (helper_functions.comp_stats): plain std
    of the time-summed spectrum after dropping the top and bottom 5 % - which
    INCLUDES the coarse channel's bandpass shape."""
    s = np.asarray(data, dtype=np.float64).sum(0)
    lo, hi = np.percentile(s, [5, 95])
    return s[(s >= lo) & (s <= hi)].std()


def level_for(data, snr):
    """Per-spectrum power that gives a ZERO-drift, channel-centred tone this
    TRUE (radiometer) integrated SNR. data: (nt, nchan) real data."""
    return snr * radiometer_noise(data) / data.shape[0]


def plant(data, foff_hz, tsamp, plants, seed=0):
    """Add tones to a (nt, nchan) array in place and return it.

    plants: list of dict(ch=start channel at t=0 (float), drift=Hz/s, snr=...)
    Channel index follows the data's own order; foff_hz is the (signed)
    channel step, so with BL's negative foff a positive drift moves to LOWER
    channel numbers. Internally the slice is flipped to ascending frequency so
    setigen's index/frequency convention cannot bite us.
    """
    import setigen as stg
    nt, nc = data.shape
    asc = foff_hz > 0
    work = data if asc else data[:, ::-1]
    df = abs(foff_hz)
    fr = stg.Frame.from_data(df=df, dt=tsamp, fch1=0.0, ascending=True,
                             data=np.ascontiguousarray(work, dtype=np.float64), seed=seed)
    for p in plants:
        ch = p["ch"] if asc else (nc - 1 - p["ch"])
        lvl = p.get("level") or level_for(data, p["snr"])
        f0 = fr.get_frequency(0) + ch * df
        f1 = f0 + p["drift"] * nt * tsamp
        pad = 8 * df + abs(p["drift"]) * tsamp
        fr.add_signal(stg.constant_path(f_start=f0, drift_rate=p["drift"]),
                      stg.constant_t_profile(level=lvl),
                      stg.box_f_profile(width=df),
                      stg.constant_bp_profile(level=1.0),
                      bounding_f_range=(min(f0, f1) - pad, max(f0, f1) + pad),
                      doppler_smearing=True, smearing_subsamples=
                      max(10, int(np.ceil(abs(p["drift"]) * tsamp / df * 4))))
    out = fr.data if asc else fr.data[:, ::-1]
    data[:] = out.astype(data.dtype)
    return data


def _track(ch0, drift_hz_s, t, foff_hz):
    return ch0 + drift_hz_s * t / foff_hz


def match(plants, hits, foff_hz, tsamp, nt, tol_ch=3.0):
    """One-to-one path match.

    plants: list of dict(ch, drift, snr); hits: list of dict(ch, drift, snr)
    with ch = the hit's channel at t=0 (turboSETI start index).
    A pair is allowed if the two straight tracks stay within tol_ch channels
    plus the planted tone's own per-spectrum smear, at both the first and the
    last spectrum. Returns list (one per plant) of matched hit index or None.
    """
    from scipy.optimize import linear_sum_assignment
    T = (nt - 1) * tsamp
    if not plants or not hits:
        return [None] * len(plants)
    big = 1e9
    cost = np.full((len(plants), len(hits)), big)
    for i, p in enumerate(plants):
        smear = abs(p["drift"]) * tsamp / abs(foff_hz)
        tol = tol_ch + smear
        for j, h in enumerate(hits):
            d0 = abs(h["ch"] - p["ch"])
            d1 = abs(_track(h["ch"], h["drift"], T, foff_hz) - _track(p["ch"], p["drift"], T, foff_hz))
            if max(d0, d1) <= tol:
                cost[i, j] = d0 + d1
    r, c = linear_sum_assignment(cost)
    out = [None] * len(plants)
    for i, j in zip(r, c):
        if cost[i, j] < big:
            out[i] = int(j)
    return out


def drift_ok(planted, reported, foff_hz, T):
    """Correct enough to tell planets apart: within 0.1 Hz/s, 10 %, or one drift bin."""
    return abs(reported - planted) <= max(0.1, 0.10 * abs(planted), abs(foff_hz) / T)


def explain_misses(plants, hits, matched, foff_hz, tsamp, nt, max_drift):
    """Reason for every plant: found / wrong-drift / suppressed / below-threshold."""
    T = (nt - 1) * tsamp
    radius = int(np.ceil(max_drift * T / abs(foff_hz)))   # turboSETI neighbour window
    used = {j for j in matched if j is not None}
    out = []
    for i, p in enumerate(plants):
        j = matched[i]
        if j is not None:
            h = hits[j]
            out.append(dict(status="found" if drift_ok(p["drift"], h["drift"], foff_hz, T)
                            else "found-wrong-drift", hit=h))
            continue
        mid = _track(p["ch"], p["drift"], T / 2, foff_hz)
        near = [(k, h) for k, h in enumerate(hits) if k not in used
                and abs(_track(h["ch"], h["drift"], T / 2, foff_hz) - mid) <= 3 + abs(p["drift"]) * T / 2 / abs(foff_hz)]
        if near:
            out.append(dict(status="wrong-path", hit=max(near, key=lambda kh: kh[1]["snr"])[1]))
            continue
        strong = [h for h in hits if abs(h["ch"] - p["ch"]) <= radius]
        if strong:
            out.append(dict(status="suppressed", hit=max(strong, key=lambda h: h["snr"])))
        else:
            out.append(dict(status="below-threshold", hit=None))
    return out


def read_dat(path, fch1_mhz, foff_mhz):
    """turboSETI .dat -> list of dict(ch, drift, snr, mhz). ch is relative to
    the file the search ran on (fch1_mhz is that file's first channel)."""
    hits = []
    for ln in open(path):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        mhz = float(p[3])
        hits.append(dict(drift=float(p[1]), snr=float(p[2]), mhz=mhz,
                         ch=(mhz - fch1_mhz) / foff_mhz))
    return hits


def selftest():
    """Offline: planting conserves power, drifts the right way for BL's negative
    foff, and the matcher refuses double credit and wrong paths."""
    rng = np.random.default_rng(1)
    nt, nc, foff, ts = 16, 8192, -2.7939677238464355, 18.253611008
    noise = rng.chisquare(64, (nt, nc)).astype(np.float32)
    d = noise.copy()
    plants = [dict(ch=2000.0, drift=+1.0, snr=50), dict(ch=6000.0, drift=-2.0, snr=50)]
    plant(d, foff, ts, plants)
    diff = d.astype(np.float64) - noise
    checks = []
    lvl = level_for(noise, 50)
    checks.append(("power conserved per spectrum",
                   np.allclose(diff[:, :4000].sum(1), lvl, rtol=0.02)))
    peak0, peak15 = np.argmax(diff[0, :4000]), np.argmax(diff[15, :4000])
    # +1 Hz/s with negative foff -> lower channel numbers, ~98 ch over 15 spectra
    checks.append(("positive drift moves to lower channels when foff<0",
                   -110 < peak15 - peak0 < -85))
    spread = int((diff[8, :4000] > lvl * 0.02).sum())
    checks.append(("drifting tone is smeared over ~6-7 channels per spectrum", 5 <= spread <= 9))
    hits = [dict(ch=2000.4, drift=+1.0, snr=48), dict(ch=2001.0, drift=+1.02, snr=30),
            dict(ch=6000.0, drift=+2.0, snr=40)]
    m = match(plants, hits, foff, ts, nt)
    checks.append(("one hit per plant, best path wins", m[0] == 0))
    checks.append(("opposite drift is NOT a match", m[1] is None))
    why = explain_misses(plants, hits, m, foff, ts, nt, max_drift=4)
    checks.append(("miss at same mid-point is labelled wrong-path or suppressed",
                   why[1]["status"] in ("wrong-path", "suppressed")))
    checks.append(("drift tolerance is not circular (0.5 Hz/s off is wrong)",
                   not drift_ok(1.0, 1.5, foff, (nt - 1) * ts)))
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(selftest() if sys.argv[1:2] in (["selftest"], []) else 2)
