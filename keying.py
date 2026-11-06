#!/usr/bin/env python3
"""keying.py - find carriers that blink: on, off, on - the way you send bits.

turboSETI adds up a drift track over the whole 5-minute scan. A carrier that is
ON for only 5 of the 16 spectra loses most of its signal-to-noise in that sum,
and can fall under the threshold even when each "on" spectrum is bright. On/off
keying is the oldest way to send information, so this is a blind spot worth
closing.

  stage 1   find bright pixels in EACH 18 s spectrum on its own (after removing
            the bandpass)
  stage 2   link them: a Hough vote over drift rate (+-4 Hz/s) - pixels from
            different spectra that fall on ONE straight drifting line
  stage 3   describe the track: how many spectra it is on in (duty cycle), the
            on/off pattern, and the SNR turboSETI would have seen integrating it
The false-alarm level for stage 2 is MEASURED by re-running with the spectra
shuffled in time (that destroys real tracks and keeps everything else).

Detection rule (fixed in selftest before any real-sky run): a track needs
>= LINK_MIN linked spectra and must beat the shuffled-time null.

    python keying.py selftest      # plant keyed carriers into REAL Green Bank noise
    python keying.py scan          # the GJ699 HI coarse channel, all 6 scans (on disk)
    python keying.py scan FILE.h5  # any BL fine-resolution file (one coarse channel at a time)
"""
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent

PIX_SIGMA = 6.0          # stage-1 threshold per pixel, in robust sigma
LINK_MIN = 3             # stage-2: a track must be seen in at least this many spectra
DRIFT_MAX = 4.0          # Hz/s


def flatten(x):
    """Divide out the bandpass (time-median per channel), robust-normalise each
    spectrum. Returns z-scores."""
    # Bandpass from the time-AVERAGED spectrum smoothed ACROSS frequency (51 ch).
    # A per-channel median over time would swallow any carrier that is on most of
    # the scan - exactly the steady signals we must not erase.
    from scipy.ndimage import median_filter
    bp = median_filter(x.mean(axis=0), size=51, mode="nearest")
    bp = np.where(bp > 0, bp, np.median(bp))
    y = x / bp - 1.0
    # A drifting tone splits its power between two neighbouring channels; sum
    # adjacent pairs so a split tone is seen whole. Pixel j now means j..j+1.
    y = y + np.concatenate([y[:, 1:], y[:, -1:]], axis=1)
    med = np.median(y, axis=1, keepdims=True)
    mad = np.median(np.abs(y - med), axis=1, keepdims=True) * 1.4826
    return (y - med) / np.where(mad > 0, mad, 1.0)


def peaks(z, sigma=PIX_SIGMA):
    """(spectrum index, channel, z) of every pixel above sigma that is a local max."""
    t, c = np.nonzero(z > sigma)
    keep = []
    for i, j in zip(t, c):
        lo, hi = max(0, j - 2), min(z.shape[1], j + 3)
        if z[i, j] >= z[i, lo:hi].max():
            keep.append((i, j, z[i, j]))
    return np.array(keep, dtype=float).reshape(-1, 3)


def link(pk, nt, tsamp, foff_hz, drift_max=DRIFT_MAX):
    """Hough vote: for each drift rate, project every peak back to its channel at
    the scan middle; count DISTINCT spectra landing in the same 1-channel bin."""
    if len(pk) == 0:
        return []
    t = (pk[:, 0] + 0.5 - nt / 2) * tsamp
    step = abs(foff_hz) / (nt * tsamp)                    # 1 channel over the scan
    drifts = np.arange(-drift_max, drift_max + step / 2, step)
    found = {}
    # 2-channel bins on two staggered grids, so peaks that round to neighbouring
    # channels still vote together
    for dr, off in ((d, o) for d in drifts for o in (0, 1)):
        c0 = (np.floor((pk[:, 1] - dr * t / foff_hz + off) / 2).astype(np.int64)) * 2 - off
        order = np.argsort(c0, kind="stable")
        cs, ts = c0[order], pk[order, 0].astype(int)
        bounds = np.flatnonzero(np.diff(cs)) + 1
        for grp_c, grp_t, grp_i in zip(np.split(cs, bounds), np.split(ts, bounds),
                                       np.split(order, bounds)):
            n = len(set(grp_t.tolist()))
            if n >= LINK_MIN:
                key = int(grp_c[0])
                if key not in found or n > found[key]["n_on"]:
                    found[key] = dict(c_mid=key, drift=float(dr), n_on=n,
                                      spectra=sorted(set(grp_t.tolist())),
                                      z=[float(v) for v in pk[grp_i, 2]])
    # merge neighbours (the same track seen at adjacent drift/channel bins)
    out = []
    for k in sorted(found):
        if out and k - out[-1]["c_mid"] <= 3 and abs(found[k]["drift"] - out[-1]["drift"]) < 0.2:
            if found[k]["n_on"] > out[-1]["n_on"]:
                out[-1] = found[k]
            continue
        out.append(found[k])
    return out


def integrated_snr(z, tr, tsamp, foff_hz):
    """The SNR a whole-scan drift integration (turboSETI-style) would give this
    track: sum the 16 z-values along it (+-1 channel), divide by sqrt(16)."""
    nt = z.shape[0]
    t = (np.arange(nt) + 0.5 - nt / 2) * tsamp
    tot = 0.0
    for k in range(nt):
        c = int(np.floor(tr["c_mid"] + tr["drift"] * t[k] / foff_hz))   # pair c..c+1
        if 0 <= c < z.shape[1]:
            tot += z[k, c]
    return tot / np.sqrt(nt)


def describe(z, tracks, tsamp, foff_hz):
    nt = z.shape[0]
    for tr in tracks:
        tr["duty"] = tr["n_on"] / nt
        tr["pattern"] = "".join("#" if k in tr["spectra"] else "." for k in range(nt))
        tr["int_snr"] = integrated_snr(z, tr, tsamp, foff_hz)
        tr["kind"] = ("steady (turboSETI sees it)" if tr["duty"] > 0.85 else
                      "KEYED - turboSETI would miss it" if tr["int_snr"] < 10 else
                      "keyed, bright enough for turboSETI too")
    return tracks


def null_rate(z, tsamp, foff_hz, n_perm=3, seed=0):
    """Tracks found with the spectra SHUFFLED in time: the chance-alignment rate."""
    rng = np.random.default_rng(seed)
    counts = []
    for _ in range(n_perm):
        zp = z[rng.permutation(z.shape[0])]
        counts.append(len(link(peaks(zp), z.shape[0], tsamp, foff_hz)))
    return float(np.mean(counts))


def load_npz(path):
    zf = np.load(path)
    return zf["data"].astype(np.float64), float(zf["tsamp"]), float(zf["foff"]) * 1e6, \
        float(zf["fch1"]) + float(zf["foff"]) * int(zf["ch0"])


def plant(x, c_mid, drift, pattern, snr_pix, tsamp, foff_hz):
    """Add a carrier ON in the spectra marked '#', per-spectrum SNR snr_pix."""
    y = x.copy()
    nt = x.shape[0]
    t = (np.arange(nt) + 0.5 - nt / 2) * tsamp
    from scipy.ndimage import median_filter
    bp = median_filter(x.mean(axis=0), size=51, mode="nearest")
    r = x / bp - 1
    pair = r + np.concatenate([r[:, 1:], r[:, -1:]], axis=1)
    rel = np.median(np.abs(pair - np.median(pair, axis=1, keepdims=True)), axis=1) * 1.4826
    for k in range(nt):
        if pattern[k] != "#":
            continue
        pos = c_mid + drift * t[k] / foff_hz
        j = int(np.floor(pos))
        fr = pos - j
        a = snr_pix * rel[k] * bp[j]
        # spread like a real drifting tone: within one 18 s spectrum the tone moves
        # drift*tsamp Hz, but at |drift| < 0.16 Hz/s that is < 1 channel
        y[k, j] += a * (1 - fr)
        y[k, j + 1] += a * fr
    return y


def cmd_selftest():
    npz = sorted((HERE / "data").glob("hi_*GJ699_0003*.npz"))
    if not npz:
        sys.exit("needs data/hi_*.npz from `python hi_starframe.py fetch` (real GBT noise)")
    x, tsamp, foff_hz, _ = load_npz(npz[0])
    nt, nc = x.shape
    rng = np.random.default_rng(5)
    # integrated SNR = n_on * snr / sqrt(16); turboSETI (threshold 10) misses it when n_on*snr < 40
    cases = [  # (name, pattern, per-spectrum SNR, drift, expect keying?, expect turboSETI?)
        ("keyed 5 of 16, SNR 7.5 each", "#..#.#...#...#..", 7.5, 0.11, True, False),
        ("keyed 3 of 16, SNR 10 each", "..#.....#....#..", 10.0, -0.35, True, False),
        ("keyed 6 of 16, SNR 9 each", "#..#.#..##...#..", 9.0, 0.21, True, True),
        ("steady carrier, SNR 9 each", "#" * 16, 9.0, 0.05, True, True),
        ("faint steady, SNR 3 each", "#" * 16, 3.0, -0.2, False, True),
        ("single 18 s blip, SNR 20", "......#.........", 20.0, 0.0, False, False),
    ]
    ok = True
    print(f"real GBT noise: {npz[0].name[3:60]}...  ({nt} x {nc})")
    for name, pat, snr, dr, want_k, want_t in cases:
        c_mid = int(rng.integers(50_000, nc - 50_000))
        y = plant(x, c_mid, dr, pat, snr, tsamp, foff_hz)
        sl = slice(c_mid - 4000, c_mid + 4000)              # a local patch is enough
        z = flatten(y[:, sl])
        tr = describe(z, link(peaks(z), nt, tsamp, foff_hz), tsamp, foff_hz)
        hit = [t for t in tr if abs(t["c_mid"] - 4000) <= 3]
        got_k = bool(hit)
        isnr = integrated_snr(z, dict(c_mid=4000, drift=dr), tsamp, foff_hz)
        got_t = isnr >= 10
        good = (got_k == want_k) and (got_t == want_t)
        ok &= good
        desc = (f"found {hit[0]['n_on']}/16 '{hit[0]['pattern']}' drift {hit[0]['drift']:+.3f}"
                if hit else "not linked")
        print(f"  {name:28s} keying: {desc:46s} turboSETI-style SNR {isnr:5.1f}  "
              f"{'OK' if good else 'UNEXPECTED'}")
    z = flatten(x[:, 200_000:400_000])
    n_null = null_rate(z, tsamp, foff_hz)
    n_real = len(link(peaks(z), nt, tsamp, foff_hz))
    print(f"  false-alarm check on 200k real channels: {n_real} tracks in time order, "
          f"{n_null:.1f} with spectra shuffled")
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def scan_one(p):
    x, tsamp, foff_hz, f0 = load_npz(p)
    nt, nc = x.shape
    edge, dc = 2000, nc // 2                               # coarse roll-off, DC spike
    z = flatten(x)
    z[:, :edge] = 0
    z[:, -edge:] = 0
    z[:, dc - 3:dc + 4] = 0
    tr = describe(z, link(peaks(z), nt, tsamp, foff_hz), tsamp, foff_hz)
    for t in tr:
        t["f_mhz"] = f0 + t["c_mid"] * foff_hz / 1e6
    return tr, null_rate(z, tsamp, foff_hz, n_perm=2)


def cmd_scan(paths):
    """Cadence logic (fixed here, before reading the result): a candidate must be
    linked in ALL THREE ON scans within +-1 kHz (Earth's drift between scans is
    < 100 Hz at L band), with nothing linked within +-10 kHz in ANY OFF scan."""
    if not paths:
        paths = sorted((HERE / "data").glob("hi_*.npz"))
    per = []
    for p in paths:
        kind = "OFF" if "_OFF_" in p.name else "ON"
        tr, nn = scan_one(p)
        keyed = [t for t in tr if t["kind"].startswith("KEYED")]
        per.append((kind, p.name, tr))
        print(f"{kind:3s} {p.name[3:62]}: {len(tr):3d} tracks in time order vs {nn:5.1f} with "
              f"spectra shuffled; {len(keyed)} keyed below turboSETI's threshold")
    ons = [tr for k, _, tr in per if k == "ON"]
    offs = [t for k, _, tr in per if k == "OFF" for t in tr]
    cands = []
    for t in (ons[0] if ons else []):
        if any(abs(o["f_mhz"] - t["f_mhz"]) * 1e6 < 10e3 for o in offs):
            continue
        if all(any(abs(u["f_mhz"] - t["f_mhz"]) * 1e6 < 1e3 for u in other) for other in ons[1:]):
            cands.append(t)
    rej = sum(1 for t in (ons[0] if ons else [])
              if any(abs(o["f_mhz"] - t["f_mhz"]) * 1e6 < 10e3 for o in offs))
    print(f"\ncadence: {len(ons[0]) if ons else 0} tracks in ON 1; {rej} within 10 kHz of an OFF "
          f"track (local); ON-only and in all 3 ONs: {len(cands)}")
    for t in cands:
        print(f"   CANDIDATE {t['f_mhz']:.6f} MHz drift {t['drift']:+.3f} {t['pattern']} "
              f"int SNR {t['int_snr']:.1f}  -> rung 2, try to kill it")
    if not cands:
        print("   none - NULL for keyed carriers in this coarse channel")
    return 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    sys.exit(cmd_selftest() if cmd == "selftest" else cmd_scan([Path(a) for a in sys.argv[2:]]))
