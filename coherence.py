#!/usr/bin/env python3
"""coherence.py - is a narrow line a MACHINE (steady) or a MASER (amplified noise)?

Nature's narrowest radio lines are masers, and a maser is amplified noise: its
power in each spectrogram pixel fluctuates like noise does, by a fraction
1/sqrt(N_eff) of ITSELF. A transmitter's carrier is a steady oscillator: its
pixels fluctuate only because the background noise rides on top of it. With

    B      background power in the pixel      S   line power
    N_eff  independent samples averaged per pixel (measured from noise pixels)

the variance of a line's pixel power across time is

    steady oscillator (machine)  :  (B^2 + 2 S B) / N_eff
    amplified noise (maser)      :  (S + B)^2      / N_eff

coherence_index = (observed - machine) / (maser - machine)  ->  ~0 machine, ~1 maser.

Honest limits, printed with every result:
  * The test only has teeth when S is comparable to or bigger than B: the
    two predictions differ by a factor ~ (1 + S/2B). Faint lines cannot be told
    apart, and the tool says so instead of guessing.
  * RFI is coherent too: "machine" means "not a maser", never "alien".
  * Slow gain drifts or scintillation add variance; a linear trend is removed,
    and with only 16 time steps the index itself is uncertain by ~35 %.
It is a REPORTED FEATURE for candidates, not a detection gate.

REAL-DATA RESULT (2026-11-06), and why it is not a veto yet: on the one real
transmitter we have - Voyager 1's X-band carrier, S/B 3.5 - the index comes out
+3.2: MORE variable than even amplified noise. The carrier's power jumps ~3 %
rms between 18 s spectra while the band-wide gain is steady to 0.1 %, so it is
not the receiver. It is something between the source and the telescope -
most likely GBT pointing jitter on a 1.4' X-band beam, or interplanetary
scintillation. A truly steady machine therefore scored "noise-like". Until a
known-steady source in the SAME band and setup calibrates that external
modulation, an index above ~1.5 means "externally modulated - inconclusive",
never "maser". The selftest (synthetic, no external modulation) passes.

    python coherence.py selftest     # synthetic machine vs maser on real noise
    python coherence.py voyager      # the real Voyager 1 carrier + sidebands
"""
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

VOYAGER = "data/Voyager1.single_coarse.fine_res.h5"
W = 21              # channels summed around the track. Must be WIDE compared with the
                    # line's smear per spectrum (Voyager: ~2.5 ch) plus track error: with
                    # W=7 the energy spilling past the window edges jittered by ~3 %,
                    # which this test reads as maser-like self-noise.


def load(path):
    import h5py, hdf5plugin  # noqa: F401
    with h5py.File(path, "r") as f:
        a = dict(f["data"].attrs)
        d = f["data"][:, 0, :].astype(np.float64)
    freqs = a["fch1"] + a["foff"] * np.arange(d.shape[1])
    return d, freqs, a


def window_series(d, centre_ch, W=W):
    """Sum W channels around a per-spectrum centre position -> power time series."""
    nt, nc = d.shape
    out = np.empty(nt)
    for k in range(nt):
        c = int(round(centre_ch[k]))
        lo, hi = max(0, c - W // 2), min(nc, c + W // 2 + 1)
        out[k] = d[k, lo:hi].sum()
    return out


def detrend(x):
    t = np.arange(len(x))
    return x - np.polyval(np.polyfit(t, x, 1), t) + x.mean()


def n_eff_from_noise(d, avoid, W=W, n=4000, seed=0):
    """Effective sample count of a W-channel pixel, from quiet windows: mean^2/var."""
    rng = np.random.default_rng(seed)
    nt, nc = d.shape
    vals = []
    while len(vals) < n:
        c = int(rng.integers(W, nc - W))
        if any(abs(c - a) < 300 for a in avoid):
            continue
        x = d[:, c - W // 2:c + W // 2 + 1].sum(1)
        x = detrend(x)
        v = x.var(ddof=1)
        if v > 0:
            vals.append(x.mean() ** 2 / v)
    return float(np.median(vals))


def index(series, B, n_eff):
    x = detrend(series)
    P = x.mean()
    S = max(P - B, 0.0)
    v_obs = x.var(ddof=1)
    v_machine = (B ** 2 + 2 * S * B) / n_eff
    v_maser = (S + B) ** 2 / n_eff
    sep = v_maser / v_machine
    ci = (v_obs - v_machine) / (v_maser - v_machine) if v_maser > v_machine else float("nan")
    return dict(S_over_B=S / B, coherence_index=ci, separation=sep,
                verdict=("too faint to tell (S/B < 1)" if S / B < 1 else
                         "EXTERNALLY MODULATED - inconclusive (more variable than a maser)"
                         if ci > 1.5 else
                         "MACHINE-like (steady)" if ci < 0.3 else
                         "MASER-like (noise-like)" if ci > 0.7 else "ambiguous"))


def background(d, centre_ch, W=W, offsets=(-90, -60, 60, 90)):
    """Mean background power of a W-channel window, from neighbours either side."""
    return float(np.mean([window_series(d, centre_ch + o, W).mean() for o in offsets]))


def cmd_selftest():
    """Plant a steady tone and a noise-like 'maser' line into REAL Voyager noise."""
    d, freqs, a = load(VOYAGER)
    nt, nc = d.shape
    carrier = int(np.argmax(d.sum(0)))
    avoid = [carrier, carrier - 8000, carrier + 8000]
    n_eff = n_eff_from_noise(d, avoid)
    rng = np.random.default_rng(7)
    ok = True
    print(f"N_eff per {W}-channel pixel, measured from noise: {n_eff:.0f}")
    for sb in (0.3, 3.0, 30.0):
        res = {"machine": [], "maser": []}
        for trial in range(40):
            c = int(rng.integers(2000, nc - 2000))
            if abs(c - carrier) < 20000:
                continue
            centre = np.full(nt, float(c))
            B = background(d, centre)
            x = window_series(d, centre)
            S = sb * B
            # A tone mixing with noise in the receiver adds the cross term 2SB/N_eff
            # to the pixel variance - for a steady tone AND for amplified noise.
            # The real noise x already carries B^2/N_eff.
            cross = rng.normal(0.0, np.sqrt(2 * S * B / n_eff), size=nt)
            res["machine"].append(index(x + S + cross, B, n_eff)["coherence_index"])
            maser = S * rng.gamma(n_eff, 1.0 / n_eff, size=nt)   # self-noise S^2/N_eff
            cross = rng.normal(0.0, np.sqrt(2 * S * B / n_eff), size=nt)
            res["maser"].append(index(x + maser + cross, B, n_eff)["coherence_index"])
        mm, ms = np.median(res["machine"]), np.median(res["maser"])
        sep_ok = (ms - mm) > 0.5
        if sb >= 3:
            ok &= sep_ok and mm < 0.3 and ms > 0.7
        print(f"  S/B {sb:5.1f}:  machine index median {mm:+.2f}   maser index median {ms:+.2f}"
              f"   {'separated' if sep_ok else 'NOT separable (expected when faint)'}")
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_voyager():
    d, freqs, a = load(VOYAGER)
    nt, nc = d.shape
    tsamp, foff_hz = float(a["tsamp"]), float(a["foff"]) * 1e6
    avoid = []
    integ = d.sum(0)
    # Blank the spectrometer's DC spike at the centre of the coarse channel: it is
    # the brightest, steadiest, zero-drift "line" in the file and it is the
    # instrument (8419.921875 MHz here, a suspiciously round number).
    integ = integ.copy()
    integ[nc // 2 - 4:nc // 2 + 5] = np.median(integ)
    carrier = int(np.argmax(integ))
    side = int(round(22340.0 / abs(foff_hz)))      # 22.34 kHz telemetry sidebands
    lines = {"carrier": carrier, "lower sideband": carrier - side, "upper sideband": carrier + side}
    # refine each sideband to its local integrated peak
    for k in ("lower sideband", "upper sideband"):
        c = lines[k]
        lines[k] = c - 30 + int(np.argmax(integ[c - 30:c + 30]))
    avoid = list(lines.values())
    n_eff = n_eff_from_noise(d, avoid)
    drift = -0.3776                                 # Hz/s, see earth_term.py gate
    t = (np.arange(nt) + 0.5) * tsamp
    print(f"Voyager 1, N_eff per {W}-channel pixel = {n_eff:.0f} (measured from noise)")
    for name, c0 in lines.items():
        # A drifting line smears into a ~40-channel plateau when integrated, so the
        # integrated peak is NOT the track. Find the line in each spectrum near the
        # expected drift, then fit a straight track through those peaks.
        guess = c0 + drift * (t - t.mean()) / foff_hz
        span = abs(drift * t[-1] / foff_hz) + 30
        lo, hi = int(max(0, c0 - span)), int(min(nc, c0 + span))
        sm = np.array([np.convolve(d[k, lo:hi], np.ones(3), "same") for k in range(nt)])
        # restrict each spectrum's search to +-25 channels of the drift guess
        pk = []
        for k in range(nt):
            g = int(round(guess[k])) - lo
            a0, a1 = max(0, g - 25), min(hi - lo, g + 26)
            pk.append(lo + a0 + int(np.argmax(sm[k, a0:a1])))
        fit = np.polyfit(t, pk, 1)
        centre = np.polyval(fit, t)
        fitted_drift = fit[0] * foff_hz
        B = background(d, centre)
        r = index(window_series(d, centre), B, n_eff)
        print(f"  {name:15s} {freqs[int(round(centre.mean()))]:.6f} MHz  drift {fitted_drift:+.3f} Hz/s  S/B {r['S_over_B']:9.1f}  "
              f"index {r['coherence_index']:+.3f}  (maser/machine variance {r['separation']:.0f}x)  "
              f"-> {r['verdict']}")
    print("Hoped for: carrier MACHINE-like. Got: externally modulated (~3 % power jitter")
    print("between spectra; band gain steady) - see the REAL-DATA RESULT in the docstring.")
    print("'Machine' would mean 'not a maser', never 'alien'.")
    return 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    sys.exit({"selftest": cmd_selftest, "voyager": cmd_voyager}[cmd]())
