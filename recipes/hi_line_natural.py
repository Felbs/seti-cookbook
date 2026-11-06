#!/usr/bin/env python3
"""Galactic hydrogen and the maser lines - the natural narrowband signals that a
technosignature search finds *instead* of aliens, and the reason a bright narrow
line is not news by itself.

Neutral hydrogen radiates at 1420.405751768 MHz (the 21 cm hyperfine spin-flip).
It is everywhere in the Milky Way, and in a Green Bank L-band file it is a BROAD
(tens of kHz to a few hundred kHz) bump whose centre frequency is a Doppler
readout of the gas velocity along your line of sight. SETI deliberately hunts
near it - Drake's "water hole", 1420 (HI) to 1720 (OH) MHz - which means SETI
data is *full* of the galaxy's own hydrogen.

The masers are worse mimics: OH at 1612.231 / 1665.402 / 1667.359 / 1720.530 MHz,
methanol at 6668.519 / 12178.593 MHz, water at 22235.08 MHz. These are genuinely
narrowband (down to a few kHz), genuinely bright, and genuinely NOT artificial -
nature's own carriers. Any recipe claiming a technosignature has to clear them.

This recipe finds them ON PURPOSE and reports the velocity, so you can tell "the
galaxy" from "a transmitter". Run it on your own BL download with
`python natural_signals.py <file>` for the full report and figure.

Line frequencies: HI from Essen et al. 1971 / IAU adopted value; OH and maser
rest frequencies per Lovas, NIST Recommended Rest Frequencies (2004).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import recipe_api as R                                            # noqa: E402
from seti_io import HI_MHZ, OH_MHZ, C_KMS                         # noqa: E402

NAME = "hi_line_natural"
DESCRIPTION = ("finds the galaxy's own 21 cm hydrogen and the OH/CH3OH/H2O maser "
               "lines, with velocities - so you can tell nature from a beacon")
AUTHOR = "seti-cookbook"
VERSION = "1.0"
INPUT = "spectrogram"
TAGS = ["natural", "spectral-line", "hydrogen", "maser", "veto"]

# (rest MHz, name, search window MHz, expected width km/s)
LINES = [
    (HI_MHZ, "HI 21 cm", 3.0, 20.0),
    (OH_MHZ[0], "OH 1612 MHz maser", 0.5, 2.0),
    (OH_MHZ[1], "OH 1665 MHz maser", 0.5, 2.0),
    (OH_MHZ[2], "OH 1667 MHz maser", 0.5, 2.0),
    (OH_MHZ[3], "OH 1720 MHz shock maser", 0.5, 2.0),
    (6668.5192, "methanol 6.7 GHz maser", 2.0, 2.0),
    (12178.593, "methanol 12.2 GHz maser", 2.0, 2.0),
    (22235.08, "water 22 GHz maser", 5.0, 3.0),
]

DEFAULTS = dict(zmin=8.0, nmax=6)


def run(spec, params=None):
    p = dict(DEFAULTS, **(params or {}))
    out = []
    norm = None
    for rest, name, win, width_kms in LINES:
        if not (spec.f_lo - win < rest < spec.f_hi + win):
            continue
        if norm is None:
            norm = spec.bandpass_normalized()
        f = norm.freqs_mhz()
        sp = norm.integrated()
        # Matched filter: a box of the line's expected width, then a LOCAL
        # baseline (median filter several line-widths wide) so the detection
        # statistic is 'excess over the local continuum' and slow instrument
        # ripple cannot masquerade as a line.
        from scipy.ndimage import uniform_filter1d, median_filter
        w_mhz = width_kms / C_KMS * rest
        wch = max(1, int(round(w_mhz / spec.df_mhz)))
        sel = np.abs(f - rest) < win
        if sel.sum() < 4 * wch:
            continue
        sm = uniform_filter1d(sp, wch, mode="nearest")
        bw = min(len(sm) - 1 | 1, max(9, 8 * wch) | 1)
        resid = sm - median_filter(sm, bw, mode="nearest")
        med = np.median(resid)
        sig = 1.4826 * np.median(np.abs(resid - med)) or (resid.std() or 1e-9)
        z = (resid - med) / sig
        med = np.median(sm)                     # for the reported line/continuum
        i = int(np.argmax(np.where(sel, z, -np.inf)))
        if z[i] < p["zmin"]:
            continue
        v = -C_KMS * (f[i] - rest) / rest
        # width at half the peak excess
        half = 0.5 * (sm[i] - med) + med
        lo = i
        while lo > 0 and sm[lo] > half:
            lo -= 1
        hi = i
        while hi < len(sm) - 1 and sm[hi] > half:
            hi += 1
        fwhm_kms = (f[hi] - f[lo]) / rest * C_KMS
        out.append(R.Candidate(
            freq_mhz=float(f[i]), score=float(z[i]),
            drift_hz_s=0.0, t_start_s=spec.t0_s, duration_s=spec.duration_s,
            bandwidth_hz=float((f[hi] - f[lo]) * 1e6),
            kind="natural",
            label=f"{name} at {v:+.1f} km/s (topocentric), FWHM {fwhm_kms:.1f} km/s",
            provenance=dict(method="bandpass-normalise + matched-filter line search",
                            measures=["line_z", "velocity_kms"],
                            line=name, rest_mhz=rest,
                            velocity_kms_topocentric=round(float(v), 2),
                            fwhm_kms=round(float(fwhm_kms), 2),
                            peak_over_continuum=round(float(sm[i] / med), 4),
                            caveat="TOPOCENTRIC velocity - no barycentric/LSR "
                                   "correction applied here; see natural_signals.py",
                            reference="HI rest 1420.405751768 MHz (IAU); "
                                      "maser rest freqs: Lovas, NIST Recommended "
                                      "Rest Frequencies (2004)")))
        if len(out) >= p["nmax"]:
            break
    return out


def selftest():
    import seti_io
    print("hi_line_natural selftest")
    ok = []
    s = seti_io.synth("hi", ntime=64, nchan=2048, f0_mhz=1419.5, df_mhz=0.0028,
                      snr=0.8, v_kms=25, seed=7)
    c = run(s)
    v = c[0].provenance["velocity_kms_topocentric"] if c else float("nan")
    print(f"  synthetic HI at +25 km/s: {len(c)} hit(s), measured {v:+.1f} km/s")
    ok.append(bool(c) and abs(v - 25) < 8)
    # CONTROL: a band with no line in it
    nz = seti_io.synth("noise", ntime=64, nchan=2048, f0_mhz=1419.5,
                       df_mhz=0.0028, seed=61)
    n = run(nz)
    print(f"  pure noise at the HI frequency: {len(n)} false alarm(s)")
    ok.append(not n)
    # CONTROL: a band nowhere near any line -> must return nothing at all
    off = seti_io.synth("drift", ntime=64, nchan=2048, f0_mhz=1300.0,
                        df_mhz=0.001, snr=50, seed=62)
    o = run(off)
    print(f"  loud carrier at 1301 MHz (no line nearby): {len(o)} hit(s)")
    ok.append(not o)
    print("  RESULT:", "PASS" if all(ok) else "FAIL")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(selftest())
