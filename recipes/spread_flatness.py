#!/usr/bin/env python3
"""Spread-spectrum / "structured noise" detector - seti-cookbook's own novel line of
attack (NOVEL_DETECTORS.md #1 and #4; cyclo.py + entropy.py are the IQ-domain
versions, validated at -10 dB SNR and on 7/7 real digital captures).

The blind spot this attacks: human radio went from loud narrowband carriers
(1900) to noise-like spread-spectrum and OFDM (5G, WiFi, GPS, military) in ONE
century. A narrowband drift search is therefore tuned to find a civilization at
exactly our 1950s stage and no other. Anyone slightly ahead of us leaks
*structured noise*, and turboSETI is structurally blind to it.

A single-dish spectrogram cannot see the cyclostationary heartbeat directly
(that needs the raw voltages - use cyclo.py for .cs16 IQ). What it CAN see is
the envelope: a WIDE band of excess power that is spectrally FLAT (noise-like,
not a pile of tones) and stays put in time. That is the spectrogram-domain
fingerprint of a spread transmitter - and it is exactly what a peak-finding
search throws away as "no signal here".

Reported candidates are honest about their ambiguity: a wide flat elevated band
is also what a broadband RFI blanket, a receiver gain step, or a continuum
source looks like. That is what the verification layer and the ON/OFF cadence
are for.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import recipe_api as R                                            # noqa: E402

NAME = "spread_flatness"
DESCRIPTION = ("wide, spectrally FLAT excess power - the spectrogram signature "
               "of a spread-spectrum/noise-like transmitter turboSETI cannot see")
AUTHOR = "seti-cookbook"
VERSION = "1.0"
INPUT = "spectrogram"
TAGS = ["technosignature", "novel", "spread-spectrum", "entropy"]

DEFAULTS = dict(min_bw_chan=32, zmin=8.0, flat_min=0.55, max_regions=8)


def run(spec, params=None):
    p = dict(DEFAULTS, **(params or {}))
    sp = spec.integrated()
    n = len(sp)
    w = int(p["min_bw_chan"])
    if n < 4 * w:
        return []
    # Robust flat baseline (the global level), then smooth on the minimum
    # bandwidth of interest: a matched filter for a WIDE plateau, not a spike.
    med = np.median(sp)
    scale = 1.4826 * np.median(np.abs(sp - med)) or (sp.std() or 1.0)
    sm = np.convolve(sp, np.ones(w) / w, "same")
    # smoothing shrinks the noise by sqrt(w); ntime already averaged in sp
    sm_sigma = scale / np.sqrt(w)
    z = (sm - med) / (sm_sigma + 1e-12)
    z[:w] = z[-w:] = 0.0                      # convolution edges are not data

    hot = z > p["zmin"]
    out = []
    i = 0
    f = spec.freqs_mhz()
    while i < n and len(out) < p["max_regions"]:
        if not hot[i]:
            i += 1
            continue
        j = i
        while j < n and hot[j]:
            j += 1
        if j - i >= w:
            seg = sp[i:j]
            flat = R.spectral_flatness(seg)
            # a comb or a single carrier inside the window is NOT spread: its
            # flatness collapses toward 0. Require noise-like flatness.
            if flat >= p["flat_min"]:
                fc = 0.5 * (f[i] + f[j - 1])
                out.append(R.Candidate(
                    freq_mhz=float(fc), score=float(z[i:j].max()),
                    drift_hz_s=0.0, t_start_s=spec.t0_s,
                    duration_s=spec.duration_s,
                    bandwidth_hz=float((f[j - 1] - f[i]) * 1e6),
                    kind="techno",
                    label=f"wide flat excess band "
                          f"{(f[j-1]-f[i])*1e6/1e3:.1f} kHz, flatness {flat:.3f}",
                    provenance=dict(
                        method="matched-filter plateau + spectral flatness",
                        measures=["plateau_z", "spectral_flatness"],
                        flatness=round(flat, 4),
                        f_lo_mhz=round(float(f[i]), 6),
                        f_hi_mhz=round(float(f[j - 1]), 6),
                        excess_over_baseline=round(float(sp[i:j].mean() / med), 3),
                        ambiguity="also matches broadband RFI, a gain step, or a "
                                  "continuum source - needs ON/OFF cadence",
                        reference="seti-cookbook cyclo.py/entropy.py "
                                  "(CYCLO_RESULT.md, ENTROPY_RESULT.md)")))
        i = j
    return out


def selftest():
    import seti_io
    print("spread_flatness selftest")
    ok = []
    s = seti_io.synth("spread", ntime=64, nchan=2048, f0_mhz=1400.0,
                      df_mhz=0.001, snr=1.4, seed=4)
    c = run(s)
    print(f"  spread band: {len(c)} hit(s)"
          + (f" -> {c[0].label}" if c else ""))
    ok.append(bool(c))
    # CONTROL 1: a single loud narrowband carrier must NOT read as spread
    car = seti_io.synth("drift", ntime=64, nchan=2048, f0_mhz=1400.0,
                        df_mhz=0.001, drift_hz_s=0.0, snr=200, seed=6)
    c2 = run(car)
    print(f"  loud narrowband carrier (control): {len(c2)} misfire(s)")
    ok.append(not c2)
    # CONTROL 2: pure noise
    for seed in (41, 42, 43):
        nz = seti_io.synth("noise", ntime=64, nchan=2048, f0_mhz=1400.0,
                           df_mhz=0.001, seed=seed)
        n = run(nz)
        print(f"  pure noise seed {seed}: {len(n)} false alarm(s)")
        ok.append(not n)
    print("  RESULT:", "PASS" if all(ok) else "FAIL")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(selftest())
