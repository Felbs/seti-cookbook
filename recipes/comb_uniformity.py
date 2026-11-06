#!/usr/bin/env python3
"""Frequency-comb detector, spectrogram edition - one of seti-cookbook's own novel
detectors (see NOVEL_DETECTORS.md #2 and COMB_RESULT.md; comb.py is the IQ
version validated at 3/3 combs, 0/3 false fires on randomly-spaced tones).

The idea: nature makes tones, but it does not make tones at PERFECTLY UNIFORM
Hz spacing across a wide band. Molecular rotational lines are uniform in
*quantum number*, not in Hz; harmonics of a rotator grow multiplicatively;
plasma lines are not a ladder. A uniform Hz-spaced comb is an engineered
reference - and it survives Doppler shift, because the whole comb slides while
the SPACING is preserved. That makes it a message that needs no message.

How: a comb spaced df makes the power SPECTRUM ITSELF periodic with period df,
so the spectrum-of-the-spectrum (cepstrum) concentrates at one quefrency. High
cepstral peak/median = genuinely periodic; a pile of randomly-spaced tones
scatters it flat. This is a UNIFORMITY test, not a tone counter.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import recipe_api as R                                            # noqa: E402

NAME = "comb_uniformity"
DESCRIPTION = ("cepstral test for a UNIFORMLY Hz-spaced tone comb - an "
               "engineered frequency reference nature does not fake")
AUTHOR = "seti-cookbook"
VERSION = "1.0"
INPUT = "spectrogram"
TAGS = ["technosignature", "novel", "comb", "engineered-uniformity"]

DEFAULTS = dict(sig_min=14.0, qmin=8, min_teeth=5, tooth_z=6.0)


def run(spec, params=None):
    p = dict(DEFAULTS, **(params or {}))
    sp = spec.integrated()
    n = len(sp)
    if n < 64:
        return []
    x = sp - sp.mean()
    Q = np.abs(np.fft.rfft(x))
    band = Q[p["qmin"]:]
    if len(band) < 8:
        return []
    med = np.median(band) + 1e-12
    sig = float(band.max() / med)
    if sig < p["sig_min"]:
        return []

    # spacing from the spectrum's autocorrelation fundamental (free of the
    # harmonic ambiguity that fools a cepstral argmax)
    ac = np.fft.irfft(Q ** 2, n=2 * (len(Q) - 1))[:n // 2]
    lmin = max(4, p["qmin"])
    D = int(np.argmax(ac[lmin:])) + lmin
    spacing_hz = D * spec.res_hz

    # count real teeth so we can report an honest tooth count, and anchor the
    # candidate on the lowest one
    idx, z = R.peaks_z(sp, zmin=p["tooth_z"], min_sep=max(2, D // 4), nmax=256)
    teeth = sorted(idx)
    if len(teeth) < p["min_teeth"]:
        return []
    f = spec.freqs_mhz()
    # how uniform are the observed gaps, really?
    gaps = np.diff([f[i] for i in teeth]) * 1e6
    uni = float(np.std(gaps) / (np.mean(gaps) + 1e-12)) if len(gaps) > 1 else 1.0
    return [R.Candidate(
        freq_mhz=float(f[teeth[0]]), score=sig,
        drift_hz_s=0.0, t_start_s=spec.t0_s, duration_s=spec.duration_s,
        bandwidth_hz=float((f[teeth[-1]] - f[teeth[0]]) * 1e6),
        kind="techno", label=f"frequency comb, {len(teeth)} teeth @ "
                             f"{spacing_hz:.1f} Hz spacing",
        provenance=dict(method="cepstral uniformity + spectrum autocorrelation",
                        measures=["cepstral_peak_over_median", "gap_uniformity"],
                        cepstral_significance=round(sig, 2),
                        spacing_hz=round(spacing_hz, 3),
                        n_teeth=len(teeth),
                        gap_scatter_fraction=round(uni, 4),
                        tooth_freqs_mhz=[round(float(f[i]), 6) for i in teeth[:24]],
                        reference="seti-cookbook comb.py / COMB_RESULT.md "
                                  "(NOVEL_DETECTORS.md #2)"))]


def selftest():
    import seti_io
    print("comb_uniformity selftest")
    ok = []
    s = seti_io.synth("comb", ntime=32, nchan=2048, f0_mhz=1400.0, df_mhz=0.001,
                      teeth=12, snr=14, seed=3)
    c = run(s)
    print(f"  uniform comb: {len(c)} hit(s)"
          + (f" -> {c[0].label}, sig {c[0].score:.0f}" if c else ""))
    ok.append(bool(c))
    # CONTROL: the same number of tones, RANDOM spacing -> must stay silent
    rng = np.random.default_rng(5)
    base = seti_io.synth("noise", ntime=32, nchan=2048, f0_mhz=1400.0, df_mhz=0.001,
                         seed=21)
    for ch in rng.choice(np.arange(20, 2020), 12, replace=False):
        base.data[:, ch] += 14.0
    c2 = run(base)
    print(f"  12 RANDOM-spaced tones (control): {len(c2)} false fire(s)")
    ok.append(not c2)
    for seed in (31, 32):
        nz = seti_io.synth("noise", ntime=32, nchan=2048, f0_mhz=1400.0,
                           df_mhz=0.001, seed=seed)
        n = run(nz)
        print(f"  pure noise seed {seed}: {len(n)} false alarm(s)")
        ok.append(not n)
    print("  RESULT:", "PASS" if all(ok) else "FAIL")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(selftest())
