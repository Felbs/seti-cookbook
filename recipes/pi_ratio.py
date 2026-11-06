#!/usr/bin/env python3
"""The playful one: hunt for two tones whose FREQUENCY RATIO is a famous constant.

This is the recipe that is here to make you laugh and then make you think.

The joke: everybody's first idea for "how would aliens prove they're aliens?" is
to broadcast pi. It is a cliche. Nobody runs it.

The part that is not a joke: a frequency RATIO is Doppler-invariant. Every SETI
search fights the fact that the transmitter's frequency arrives shifted and
drifting by an unknown amount that depends on two planets' orbits, rotations and
the relative velocity of two stars - so you can never say "they transmitted at
exactly X". But if a beacon puts two tones at f and pi*f, that ratio survives the
trip *unchanged*, because Doppler multiplies both tones by the same factor. Same
for the drift: both tones drift in proportion. A ratio is the one thing about a
transmitted frequency that is frame-independent, which makes it the natural
carrier for "this number is on purpose".

Nature's ratios, for contrast, are not arbitrary constants: harmonics give exact
small integers (2, 3, 4...), and the HI/OH lines give fixed atomic ratios. So the
recipe vetoes integer and known-line ratios, and only shouts for the transcendental
ones. It has found precisely nothing, which puts it in excellent company with the
entire history of SETI (SETI_HISTORY.md).

Constants tested: pi, e, the golden ratio phi, sqrt(2), and 2*pi.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import recipe_api as R                                            # noqa: E402

NAME = "pi_ratio"
DESCRIPTION = ("tone PAIRS whose frequency ratio is pi/e/phi/sqrt2 - the one "
               "property of a transmitted frequency that Doppler cannot change")
AUTHOR = "seti-cookbook (deliberately playful, genuinely Doppler-invariant)"
VERSION = "1.0"
INPUT = "spectrogram"
TAGS = ["technosignature", "playful", "doppler-invariant", "long-shot"]

CONSTANTS = {
    "pi": np.pi,
    "e": np.e,
    "phi (golden ratio)": (1 + np.sqrt(5)) / 2,
    "sqrt(2)": np.sqrt(2),
    "2*pi": 2 * np.pi,
}
DEFAULTS = dict(zmin=10.0, rel_tol=1e-4, max_tones=48, nmax=8)


def run(spec, params=None):
    p = dict(DEFAULTS, **(params or {}))
    sp = spec.integrated()
    idx, z = R.peaks_z(sp, zmin=p["zmin"], min_sep=3, nmax=p["max_tones"])
    if len(idx) < 2:
        return []
    f = spec.freqs_mhz()
    tones = sorted((float(f[i]), float(z[i])) for i in idx)
    out = []
    for a in range(len(tones)):
        for b in range(a + 1, len(tones)):
            f1, z1 = tones[a]
            f2, z2 = tones[b]
            if f1 <= 0:
                continue
            ratio = f2 / f1
            # nature's own ratios: exact small integers (harmonics of a rotator,
            # a maser ladder). Not interesting, and not a message.
            if abs(ratio - round(ratio)) < p["rel_tol"] * ratio and round(ratio) >= 2:
                continue
            for name, k in CONSTANTS.items():
                if abs(ratio - k) <= p["rel_tol"] * k:
                    out.append(R.Candidate(
                        freq_mhz=f1, score=float(min(z1, z2)),
                        drift_hz_s=0.0, t_start_s=spec.t0_s,
                        duration_s=spec.duration_s,
                        kind="techno",
                        label=f"tone pair at ratio {name} "
                              f"({f1:.6f} / {f2:.6f} MHz, ratio {ratio:.9f})",
                        provenance=dict(
                            method="peak-find both tones, test the ratio against "
                                   "Doppler-invariant constants",
                            measures=["tone_z", "frequency_ratio"],
                            f_lo_mhz=round(f1, 6), f_hi_mhz=round(f2, 6),
                            ratio=round(ratio, 9), constant=name,
                            constant_value=round(k, 9),
                            rel_error=round(abs(ratio - k) / k, 10),
                            why_this_works="Doppler multiplies both tones by the "
                                           "same (1+v/c), so a ratio is "
                                           "frame-invariant - unlike any single "
                                           "frequency",
                            integer_ratios_vetoed=True)))
                    break
            if len(out) >= p["nmax"]:
                return out
    return out


def selftest():
    import seti_io
    print("pi_ratio selftest")
    ok = []
    # A pair of tones at exactly ratio pi, inside one wideband spectrogram.
    nchan, f0, df = 4096, 1000.0, 1.0            # 1000-5096 MHz, 1 MHz channels
    s = seti_io.synth("noise", ntime=16, nchan=nchan, f0_mhz=f0, df_mhz=df, seed=71)
    f1 = 1000.0
    f2 = f1 * np.pi                              # 3141.59... MHz, in band
    for fx in (f1, f2):
        s.data[:, int(round((fx - f0) / df))] += 60.0
    c = run(s, dict(rel_tol=2e-4))
    print(f"  planted pi-ratio pair: {len(c)} hit(s)"
          + (f" -> {c[0].label}" if c else ""))
    ok.append(bool(c))
    # CONTROL: a harmonic pair (ratio exactly 3) must be vetoed as natural
    s2 = seti_io.synth("noise", ntime=16, nchan=nchan, f0_mhz=f0, df_mhz=df, seed=72)
    for fx in (1200.0, 3600.0):
        s2.data[:, int(round((fx - f0) / df))] += 60.0
    c2 = run(s2)
    print(f"  harmonic pair (ratio 3.000, control): {len(c2)} false fire(s)")
    ok.append(not c2)
    # CONTROL: pure noise
    for seed in (81, 82):
        nz = seti_io.synth("noise", ntime=16, nchan=nchan, f0_mhz=f0, df_mhz=df,
                           seed=seed)
        n = run(nz)
        print(f"  pure noise seed {seed}: {len(n)} false alarm(s)")
        ok.append(not n)
    # And the honest confession: on real data so far, nothing.
    print("  (found on real sky data so far: nothing. As with all SETI.)")
    print("  RESULT:", "PASS" if all(ok) else "FAIL")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(selftest())
