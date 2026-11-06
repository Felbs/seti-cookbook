#!/usr/bin/env python3
"""Dispersed broadband burst hunter - the recipe that finds NATURE, on purpose.

This is the FRB/pulsar signature, and it is in this cookbook for two reasons.
First, the honest one: the overwhelming majority of "signals" a SETI pipeline
turns up are natural, and you cannot claim a technosignature until you can
recognise the astrophysics (SETI_HISTORY.md part B). Second, the sneaky one:
NOVEL_DETECTORS.md #8 argues an engineered wideband impulse would look almost
exactly like this - so measuring the dispersion is also how you would spot one
that is *too* clean, or that repeats on a mathematically loud schedule.

The physics: a broadband pulse crossing the ionised interstellar medium arrives
LATER at lower frequencies, because the cold-plasma group velocity is frequency
dependent. The delay is

    t(nu) = 4.148808e3 * DM * nu_MHz^-2  seconds
          = 4.148808 ms * DM * nu_GHz^-2

with DM the dispersion measure (integrated free-electron column, pc cm^-3). In a
waterfall this is an unmistakable quadratic-in-wavelength SWEEP from high to low
frequency: milliseconds at 1.4 GHz for DM~30 (a nearby pulsar), a second or more
for DM~1000 (a cosmological FRB). Guess the DM correctly and the sweep collapses
to a vertical line; guess wrong and it stays smeared. That is all a dedispersion
search is.

Cited: Lorimer et al. 2007 (Science 318, 777) for the first FRB; Gajjar et al.
2018 (ApJ 863, 2) for Breakthrough Listen's 21 bursts from FRB 121102 at 4-8 GHz
in *public* data; CHIME/FRB Collaboration 2021 (ApJS 257, 59) for the catalogue.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import recipe_api as R                                            # noqa: E402
from seti_io import DM_CONST                                      # noqa: E402

NAME = "dispersion_sweep"
DESCRIPTION = ("dedispersion search for a broadband millisecond burst: FRBs, "
               "pulsar single pulses, and any engineered wideband impulse")
AUTHOR = "seti-cookbook"
VERSION = "1.0"
INPUT = "spectrogram"
TAGS = ["natural", "frb", "pulsar", "dispersion", "transient"]

DEFAULTS = dict(dm_max=2000.0, dm_min=0.0, ndm=None, zmin=8.0, nmax=8)


def _dm_step(spec):
    """DM spacing that smears the pulse by at most one time sample across the
    band - the standard dedispersion-trial criterion."""
    f = spec.freqs_mhz()
    span = DM_CONST * (f[0] ** -2 - f[-1] ** -2)
    return spec.dt_s / span if span > 0 else 1.0


def run(spec, params=None):
    p = dict(DEFAULTS, **(params or {}))
    if spec.ntime < 16 or spec.nchan < 8:
        return []
    step = _dm_step(spec)
    ndm = p["ndm"] or int((p["dm_max"] - p["dm_min"]) / max(step, 1e-9)) + 1
    ndm = int(min(max(ndm, 2), 4096))
    dms = np.linspace(p["dm_min"], p["dm_max"], ndm)

    best = []
    for dm in dms:
        ts = spec.dedisperse(float(dm)).timeseries()
        # the tail of the series is contaminated by the roll-in of empty rows
        f = spec.freqs_mhz()
        lost = int(np.ceil(DM_CONST * dm * (f[0] ** -2 - f[-1] ** -2) / spec.dt_s))
        good = ts[:max(8, len(ts) - lost)]
        z = R.robust_z(good)
        i = int(np.argmax(z))
        best.append((float(z[i]), float(dm), i))
    best.sort(reverse=True)

    # One burst appears at many nearby trial DMs; the highest-z trial IS the
    # burst's DM. Deduplicate by ARRIVAL TIME - one burst, one candidate.
    out = []
    used_t = []
    for z, dm, i in best:
        if z < p["zmin"]:
            break
        if any(abs(i - j) <= max(3, spec.ntime // 50) for j in used_t):
            continue
        used_t.append(i)
        f = spec.freqs_mhz()
        sweep_s = DM_CONST * dm * (f[0] ** -2 - f[-1] ** -2)
        out.append(R.Candidate(
            freq_mhz=float(0.5 * (f[0] + f[-1])), score=z, dm=dm,
            t_start_s=float(spec.t0_s + i * spec.dt_s),
            duration_s=spec.dt_s,
            bandwidth_hz=float((f[-1] - f[0]) * 1e6),
            kind="natural",
            label=f"dispersed broadband burst, DM {dm:.1f} pc/cm3 "
                  f"(sweeps {sweep_s*1e3:.1f} ms across the band)",
            provenance=dict(method="brute-force dedispersion + robust peak",
                            measures=["dm", "dedispersed_z"],
                            dm_pc_cm3=round(dm, 3),
                            dm_trials=ndm, dm_step=round(step, 4),
                            sweep_across_band_ms=round(sweep_s * 1e3, 3),
                            delay_constant="t = 4.148808e3 * DM * nu_MHz^-2 s",
                            interpretation="DM>~50 at |b|>10deg suggests "
                                           "extragalactic (FRB); DM<50 suggests "
                                           "a galactic pulsar/RRAT",
                            reference="Lorimer+2007 Science 318:777; "
                                      "Gajjar+2018 ApJ 863:2 (BL/FRB121102); "
                                      "CHIME/FRB 2021 ApJS 257:59")))
        if len(out) >= p["nmax"]:
            break
    return out


def selftest():
    import seti_io
    print("dispersion_sweep selftest")
    ok = []
    s = seti_io.synth("frb", ntime=400, nchan=512, f0_mhz=4000.0, df_mhz=4.0,
                      dt_s=0.001, dm=560, snr=35, seed=5)
    c = run(s)
    got = c[0].dm if c else float("nan")
    print(f"  synthetic FRB (truth DM 560): {len(c)} hit(s), best DM {got:.1f}, "
          f"z {c[0].score if c else 0:.0f}")
    ok.append(bool(c) and abs(got - 560) < 60)
    for seed in (51, 52):
        nz = seti_io.synth("noise", ntime=400, nchan=512, f0_mhz=4000.0,
                           df_mhz=4.0, dt_s=0.001, seed=seed)
        n = run(nz, dict(ndm=64))
        print(f"  pure noise seed {seed}: {len(n)} false alarm(s)")
        ok.append(not n)
    print("  RESULT:", "PASS" if all(ok) else "FAIL")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(selftest())
