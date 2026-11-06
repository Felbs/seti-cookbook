#!/usr/bin/env python3
"""The classic: a Doppler-drifting narrowband tone. This is what SETI has hunted
since Project Ozma (1960) and what turboSETI does properly with a Taylor tree.
Included as the reference recipe - short, obvious, and the thing every new idea
should be compared against.

Why the drift matters: the transmitter sits on a rotating, orbiting planet, so
its carrier sweeps in frequency. Earth's rotation ALONE gives ~0.05-0.3 Hz/s at
L band. A signal with EXACTLY zero drift over minutes is bolted to the ground
next to you (that is how nearly every false alarm dies), and a signal with a
drift too large to be planetary is usually a satellite.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))   # repo root
import recipe_api as R                                            # noqa: E402

NAME = "narrowband_drift"
DESCRIPTION = ("brute-force de-Doppler: shift-and-sum over trial drift rates, "
               "then peak-find. The reference technosignature recipe.")
AUTHOR = "seti-cookbook"
VERSION = "1.0"
INPUT = "spectrogram"
TAGS = ["technosignature", "narrowband", "drift", "reference"]

# drift rates a planet-borne transmitter can plausibly show, Hz/s
DEFAULTS = dict(max_drift_hz_s=4.0, zmin=10.0, ndrift=None, nmax=32)


def run(spec, params=None):
    p = dict(DEFAULTS, **(params or {}))
    d = np.asarray(spec.data, np.float64)
    nt, nch = d.shape
    res_hz = spec.res_hz
    t = spec.times_s() - spec.t0_s
    span = t[-1] if nt > 1 else 1.0

    # Trial drifts, spaced so consecutive trials differ by 1 channel of total
    # smear across the observation - the same criterion a Taylor tree uses.
    step = res_hz / max(span, 1e-9)
    n = p["ndrift"] or int(np.ceil(2 * p["max_drift_hz_s"] / step)) + 1
    n = int(min(max(n, 5), 4001))
    drifts = np.linspace(-p["max_drift_hz_s"], p["max_drift_hz_s"], n)

    best_z = np.full(nch, -np.inf)
    best_drift = np.zeros(nch)
    for r in drifts:
        shift = np.round(r * t / res_hz).astype(int)
        acc = np.zeros(nch)
        for i in range(nt):
            acc += np.roll(d[i], -shift[i])
        z = R.robust_z(acc)
        upd = z > best_z
        best_z[upd] = z[upd]
        best_drift[upd] = r

    # One drifting tone leaves a SMEAR |drift|*span/res channels wide, so a
    # naive peak-finder reports the same signal dozens of times. Suppress
    # neighbours inside the smear width of the accepted hit - one signal, one
    # candidate.
    raw, _ = R.peaks_z(best_z, zmin=p["zmin"], min_sep=2, nmax=4096)
    idx = []
    for i in raw:
        smear = max(4, int(abs(best_drift[i]) * span / res_hz) + 2)
        if all(abs(i - j) > smear for j in idx):
            idx.append(i)
        if len(idx) >= p["nmax"]:
            break

    f = spec.freqs_mhz()
    out = []
    for i in idx:
        out.append(R.Candidate(
            freq_mhz=float(f[i]), score=float(best_z[i]),
            drift_hz_s=float(best_drift[i]), t_start_s=spec.t0_s,
            duration_s=spec.duration_s, bandwidth_hz=res_hz,
            kind="techno", label="narrowband drifting tone",
            provenance=dict(method="shift-and-sum de-Doppler",
                            measures=["drift", "robust_z"],
                            n_drift_trials=n, drift_step_hz_s=round(step, 5),
                            # published so recipe_api.railed_drift_reason() can
                            # tell "drift measured" from "drift hit the grid edge"
                            drift_max_hz_s=float(p["max_drift_hz_s"]),
                            zmin=p["zmin"],
                            reference="Taylor 1974 tree de-dispersion; "
                                      "turbo_seti (Enriquez & Price 2019)")))
    return out


def selftest():
    import seti_io
    print("narrowband_drift selftest")
    ok = []
    s = seti_io.synth("drift", ntime=64, nchan=2048, f0_mhz=1400.0,
                      df_mhz=2.79e-6, drift_hz_s=-0.35, snr=20, seed=1)
    c = run(s)
    hit = c and abs(c[0].drift_hz_s + 0.35) < 0.2
    print(f"  drifting tone: {len(c)} hit(s), best drift "
          f"{c[0].drift_hz_s if c else float('nan'):+.3f} Hz/s (truth -0.350) "
          f"score {c[0].score if c else 0:.0f}")
    ok.append(bool(hit))
    for seed in (11, 12, 13):
        nz = seti_io.synth("noise", ntime=64, nchan=2048, f0_mhz=1400.0,
                           df_mhz=2.79e-6, seed=seed)
        n = run(nz)
        print(f"  pure noise seed {seed}: {len(n)} false alarm(s)")
        ok.append(not n)
    print("  RESULT:", "PASS" if all(ok) else "FAIL")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(selftest())
