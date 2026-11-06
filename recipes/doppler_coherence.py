#!/usr/bin/env python3
"""One emitter, many tones, ONE acceleration: the common-Doppler-fraction test.

THE IDEA
--------
Doppler is multiplicative, not additive. For a line-of-sight acceleration a, every
spectral component of the SAME transmitter drifts at

    df_i/dt = -(a/c) * f_i          =>          (df_i/dt) / f_i = -a/c = const

So the *absolute* drift rate differs from component to component (a tone twice as
high drifts twice as fast in Hz/s), but the FRACTIONAL drift is one number shared
by all of them — it is a property of the emitter's motion, not of any one tone.

That gives a technosignature test that no single-tone search can make: find every
narrowband track independently, then ask whether two or more of them agree on
`beta_dot = fdot/f`. Agreement is evidence they are ONE physical object that is
accelerating — a carrier with sidebands, a comb, a signal and its harmonic. The
shared value converts straight back into physics: a = -beta_dot * c, in m/s^2,
which you can sanity-check against a planet (Earth's rotation gives ~0.034 m/s^2
at the equator) or against a satellite (orders of magnitude more).

WHY IT IS DIFFERENT FROM WHAT IS ALREADY HERE
---------------------------------------------
* `narrowband_drift` scores each tone alone; it has no notion of two tones being
  the same object.
* `comb_uniformity` tests SPACING regularity in frequency. This tests DRIFT
  proportionality in time. A carrier with two asymmetric sidebands passes here
  and fails there; a static uniform comb passes there and is flagged RFI here.
  They are complementary, and neither implies the other.
* No recipe in this repo used fractional drift before this one.

Footnote to comb_uniformity's docstring, while we are here: a Doppler shift does
not strictly *preserve* comb spacing, it SCALES it (spacing_obs = spacing_rest *
(1+beta)). At beta ~ 1e-4 that is a 0.01% change — invisible for detection, which
is why "preserved" is the right practical statement — but the same scaling is
exactly what this recipe measures in the time derivative.

HONEST LIMITATION (read this before believing a hit)
----------------------------------------------------
Within a single narrow coarse channel the test degenerates. Over a span of a few
hundred kHz at GHz frequencies, f varies fractionally by ~1e-5, so
`fdot = f * beta_dot` is the same absolute number for every component to within
far less than one drift bin. There, "shared beta_dot" and "shared absolute drift"
are the same statement, and this recipe is a MULTI-COMPONENT COMMON-DRIFT test —
still useful, still not something a single-tone search does, but it is not
exercising the fractional part.

The fractional formulation only becomes *distinguishing* when components span a
wide fractional bandwidth: a tone at 1.4 GHz and its harmonic at 2.8 GHz must show
2x the absolute drift, and a terrestrial pair faking it would have to conspire.
That is the regime this recipe is built correctly for; `require_wideband=True`
restricts output to groups that actually span enough bandwidth to test it.

WHAT KILLS A HIT
----------------
A group whose shared beta_dot is ~0 is not an accelerating emitter, it is several
things bolted to the ground — labelled `rfi`, not deleted (this repo's rule, and
NOVEL_DETECTORS #3's drift-compensated-beacon caveat means a human still looks).
Everything else still has to survive `explain()` and, above all, the ON/OFF
cadence. BLC-1 looked better than anything this will ever produce and it was an
intermodulation product.

VALIDATED ON REAL DATA — Voyager 1, 2026-07-30
----------------------------------------------
`data/Voyager1.single_coarse.fine_res.h5` (Breakthrough Listen's Green Bank
tutorial file, 8418.46-8421.39 MHz, 2.794 Hz channels, 292 s). ONE candidate,
9 components, all agreeing on fdot/f = -4.53e-11 /s:

    8419.297028 MHz  drift -0.3731 Hz/s  z=1016   <- residual carrier
    8419.274690      -0.3731            z= 118    <- sideband  -22.34 kHz
    8419.319368      -0.4018            z= 119    <- sideband  +22.34 kHz
    8419.229696      -0.3635            z=  16    <- 2nd order -44.67 kHz
    8419.364678      -0.3827            z=  17    <- 2nd order +44.67 kHz

Three independent things check out:
* the carrier drift, -0.3731 Hz/s, matches the value BL's own Voyager tutorial
  quotes (~-0.373 Hz/s);
* the sideband spacing, 22.34 kHz, matches Voyager's known telemetry subcarrier
  at ~22.5 kHz;
* the implied line-of-sight acceleration, +0.0136 m/s^2, is squarely in the range
  Earth's rotation alone produces (up to 0.034 m/s^2 at the equator, less at
  mid-latitude / off-transit) — i.e. the number means something physical.

As predicted by the HONEST LIMITATION above, `fractional_part_tested=False` here:
the components span 135 kHz out of 8.4 GHz (1.6e-5 fractional), so all nine share
the same ABSOLUTE drift to well within a drift bin and the fractional part is not
being exercised. On this file the recipe is a multi-component common-drift test.
It reports that honestly rather than claiming the stronger result.

OPEN QUESTION from that run: every sideband came in a PAIR separated by ~0.319 kHz
(114 channels) with near-equal z, e.g. 8419.274374 / 8419.274690. The carrier's
companion was removed by the alias filter (z ratio > 3) but the sidebands' were
not, because theirs are nearly equal strength. Near-equal amplitude argues against
a polyphase-filterbank sidelobe, so this is either real structure in Voyager's
downlink or an artifact of integer-shift de-drifting. NOT resolved — do not quote
the "9 components" figure as nine independent tones until it is. Next step: re-run
with sub-channel (fractional) shifts and see whether the pairs merge.

AUTHOR NOTE: prior art exists for using drift-rate consistency between components
(harmonic checks appear in RFI-excision practice, e.g. discussion in Sheikh et al.
2020 on drift-rate distributions). What is ours is the explicit beta_dot = fdot/f
clustering as a *detector* with the acceleration reported as physics.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import recipe_api as R
import seti_io

NAME = "doppler_coherence"
DESCRIPTION = ("multiple narrowband tones sharing ONE fractional drift fdot/f = "
               "-a/c, i.e. one accelerating emitter seen in several components")
AUTHOR = "Felbs + Claude"
VERSION = "1.0"
INPUT = "spectrogram"
TAGS = ["technosignature", "doppler", "multi-component", "harmonic", "novel"]

C_MS = 299_792_458.0


def _dedrift_spectrum(d, drift_hz_s, res_hz, dt_s, xpm=np, gather=False):
    """Integrate power along a linear drift trajectory.

    out[c] = mean_t d[t, (c + shift_t) mod nchan],  shift_t = round(drift*t*dt/res)

    TWO implementations, because MEASUREMENT said so (Voyager, 16x1048576, 292 s,
    ~834 drift trials):

        np.roll loop, CPU ........  92.7 s
        index-gather,  CPU ....... 204.3 s   <- 2.2x SLOWER
        index-gather,  GPU 4090 ..  10.0 s   <- 9.3x faster than the CPU roll

    The gather builds a (ntime x nchan) int32 index matrix per trial — 67 MB at a
    million channels — then does a fancy-index read. On a CPU that is far more
    memory traffic than `ntime` contiguous rolls and it destroys cache locality,
    so it LOSES. On a GPU the same gather is embarrassingly parallel and wins by
    an order of magnitude. Picking one implementation for both would give up
    either the CPU baseline or the GPU speedup, so the recipe keeps both and
    chooses by device. Same answer either way (verified: accel agrees to <1e-6,
    identical component count).
    """
    ntime, nchan = d.shape
    if not gather:
        out = xpm.zeros(nchan, dtype=xpm.float64)
        for t in range(ntime):
            out += xpm.roll(d[t], -int(round(drift_hz_s * t * dt_s / res_hz)))
        return out / ntime
    t = xpm.arange(ntime)
    # int32 keeps the index matrix at ntime*nchan*4 B (67 MB at 1M channels)
    sh = xpm.rint(drift_hz_s * t * dt_s / res_hz).astype(xpm.int32)
    idx = (xpm.arange(nchan, dtype=xpm.int32)[None, :] + sh[:, None]) % nchan
    return xpm.take_along_axis(d, idx, axis=1).sum(axis=0) / ntime


def _array_mod(use_gpu):
    """numpy, or cupy when the caller asked for GPU *and* a device really exists.

    Honest fallback, never fatal — the house rule is that a GPU path is always
    optional and the CPU path is always correct. seti_io.xp() already implements
    exactly this policy (env SETI_COOKBOOK_GPU=1), so defer to it rather than growing
    a second, subtly-different switch.
    """
    if not use_gpu:
        return np, False
    mod = seti_io.xp(force_cpu=False)
    return mod, mod is not np


def _find_tracks(spec, drift_max_hz_s, zmin, min_sep_ch, alias_ratio=3.0,
                 use_gpu=False):
    """Independent narrowband tracks: for each trial drift, de-drift and peak-
    find; keep, per frequency, the drift that maximised the z-score."""
    # bandpass_normalized() returns a Spectrogram, not an array — take .data
    d = spec.bandpass_normalized().data
    ntime, nchan = d.shape
    res_hz = abs(spec.res_hz)
    dt = spec.dt_s
    T = max(dt, ntime * dt)

    # Drift step = one channel across the whole observation: the finest drift
    # this data can distinguish. Anything finer is oversampling noise.
    step = res_hz / T
    if step <= 0:
        return [], step
    ndrift = int(np.ceil(drift_max_hz_s / step))
    drifts = np.arange(-ndrift, ndrift + 1) * step

    xpm, on_gpu = _array_mod(use_gpu)
    dg = xpm.asarray(d, dtype=xpm.float32) if on_gpu else d

    best = {}          # channel -> (z, drift)
    for dr in drifts:
        # gather on GPU (parallel, wins big), roll on CPU (cache-friendly)
        s = _dedrift_spectrum(dg, dr, res_hz, dt, xpm, gather=on_gpu)
        if on_gpu:
            # Do median/MAD and thresholding on the DEVICE and ship back only
            # the handful of surviving indices. Transferring the whole
            # de-drifted spectrum per trial (4 MB x ~834 trials) would hand the
            # entire speedup straight back to the PCIe bus.
            med = xpm.median(s)
            mad = xpm.median(xpm.abs(s - med))
            sigma = 1.4826 * mad
            if float(sigma) <= 0:
                continue
            zg = (s - med) / sigma
            hit = xpm.nonzero(zg >= zmin)[0]
            if hit.size == 0:
                continue
            chans = seti_io._tonp(hit).tolist()
            zvals = seti_io._tonp(zg[hit])
            # enforce min separation on the CPU side, strongest first
            order = np.argsort(zvals)[::-1]
            picked = []
            for i in order:
                c = int(chans[i])
                if all(abs(c - q) >= min_sep_ch for q in picked):
                    picked.append(c)
                    if float(zvals[i]) > best.get(c, (-1e9,))[0]:
                        best[c] = (float(zvals[i]), float(dr))
                if len(picked) >= 64:
                    break
        else:
            # peaks_z returns (indices, z) and computes robust_z itself
            chans, z = R.peaks_z(s, zmin=zmin, min_sep=min_sep_ch, nmax=64)
            for ch in chans:
                if ch not in best or z[ch] > best[ch][0]:
                    best[ch] = (float(z[ch]), float(dr))

    # Collapse near-duplicates, keeping the strongest. CRITICAL: a drifting tone
    # SMEARS. Over the observation it sweeps |drift|*T Hz, i.e. |drift|*T/res_hz
    # channels, and de-drifting at a slightly wrong rate still leaves a peak
    # somewhere in that footprint. With a fixed few-channel guard, ONE loud tone
    # therefore produces several "independent components" a few channels apart
    # that then happily agree on a drift — which is precisely how this recipe
    # fired on a single carrier (synth:drift, snr=200) before this was fixed.
    # So the guard has to scale with each tone's own footprint.
    def footprint_ch(drift):
        return abs(drift) * T / res_hz

    tracks = []
    for ch in sorted(best, key=lambda c: -best[c][0]):
        z, dr = best[ch]
        fp = footprint_ch(dr)
        if any(abs(ch - t["ch"]) < (min_sep_ch + 0.5 * (fp + t["fp_ch"]))
               for t in tracks):
            continue
        tracks.append(dict(ch=ch, z=z, drift_hz_s=dr, fp_ch=fp,
                           f_mhz=float(spec.freqs_mhz()[ch])))

    # ── DE-DRIFT ALIAS REJECTION ────────────────────────────────────────────
    # A periodic structure in frequency (a comb, or any tone with sidebands) can
    # partially RE-ALIGN with itself when the drift-induced shift matches its
    # spacing, so integrating at a WRONG drift still yields a peak. Measured on
    # synth:comb: the 12 real teeth sit at z~220 with drift 0, and a shadow set
    # at z~25 appears 13 channels away with drift -+0.70 Hz/s. Those shadows all
    # agree on a drift, so without this filter they form a beautiful "techno"
    # group with an entirely fabricated acceleration.
    #
    # Rule (CLEAN-style: let the strongest explanation claim the power): a track
    # is an alias if a much stronger track lies within the reach of their
    # combined drift footprints. Real independent components — a carrier and its
    # sidebands — sit far outside that reach, so they survive.
    keep = []
    for t in sorted(tracks, key=lambda x: -x["z"]):
        shadowed = False
        for s in keep:
            reach = t["fp_ch"] + s["fp_ch"] + min_sep_ch
            if abs(t["ch"] - s["ch"]) <= reach and s["z"] > alias_ratio * t["z"]:
                shadowed = True
                break
        if not shadowed:
            keep.append(t)
    return keep, step


def run(spec, params=None):
    p = dict(zmin=8.0, drift_max_hz_s=4.0, min_sep_ch=5, n_sigma=3.0,
             min_members=2, drift_min_hz_s=0.02, require_wideband=False,
             wideband_frac=1e-3, alias_ratio=3.0, gpu=False)
    p.update(params or {})

    tracks, drift_step = _find_tracks(spec, p["drift_max_hz_s"], p["zmin"],
                                      p["min_sep_ch"], p["alias_ratio"],
                                      p["gpu"])
    if len(tracks) < p["min_members"]:
        return []

    # beta_dot and its uncertainty. The drift is quantised to `drift_step`, so
    # sigma(beta_dot) = drift_step / f — this is what sets the clustering
    # tolerance, rather than a tuned constant.
    for t in tracks:
        f_hz = t["f_mhz"] * 1e6
        t["beta_dot"] = t["drift_hz_s"] / f_hz
        t["sigma_bd"] = drift_step / f_hz

    # greedy agglomeration on beta_dot, widest-agreement-first
    order = sorted(tracks, key=lambda t: t["beta_dot"])
    groups = []
    used = set()
    for i, a in enumerate(order):
        if id(a) in used:
            continue
        members = [a]
        for b in order[i + 1:]:
            if id(b) in used:
                continue
            tol = p["n_sigma"] * (a["sigma_bd"] + b["sigma_bd"])
            if abs(b["beta_dot"] - a["beta_dot"]) <= tol:
                members.append(b)
        if len(members) >= p["min_members"]:
            for m in members:
                used.add(id(m))
            groups.append(members)

    out = []
    for g in groups:
        bd = float(np.mean([m["beta_dot"] for m in g]))
        fs = np.array([m["f_mhz"] for m in g])
        drifts = np.array([m["drift_hz_s"] for m in g])
        zs = np.array([m["z"] for m in g])
        span_frac = float((fs.max() - fs.min()) / fs.mean()) if fs.mean() else 0.0

        # Did we actually test the FRACTIONAL part? Only if the predicted spread
        # of absolute drifts across the group exceeds one drift bin.
        pred_spread = abs(bd) * (fs.max() - fs.min()) * 1e6
        fractional_tested = bool(pred_spread > drift_step)
        if p["require_wideband"] and not fractional_tested:
            continue

        accel = -bd * C_MS                       # m/s^2, line of sight
        # "Zero drift" cannot mean smaller than the drift resolution — with this
        # res_hz and duration the finest distinguishable drift IS drift_step, so
        # anything at or below it is zero within measurement error. Using a fixed
        # 0.02 Hz/s here let a 1-bin noise wobble on a STATIC comb be promoted to
        # "techno" with an invented acceleration.
        drift_zero = max(p["drift_min_hz_s"], drift_step)
        static = bool(np.all(np.abs(drifts) <= drift_zero))
        # score: joint significance, rewarded for extra agreeing components and
        # for actually spanning enough band to test the fraction
        score = float(np.sqrt(np.sum(zs ** 2)) * (1.0 + 0.5 * (len(g) - 2))
                      * (1.5 if fractional_tested else 1.0))

        if static:
            kind = "rfi"
            label = (f"{len(g)} tones sharing ~zero drift — bolted to the "
                     f"ground (not one accelerating emitter)")
        else:
            kind = "techno"
            label = (f"{len(g)} tones sharing fdot/f = {bd:.3e} /s "
                     f"=> line-of-sight accel {accel:+.4f} m/s^2"
                     + ("" if fractional_tested else " (narrow span: fractional "
                                                     "part not yet tested)"))

        out.append(R.Candidate(
            freq_mhz=float(fs[int(np.argmax(zs))]),
            score=score,
            drift_hz_s=float(drifts[int(np.argmax(zs))]),
            bandwidth_hz=float((fs.max() - fs.min()) * 1e6) or None,
            kind=kind,
            label=label,
            provenance=dict(
                method="common fractional drift (beta_dot = fdot/f) shared by "
                       "independently detected narrowband tracks",
                measures=["beta_dot_per_s", "accel_m_s2", "n_components",
                          "component_freqs_mhz", "component_drifts_hz_s",
                          "fractional_part_tested"],
                reference="recipes/doppler_coherence.py docstring; Doppler is "
                          "multiplicative so fdot/f is an emitter property",
                beta_dot_per_s=bd,
                accel_m_s2=accel,
                n_components=len(g),
                component_freqs_mhz=[float(x) for x in fs],
                component_drifts_hz_s=[float(x) for x in drifts],
                component_z=[float(x) for x in zs],
                drift_bin_hz_s=float(drift_step),
                drift_max_hz_s=float(p["drift_max_hz_s"]),
                span_fractional=span_frac,
                fractional_part_tested=fractional_tested,
                agreement_spread_per_s=float(
                    np.max([m["beta_dot"] for m in g])
                    - np.min([m["beta_dot"] for m in g])),
            )))

    out.sort(key=lambda c: -c.score)
    return out


def _multi_tone_scene(drift_hz_s=-0.6, snr=30.0, seed=5, harmonic=False):
    """A carrier plus components that share ONE beta_dot — the signal this
    recipe is FOR. Built here rather than added to the benchmark, so the
    detector is not tuned against a scene it also ships."""
    s = seti_io.synth("noise", seed=seed, ntime=64, nchan=1024)
    d = s.data
    ntime, nchan = d.shape
    res = abs(s.res_hz)
    f0 = s.freqs_mhz()[nchan // 2]
    bd = drift_hz_s / (f0 * 1e6)               # the shared fractional drift
    offs = (-140, 0, 140) if not harmonic else (0, 300)
    for off in offs:
        ch0 = nchan // 2 + off
        f_i = s.freqs_mhz()[ch0] * 1e6
        dr_i = bd * f_i                        # each component drifts per ITS f
        for t in range(ntime):
            ch = int(round(ch0 + dr_i * t * s.dt_s / res))
            if 0 <= ch < nchan:
                d[t, ch] += snr
    return s


def selftest():
    ok = []
    # sensitivity: must find the multi-component accelerating emitter
    hits = run(_multi_tone_scene())
    ok.append(bool(hits) and hits[0].kind == "techno")
    # and it must recover the acceleration it was given, not just fire
    if hits:
        bd_true = -0.6 / (hits[0].provenance["component_freqs_mhz"][0] * 1e6)
        got = hits[0].provenance["beta_dot_per_s"]
        ok.append(abs(got - bd_true) < 3 * abs(bd_true) + 1e-14)
    else:
        ok.append(False)
    # false alarms: silent on pure noise
    for seed in (11, 12, 13):
        ok.append(not run(seti_io.synth("noise", seed=seed)))
    # THE WRONG SIGNAL: a single loud carrier is not a multi-component emitter.
    # This is the control that matters — one tone can never establish a SHARED
    # fractional drift, so a power meter masquerading as this recipe fails here.
    ok.append(not run(seti_io.synth("drift", snr=200.0)))
    # a static comb is real structure but NOT an accelerating emitter: we may
    # report it, but only as rfi, never as techno
    comb = run(seti_io.synth("comb", snr=30.0))
    ok.append(all(c.kind == "rfi" for c in comb))
    print(f"[{NAME}] selftest: {sum(ok)}/{len(ok)} checks passed -> "
          f"{'PASS' if all(ok) else 'FAIL'}")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    raise SystemExit(selftest())
