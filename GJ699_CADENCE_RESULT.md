# Hunt result: Barnard's Star, the full fine-resolution ON/OFF cadence — 2026-11-06

> **Evidence level E2 for slow drifts (≤ 0.15 Hz/s), E1 above 0.5 Hz/s** — see [EVIDENCE.md](EVIDENCE.md).

**What this replaces.** The earlier sweep called Barnard's Star "all RFI" from a
mid-resolution file whose drift steps (~9.8 Hz/s) were too coarse to judge
drift at all (see `drift_audit.py`, `SWEEP_RESULTS.md`). This is the proper search.

**Data.** Green Bank, 2016-03-10 (MJD 57457), the six fine-resolution (2.84 Hz × 18 s)
files of one ON/OFF cadence: ON, OFF, ON, OFF, ON, OFF, 288 s each, 1.1–1.9 GHz,
12.8 GB each (77 GB total).

**Search.** turboSETI with `max_drift=4` Hz/s (a planet's whole range),
`min_drift=0` (keep zero drift), SNR ≥ 10; then turboSETI's
`find_event_pipeline` (SNR cut 10). Run on a home PC at lowest priority beside
another experiment, which lost zero samples throughout.

| stage | count |
|---|---|
| hits per scan | 2,357–2,712 (15,237 total) |
| filter 1: in at least one ON | 7,630 |
| filter 2: in ONs, in no OFF | 1,029 |
| **filter 3: in all three ONs, in no OFF** | **9 events** |

**Every one of the 9 is explained** (`python cadence_verdict.py data/gj699_cadence/events_f3.csv`):

| frequency (MHz) | SNR in ON 1/2/3 | drift (Hz/s) | explanation |
|---|---|---|---|
| 1522.118110 | 19 / 16 / 14 | +0.0105 | zero drift; aeronautical telemetry band |
| 1536.971770 | 11 / 10 / 11 | +0.0105 | zero drift; Inmarsat downlink |
| 1538.888094 | 246 / 492 / 727 | +0.0105 | zero drift; Inmarsat downlink |
| 1550.955333 | 54 / 21 / 11 | +0.0105 | zero drift; Inmarsat downlink |
| 1554.131723 | 12 / 10 / 10 | +0.0105 | zero drift; Inmarsat downlink |
| 1574.980015 | 79 / 59 / 31 | +0.0105 | zero drift; GPS L1 band |
| 1575.032015 | 99 / 84 / 41 | +0.0105 | zero drift; GPS L1 band |
| 1575.860014 | 47 / 72 / 31 | +0.0105 | zero drift; GPS L1 band |
| 1576.055014 | 27 / 28 / 17 | +0.0105 | zero drift; GPS L1 band |

"Zero drift" here is decisive: any signal from Barnard's Star must carry Earth's
own −0.12 Hz/s at these frequencies (`earth_term.py`); these sit one drift step
(+0.0105) from zero — they are locked to Earth. Geostationary satellites hold
still in the sky, and Barnard's Star (Dec +4.7°) is ~10° from the geostationary
belt as seen from Green Bank, so in the ON beams they look "on target".

**The looser filter, checked too.** Of the 1,029 filter-2 hits, 1,024 sit in a
listed interference band, at zero drift, or on a natural line. Five did not:
four were single-scan hits (1426.63, 1646.66, 1769.99 MHz and one 1675 MHz
hit); the fifth was the same 1674.9977 MHz track in ON 2 and ON 3 (drift
−0.33/−0.20, phase-connected −0.30 Hz/s). Looking directly at all six files,
that track is plainly present in **both later OFF scans (SNR 31 and 33)** —
turboSETI simply had not recorded it there. Local; killed.

**Verdict: null — for slowly drifting signals.** No narrowband (≲ 3 Hz)
signal from Barnard's Star at SNR ≥ 10 in 1.1–1.9 GHz during this cadence.
Two corrections, both found by the peer-project review the same day:

- **Drift range actually searched: ±3.78 Hz/s, not ±4.** turboSETI 2.3.2 builds
  its drift blocks in steps of (N−1) integrations and stops short of
  `max_drift`; the largest drift in our own hit lists is 3.7836 Hz/s.
- **Sensitivity falls with drift.** Above ~0.15 Hz/s (one 2.8 Hz channel per
  18 s spectrum) a tone smears across channels inside each spectrum, and
  turboSETI sums one channel per spectrum. `drift_audit.py`, now planting
  realistically smeared tones, finds every slow tone (≤ 0.3 Hz/s) and misses
  every tone at ≥ 1 Hz/s at a strength that is easy without drift. So the
  limit below holds best for |drift| ≲ 0.15 Hz/s and degrades with drift
  (measured below). Earth alone gives −0.12 Hz/s here, so a transmitter whose
  own drift is small sits inside the good zone; a fast-spinning planet's would not.

**Sensitivity — measured by injection-recovery (2026-11-06, `recovery.py`).**
Smeared setigen tones were planted into one real coarse channel of ON scan 0003
(1418.6–1421.6 MHz) and searched with exactly these turboSETI settings
(max_drift 4, min_drift 0, snr 10); 1,216 plants over three runs, scored
one-to-one along each tone's track. Strength is TRUE radiometer SNR, i.e. what
the radiometer equation (and `reach.py`) means by SNR.

| drift (Hz/s) | 90 % found per scan at true SNR | S (Jy) | EIRP at 5.96 ly |
|---|---|---|---|
| ≤ 0.15 | ~70 | ~49 | ~20 GW |
| 0.5 | ≤ 200 | ≤ 140 | ≤ 56 GW |
| 1 | ~300 | ~210 | ~84 GW |
| 2 | ~600 | ~420 | ~170 GW |
| 3.5 | ~1200 | ~840 | ~340 GW |

(SEFD 10 Jy, 2.84 Hz channels, 288 s per scan, 2 polarisations, 1 Hz signal.)
Requiring the tone in all three ONs raises the slow-drift figure to about true
SNR 80: **S ≲ 56 Jy, EIRP ≲ 22 GW for slow drifts.** No plant was lost to a wrong
drift or a stronger neighbour; every miss was too faint.

**Correction (2026-11-06).** This file first quoted "S ≲ 7 Jy, ~3 GW" from the
radiometer equation at SNR 10. That is ~8× too optimistic, for two reasons found
by the injections:

1. **turboSETI's SNR is not radiometer SNR.** turboSETI 2.3.2 divides by the
   5–95 %-trimmed spread of the whole time-summed coarse channel
   (`helper_functions.comp_stats`), which includes the channel's bandpass
   shape. Measured on 62 coarse channels of ON 0003 and OFF 0004 across
   1.15–1.86 GHz, that spread is 3.2–6.1× (typically ~4.1×) the true
   channel-to-channel noise — and the true noise matches the radiometer
   equation to ~5 % (0.023 vs 0.0248 of the mean). So "SNR 10" in turboSETI
   is a true SNR of ~40.
2. **Drift and the drift tree.** Taylor-tree discretisation and smearing cost
   another ~1.7× even at slow drift. (turboSETI also never searches drift
   exactly 0 — `abs(drift_rate) > min_drift` is strict — irrelevant for a sky
   source, which always carries Earth's term.)

The null itself is unchanged. Not covered: the GBT's notched 1.20–1.34 GHz
radar band, |drift| > 3.78 Hz/s, one night only. The completeness was measured
in one coarse channel; the noise ratio, which sets most of it, was checked
band-wide.

Prior work to compare with: FAST's Barnard's Star search (2023, AJ 166, 190).
