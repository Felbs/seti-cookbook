# Hunt result: Barnard's Star, the full fine-resolution ON/OFF cadence — 2026-11-06

**What this replaces.** The earlier sweep called Barnard's Star "all RFI" from a
mid-resolution file whose drift steps (~9.8 Hz/s) were too coarse to judge
drift at all (see `drift_audit.py`, `SWEEP_RESULTS.md`). This is the proper search.

**Data.** Green Bank, 2016-03-12, the six fine-resolution (2.84 Hz × 18 s)
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
  limit below holds for |drift| ≲ 0.15 Hz/s and degrades roughly as
  |drift| / 0.15 above that (≈ 6× worse at 1 Hz/s, ≈ 25× at 3.8 Hz/s).
  Earth alone gives −0.12 Hz/s here, so a transmitter whose own drift is
  small sits inside the good zone; a fast-spinning planet's would not.

With BL's SEFD of 10 Jy: S ≲ 7 Jy, an EIRP of order 3 GW at 5.96 ly, for slow
drifts — a radiometer estimate, not injection-calibrated. Not covered: the
GBT's notched 1.20–1.34 GHz radar band, |drift| > 3.78 Hz/s, one night only.
Prior work to compare with: FAST's Barnard's Star search (2023, AJ 166, 190).
