# Cross-epoch matching and epoch stacking in the star frame — 2026-11-06

> Results below the pre-registration are added after the runs. This section
> was committed BEFORE any stack was computed on real data.

## Pre-registration: star-frame stack (`stack.py`)

Fixed before looking (values in `stack.PREREG`):

| item | value |
|---|---|
| beacon model | steady tone at the hydrogen line, 1420.405751768 MHz, in the star's rest frame |
| window | star frame, 1420.405751768 MHz ± 3.0 kHz (RV error ±0.5 km/s = ±2.4 kHz, plus margin) |
| star model | targets.json (SIMBAD): position propagated with proper motion and parallax; RV(t) including the secular term |
| intrinsic drift grid | −0.05 … +0.05 Hz/s, step 0.01 (11 rows), on top of the predicted Earth drift at mid-scan |
| RFI mask | channels where an OFF of the same epoch exceeds z = 6 in any drift row, ±3 channels; DC ±1 channel |
| statistic | max over the window and the drift grid of the stacked z |
| threshold | max(6.0, 99th percentile of the same statistic at ≥ 300 empty windows of the same stack) |
| detection | above threshold in the stack AND not above it at that place in the OFF-as-ON control stack |
| gain | measured from setigen plants (Doppler-smeared, true radiometer SNR per scan) at 2–12 per scan in every epoch |
| false alarms | OFF-as-ON control: the same pipeline with the OFF scans playing the ONs, fraction of empty windows above threshold |

**Run 1 — GJ699 PSEUDO-epochs.** The archive holds only ONE L-band date of
Barnard's Star (2016-03-10; see the census below), so the three ON scans of
that session are stacked as three pseudo-epochs: P1 = ON 0003 (OFF 0004 for
its mask), P2 = ON 0005 (OFF 0006), P3 = ON 0007 (OFF 0008), each at its own
real mid-scan time. This tests the machinery on real noise and real RFI; it
is NOT a multi-date result. Data: the clean 749,568-channel part of coarse
channel 155 of each `data/hi_*GJ699*.npz` slice (KNOWN_ISSUES #8).

**Run 2 — real dates, if fetched in time.** GJ725 (Struve 2398, 11.5 ly),
L band, sessions MJD 57587, 57667 and 57804 (2016-07-26, 2016-10-14,
2017-03-01), the coarse channel that holds the star-frame hydrogen line, all
three ON scans of each session; same window rule, threshold rule and drift
grid. GJ725 A and B share one archive name and one pointing; the star model
uses the first targets.json row (`Star.from_targets("GJ725")`).
