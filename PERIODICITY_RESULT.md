# Periodic pulse trains: Barnard's Star mid-resolution cadence — 2026-11-07

> **Evidence level: tool E2 (completeness measured by plants in real data, three bands of
> one scan, two search codes); the Barnard's Star search is E1** (full ON/OFF cadence,
> every ON-only event explained or listed below, but completeness was measured in 48 of
> 256 coarse channels of one scan only). See [EVIDENCE.md](EVIDENCE.md).

## What it is

`periodicity.py` looks for a pulse train in one frequency channel: on for a fraction d of
every period P. For each channel of a dynamic spectrum it folds the time series at every
trial period of a pre-registered range and scores a boxcar matched filter over pulse
phase and width (method from BLIPSS, Suresh et al. 2023, AJ 165, 255, MIT; reimplemented in
numpy, no code copied). Steps: divide each spectrum by its median over its own coarse
channel (common gain); per channel subtract the median and divide by the robust spread;
blank DC channels; fold; score `snr = window sum / sqrt(n (1 − n/N))` (mean 0, spread 1 on
white noise at any fixed trial: selftest −0.003 / 0.997); per channel merge neighbouring
periods (friends-of-friends), keep the best per cluster, label harmonics (P ≈ p0/n) and
sub-harmonics (P ≈ n·p0); merge neighbouring channels with the same period (wider than
4 channels = BROAD). Across the cadence every ON fundamental ≥ 8 is **measured** in every
scan at its channel ±1 and period ± tolerance (a forced fold, not "was a candidate logged",
KNOWN_ISSUES #9), giving a code such as `101010` (1 = forced snr ≥ 8 in an ON, ≥ 6 in an
OFF: asymmetric bars as in BLIPSS, an OFF vetoes more easily).

## Which data can host it

A period needs ≥ 9 samples per period (bins) and ≥ 3 periods per scan:
`tsamp × 9 ≤ P ≤ N tsamp / 3` (`python periodicity.py reach`).

| product (GBT, 2016) | samples per scan | reachable periods |
|---|---|---|
| fine-res .0000 (2.84 Hz × 18 s) | 16 | **none** (162 s > 96 s; refused) |
| mid-res .0002 (2.86 kHz × 1.058 s) | 277 | **9.5 – 97.7 s** |
| high-time-res .0001 (366 kHz × 349 µs) | ~824,000 | 3 ms – 96 s, but ~1.3 GB per scan |

Pre-registered (in the code, `PREREG`, before the first real search): mid-res .0002,
periods 10–70 s (BLIPSS's range), ≥ 3 periods per scan, 9–10 phase bins, duty cycle up to
0.5 (boxcar 1–5 bins: duty resolution ~0.1), ON bar 8, OFF bar 6, candidates logged
from 6. Period step 0.5·p²/(N·nb): **502 trial periods, 24,694 trials per channel** (periods
× phases × widths). The six GJ699 .0002 files (2016-03-10, 200 MB each, 1.2 GB) were
downloaded whole.

## Injection gate (real data)

Boxcar pulse trains planted into the real ON scan 0003, one per channel in channels with
no candidate of their own, at each channel's **measured** robust noise; duty 0.1, period
log-uniform 10–70 s, random phase; searched with the production settings. Strength is the
**integrated** matched-filter SNR (on-pulse excess over the noise of n on-pulse samples);
BLIPSS's "pulse SNR" (per-sample peak) is in brackets. 40 plants per level per band.
`results/periodic_inject/`.

| integrated SNR (pulse SNR) | 1395–1442 MHz (cc 150–165) | 1744–1788 MHz (cc 30–45) | 1504–1548 MHz (cc 112–127, Inmarsat/AMT) | pooled | riptide-ffa, cc 150–165 |
|---|---|---|---|---|---|
| 6 (1.2) | 0/40 | 0/40 | 0/40 | 0.00 | 0/40 |
| 7 (1.4) | 0/40 | 1/40 | 0/40 | 0.01 | 0/40 |
| 8 (1.6) | 3/40 | 3/40 | 2/40 | 0.07 | 0/40 |
| 9 (1.8) | 7/40 | 8/40 | 9/40 | 0.20 | 0/40 |
| 10 (2.0) | 17/40 | 18/40 | 21/40 | 0.47 | 0/40 |
| 12 (2.4) | 36/40 | 36/40 | 35/40 | 0.89 | 6/40 |
| 15 (3.0) | 40/40 | 40/40 | 40/40 | 1.00 | 13/40 (+1 harmonic) |
| 20 (4.0) | 40/40 | 40/40 | 40/40 | 1.00 | 33/40 (+1) |
| 30 (6.0) | 40/40 | 40/40 | 40/40 | 1.00 | 40/40 |
| 100 (20, the BLIPSS gate) | 40/40 | 40/40 | 40/40 | 1.00 | 39/40 (+1) |

**90 % recovery at integrated SNR 12.2** (12.0 / 12.0 / 12.5 in the three bands; log
interpolation). BLIPSS's own gate, pulse SNR 20 at duty 0.1, is integrated SNR ~100 here:
120 of 120 found. The numpy statistic loses ~23 % to phase-bin and period-grid
quantisation (selftest median 0.77 of the planted SNR), so 50 % recovery sits at ~10, not 8.

**Second search code.** riptide-ffa 0.2.7 (`pip install riptide-ffa` into a separate virtualenv;
installed and ran on a 4-core CPU without AVX2) on the same planted data at the same bar of 8 on
its own S/N: 90 % at integrated SNR 23.8, ~2× ours: riptide's S/N runs below ours on the
same trains (probably because it renormalises the series with the pulses in it and defines
its boxcar S/N differently; not traced in its code), so at a common bar it is less sensitive; with 0 unplanted candidates ≥ 8 in
this band in both codes, the two disagree only in scale. riptide's detections matched the
planted period except 3 found at a harmonic. riptide is ~5× faster (0.6 ms per channel).

## Null rate

Per scan: 261,888 channels × 24,694 trials = **6.5 × 10⁹ trials**. For independent
Gaussian trials one would expect ≤ 6 above 6 and ~4 × 10⁻⁶ above 8 per scan (the trials
are correlated, so these are upper bounds). Measured (fundamental groups ≥ 8 per scan, all
six real scans, plus scan 0004 with each block's spectra shuffled in time):

| data | groups ≥ 8 | in a listed RFI band | 1626.5–1648 MHz (next to Iridium, not listed) | elsewhere |
|---|---|---|---|---|
| ON 0003 | 3,899 | 3,778 | 121 | **0** |
| OFF 0004 | 4,107 | 3,970 | 136 | **1** (1717 MHz) |
| ON 0005 | 4,170 | 4,078 | 92 | **0** |
| OFF 0006 | 4,217 | 4,133 | 84 | **0** |
| ON 0007 | 4,042 | 3,941 | 101 | **0** |
| OFF 0008 | 3,511 | 3,471 | 40 | **0** |
| OFF 0004 **time-shuffled** | 4,531 | 4,351 | 178 | 2 (incl. 1717 MHz) |

80–85 % of them are in the Iridium allocation (1610–1626.5 MHz), the rest in the
Inmarsat, GOES, GNSS, radar and AMT bands. **Time-shuffling the spectra leaves the count
unchanged**: these channels are impulsive (bursty), and a few huge samples give a large
fold S/N at any period. So in RFI bands this statistic measures burstiness, not
periodicity (KNOWN_ISSUES #13). In the rest of the band, the null is what Gaussian noise
predicts: 1 group ≥ 8 in six scans × ~80,000 quiet channels (~225 MHz), and that one is also in the
shuffled data. The injection gate's three bands (cc 150–165, 30–45 quiet; cc 112–127 with
879 unplanted groups ≥ 8) all gave the same 90 % point, because plants went into channels
without candidates of their own.

## Barnard's Star cadence result

Six mid-res scans, ON OFF ON OFF ON OFF (2016-03-10, 1.1–1.9 GHz). Every ON
fundamental ≥ 8 (11,230 after merging the same signal seen from two ONs) was folded in all
six scans (`results/periodic_gj699/`, 496 s of forced folds with a coarse-channel cache;
the six searches took 334–509 s each with 2 workers on a 4-core machine):

| verdict | ON scans as ON | null: OFF scans play ON |
|---|---|---|
| IN-OFF (present in an OFF at ≥ 6) | 9,850 | 9,537 |
| PARTIAL (not in every ON) | 1,296 | 1,428 |
| **ON-ONLY (`101010`)** | **84** | **22** (`010101`) |

The 84 ON-only events: 70 in the Iridium allocation, 5 at 1373.875-1373.880 MHz (one
channel cluster, several period trials) in the air-surveillance radar band, among IN-OFF
neighbours at 1373.53-1373.93 MHz, 2 Inmarsat, 1 GNSS L1, and **6 at 1626.70–1626.72 MHz**, within
four channels of each other, 0.2 MHz above the listed Iridium edge (1626.5 MHz) — one
emitter at the edge of the Iridium band, not in our band list (post hoc, so this is E1). None is in
quiet spectrum. With the OFF scans playing ON the same filter passes 22 OFF-only events,
all in the Iridium, Inmarsat and AMT bands: bursty RFI happens to line up with one role
or the other at a rate of tens per cadence. About half (40 of 84) show snr 50–240 in
one ON and ~8–20 in another; an emitter that is bursty for 30 minutes gives exactly that.

**Verdict: no periodic pulse train with 10 ≤ P ≤ 70 s from Barnard's Star in quiet
spectrum** (outside the listed bands and 1610–1700 MHz) above integrated SNR ~12 per scan
for single-scan detection (90 %, duty 0.1, measured in three bands of one scan); the
completeness of the three-ON `101010` cadence code is not measured (≈ 75 % at 12 if scans
are independent). Inside the RFI bands the
search is blind (impulsive channels). The 84 ON-only events are all in or at the edge of
satellite and radar allocations and are not listed as candidates; the cadence null (22)
shows the ON/OFF test alone cannot separate them. This is E1: completeness was measured
in 48 of 256 coarse channels of one scan, and no flux/EIRP conversion was made.

## Not done

- Completeness outside 3 × 16 coarse channels of scan 0003; duty cycles other than 0.1;
  drifting pulse trains (mid-res channels are 2.86 kHz: a −0.12 Hz/s drift moves ~35 Hz
  in one scan, negligible; ~10 Hz/s would leave the channel).
- Cadence-level injection (a train planted in the three ONs and recovered as `101010`)
  only in the synthetic selftest, not in the real files.
- No conversion of SNR to flux or EIRP for pulsed signals (the radiometer equation for a
  duty-d train needs the mid-res SEFD per channel; not checked).
- riptide's own clustering/harmonic code (BLIPSS uses it) was not used; only its
  periodogram, mapped onto our period grid.
- `--no-common` (no common-mode division) and `--backend riptide` on the full cadence.
- A per-candidate time-shuffle (or kurtosis) test that would make the search work inside
  RFI bands (KNOWN_ISSUES #13); a second session of Barnard's Star; other targets.

    python periodicity.py selftest
    python periodicity.py reach
    python periodicity.py cadence data/gj699_midres --out-dir results/periodic_gj699 --workers 2
    python periodicity.py inject FILE.0002.h5 --coarse 150-165 --snrs 6,7,8,9,10,12,15,20,30,100 --per-level 40 --seed 1 --json r.json
    python periodicity.py search FILE.0002.h5 --coarse 150-165 --shuffle --out-dir out/   # time-shuffled null
