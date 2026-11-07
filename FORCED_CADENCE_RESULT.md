# Forced cadence check: Barnard's Star, every ON-only hit measured in all six scans (2026-11-06)

> **Evidence level E2** for the forced RFI/absence calls on slow drifts (≤ 0.5 Hz/s):
> checked with signals planted into the real data. **E1** at 1–3 Hz/s, where the
> measured sensitivity falls off (table below). Rung: nothing passes, so no rung-2 rows.
> See [EVIDENCE.md](EVIDENCE.md).

**What the tool does.** `cadence_forced.py` does not ask "did turboSETI log a hit
in the OFF scan?". It *measures* every scan. From the hit's frequency and drift it
predicts the signal's channel track in each scan, reads only that window, adds up the
power along the track (the box is as wide as the tone's Doppler smear), and compares
the best offset in a small search window with the same track shape at about 2,000
other offsets nearby. Every scan gets a number: a true radiometer SNR
(`inject.radiometer_noise` units) and the bar it had to reach. "Absent" means
"measured, and below the bar". A scan that could not be measured is never counted as absent.

**Data.** The same 2016-03-10 Green Bank cadence as [GJ699_CADENCE_RESULT.md](GJ699_CADENCE_RESULT.md)
(ON OFF ON OFF ON OFF, 2.84 Hz × 18 s, 1.1–1.9 GHz). Input: turboSETI's
find_event tables (filter 2: 1,029 rows; filter 3: 9 events × 3 rows). Run on the
home PC at lowest priority; 1,029 rows in 826 s.

## Results

| input | RFI | PASSES | ONE-OFF | PARTIAL | INCOMPLETE | NOT-RECOVERED | UNTESTED |
|---|---|---|---|---|---|---|---|
| filter 3 (27 rows = 9 events) | **27** | 0 | 0 | 0 | 0 | 0 | 0 |
| filter 2 (1,029 rows) | **767** | **0** | 115 | 7 | 0 | 140 | 0 |
| filter 2, OFF free-drift look off | 746 | 0 | 134 | 9 | 0 | 140 | 0 |
| filter 2, Earth curvature modelled | 767 | 0 | 116 | 6 | 0 | 140 | 0 |

- **The 9 filter-3 "survivors" (GEO satellites)** are present in **all three OFF
  scans** at true SNR 29–673 (z 20–600), on the same track. turboSETI listed no OFF hit
  there because each is part of a comb of carriers ~300 Hz apart (Inmarsat at 1538.888
  MHz sits beside 1538.8894 MHz, and the OFF hit list picked the neighbour). Before this,
  they were explained only by zero drift and band. Now the cadence test itself rejects them.
- **1674.9977 MHz (ON2 + ON3, −0.30 Hz/s): RFI from both hits.** On the predicted
  sky-like track it is present in OFF 0006 (true SNR 45, z 24). The free-drift look
  finds it in every OFF (z 34–39). It does not keep one drift. Measured scan by
  scan it is −0.33, −0.53, −0.20 and +0.21 Hz/s (`probe1675.py`). That is why OFF
  scans also get a free-drift look: anything at ±0.75 Hz/s of the track within ±100 Hz.
  This look is the only change that is not symmetric between ON and OFF, and it can
  only make a call stricter (more RFI). Without it, the ON3 hit's RFI call hangs on one
  OFF at z 6.4 (bar 6).
- **The other out-of-band rows from the July/November review:** 1426.63 MHz is in
  OFF 0004 at true SNR 52 → RFI. 1769.99 MHz is in OFFs 0004/0006/0008 at SNR 35/48 →
  RFI. 1646.655 MHz is ONE-OFF: true SNR 154 in ON1, below the bar (~23–30) in every
  other scan, which fits a one-time transmission. The 1420.59–1420.61 MHz rows (Galactic
  HI region) are NOT-RECOVERED: the tool sees no steady narrowband track even in their
  own scan.
- **ONE-OFF / PARTIAL** (122 rows): all but one are inside listed interference bands,
  mostly Iridium 1616–1626.5 MHz bursts, plus AMT and GOES.
- **NOT-RECOVERED** (140 rows): 108 of them drift faster than 1 Hz/s. Their own-scan
  z is ≤ 6.0 (median 3.9), so they are not a steady narrowband track this check can
  follow. This is a sensitivity limit of the check, not a claim about the signal. It is
  never counted as PASSES.
- **Earth-term curvature** (+3.8e-6 Hz/s² at 1.5 GHz) bends a track by at most 2.5
  channels over the cadence. The search window absorbs it by default. Modelling it
  moves 1 of 1,029 rows (PARTIAL → ONE-OFF). It does not matter here.

**No row PASSES.** There is nothing "not yet explained" left from this cadence.

## Controls (every run)

**Negative: random positions** (same coarse channels and drift mix, ≥ 2 kHz from any
turboSETI hit). These show how often "present" fires on empty sky:

| where | in own scan | other ON | OFF (with free-drift look) |
|---|---|---|---|
| quiet HI coarse channel (npz, 200 positions) | 0 % | 0 % | 0 % |
| filter-2 coarse channels (300) | 4.7 % | 8.8 % | 16.7 % (8.7 % without) |
| filter-3 coarse channels: GPS L1, Inmarsat (100) | 11 % | 30 % | 52 % |

In crowded bands, many "random" positions sit on real carriers that turboSETI did not
list. There, absence from an OFF is hard to prove, and a real signal in those bands
would often be called RFI. The tool reports this. It does not hide it.

**Positive: smeared tones planted with `inject.plant` into the six real GJ699 slices**
(1418.6–1421.6 MHz). A plant in the three ONs must come back PASSES; a plant in all six
scans must come back RFI. The hit is given with a 0.4-quantum drift error.
4 plants per cell, 36 + 96 in total.

| drift (Hz/s) | true SNR 15 | 20 | 30 | 40 | 50 | 80 |
|---|---|---|---|---|---|---|
| −0.12 | 8/8 | 8/8 | 8/8 | 6/6 | 8/8 | 6/6 |
| −0.5 | – | – | – | 6/6 | – | 6/6 |
| ±1.0 | 0/8 | 5/8 | 8/8 | 6/6 | 8/8 | 6/6 |
| +3.0 | 0/8 | 0/8 | 1/8 | – | 6/8 | – |

All 36 misses came out NOT-RECOVERED, ONE-OFF or PARTIAL (31, 2 and 3). No plant was
called the wrong way round: an ON-only plant was never called RFI, and an every-scan
plant was never called PASSES. The true SNR read back matches the planted value to
about 10 %.

## Not verified

- HTTP range reads were checked on one scan and one window: the bytes match the local
  slice exactly, and the read took 30 s. A whole cadence over HTTP was not run (download limit).
- The drift-error model (`drift_error()`) was fitted to synthetic plants. It is
  generous on purpose, but no test proves it covers every real signal.
- One cadence, one night. The "absent" bar is z ≥ 6, about true SNR 8–30 depending on
  drift and band (per-scan `bar_*` columns).

Files: `results_forced/forced_GJ699_f2.csv`, `forced_GJ699_f3.csv` (per scan: z, true
SNR, bar, present, reason), `forced_GJ699_summary.json` (counts, controls, options).

    python cadence_forced.py selftest
    python cadence_forced.py run --events events_f2.csv --files h5s.lst --dats dat/*.dat --out forced_f2.csv
    python cadence_forced.py control --files data/hi_*GJ699*.npz
    python cadence_verdict.py events_f3.csv --forced forced_f3.csv
