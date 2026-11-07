# Forced cadence check: Barnard's Star, every ON-only hit measured in all six scans (2026-11-06)

> **Evidence level.** **E2** for forced calls in *quiet* spectrum at slow drifts
> (≤ 0.5 Hz/s): there, signals planted into the real data were never called the
> wrong way. **E1** in *crowded* bands (where almost all of this cadence's hits sit):
> there an ON-only, sky-like plant came back **RFI 6 times in 60** and **OFF-OCCUPIED
> 32 times in 60** (table below), so an RFI call in a crowded band does not by itself
> show that the hit is local. **E1** at 1–3 Hz/s, where sensitivity falls off.
> Rung: nothing passes, so no rung-2 rows. See [EVIDENCE.md](EVIDENCE.md).

**What the tool does.** `cadence_forced.py` does not ask "did turboSETI log a hit
in the OFF scan?". It *measures* every scan. From the hit's frequency and drift it
predicts the signal's channel track in each scan, reads only that window, adds up the
power along the track (the box is as wide as the tone's Doppler smear), and compares
the best offset in a small search window with the same track shape at about 2,000
other offsets nearby. Every scan gets a number: a true radiometer SNR
(`inject.radiometer_noise` units) and the bar it had to reach. "Absent" means
"measured, and below the bar". A scan that could not be measured is never counted as absent.

**Two kinds of "something in the OFF".**
- **RFI**: present *on the predicted track* in an OFF, and a free-drift search in that
  same small window picks a drift that agrees with the track (within 2 × the drift
  error). This is evidence of the same emitter.
- **OFF-OCCUPIED**: an OFF has power near the track, but only through the OFF
  free-drift look (any drift within ±0.75 Hz/s, within ±100 Hz), or the on-track box
  fires while the best drift there is a different one (a carrier crossing the box).
  The cadence is failed, but whether it is the same emitter or a neighbour is
  **ambiguous**.

**Data.** The same 2016-03-10 Green Bank cadence as [GJ699_CADENCE_RESULT.md](GJ699_CADENCE_RESULT.md)
(ON OFF ON OFF ON OFF, 2.84 Hz × 18 s, 1.1–1.9 GHz). Input: turboSETI's
find_event tables (filter 2: 1,029 rows; filter 3: 9 events × 3 rows). Run on the
home PC at lowest priority; 1,029 rows in 1,056 s.

## Results

| input | RFI | OFF-OCCUPIED | PASSES | ONE-OFF | PARTIAL | INCOMPLETE | NOT-RECOVERED | UNTESTED |
|---|---|---|---|---|---|---|---|---|
| filter 3 (27 rows = 9 events) | **27** | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| filter 2 (1,029 rows) | **682** | 85 | **0** | 115 | 7 | 0 | 140 | 0 |
| filter 2, OFF free-drift look off | 682 | 64 | 0 | 134 | 9 | 0 | 140 | 0 |
| filter 2, Earth curvature modelled | 687 | 80 | 0 | 116 | 6 | 0 | 140 | 0 |

- **The 9 filter-3 "survivors" (GEO satellites)** are RFI on the same track: all 81
  OFF measurements (27 rows × 3 OFFs) are on-track at a consistent drift, true SNR
  29–673. turboSETI listed no OFF hit there because each is part of a comb of carriers
  ~300 Hz apart (Inmarsat at 1538.888 MHz sits beside 1538.8894 MHz, and the OFF hit
  list picked the neighbour). They are also explained by zero drift and band.
- **1674.9977 MHz (ON2 + ON3, −0.30 Hz/s): OFF-OCCUPIED from both hits.** It does not
  keep one drift: measured scan by scan it is −0.33, −0.53, −0.20 and +0.21 Hz/s
  (`probe1675.py`). In OFF 0006 the on-track box holds true SNR 45, but the best drift
  there is not the track's; the free-drift look finds power in every OFF (z 34–39, at
  drifts +0.2 to +0.3 Hz/s, within 100 Hz). So every OFF is occupied, but the tool no
  longer claims it has proved "same emitter".
- **The other out-of-band rows from the July/November review:** 1426.63 MHz is RFI
  (same track in OFF 0004, true SNR 52). 1769.99 MHz is OFF-OCCUPIED (on-track boxes in
  OFFs 0004/0006 at SNR 35/42, but at another drift). 1646.655 MHz is ONE-OFF: true SNR
  154 in ON1, below the bar (~23–30) in every other scan, which fits a one-time
  transmission. The 1420.59–1420.61 MHz rows (Galactic HI region) are NOT-RECOVERED.
- **OFF-OCCUPIED** (85 rows): 82 are inside listed interference bands. The other 3
  are 1769.99 MHz and the two 1674.9977 MHz hits. How the OFFs were occupied: crossing
  and free-drift look 51, free-drift look only 21, crossing only 13.
- **ONE-OFF / PARTIAL** (122 rows): all but one are inside listed interference bands,
  mostly Iridium 1616–1626.5 MHz bursts, plus AMT and GOES. The exception is 1646.655 MHz.
- **NOT-RECOVERED** (140 rows) = **the forced check could not test them.** 108 drift
  faster than 1 Hz/s. Their own-scan z is ≤ 6.0 (median 3.9): no steady narrowband
  track this check can follow. A brute-force look at two of them found a broadband
  burst (1555.504 MHz) and a carrier that switches on part-way through (1530.504 MHz).
  This is a limit of the check, not a finding about the signal, and it is **not** an
  explanation: `cadence_verdict.py --forced` prints it as information only.
- **Earth-term curvature** (+3.8e-6 Hz/s² at 1.5 GHz) bends a track by at most 2.5
  channels over the cadence. The search window absorbs it by default. Modelling it
  moves 5 rows between RFI and OFF-OCCUPIED and 1 row from PARTIAL to ONE-OFF. No row
  changes between passing and failing.

**What is left.** 0 rows PASS. 140 rows could not be tested by this check; every one
of them is explained by the band or zero-drift filters of `cadence_verdict.py`. With
`--forced`, `cadence_verdict.py` leaves **0 of 1,029** filter-2 events unexplained
(9 of 9 filter-3 events explained). Without `--forced` it leaves 5. Those 5 now rest on
forced verdicts: 1426.63 MHz (RFI, same track), 1646.655 MHz (ONE-OFF), and
1769.99 MHz plus the two 1674.9977 MHz hits (OFF-OCCUPIED, i.e. the cadence fails but
"same emitter" is not proved). No event rests on a NOT-RECOVERED or INCOMPLETE row.

## Controls (every run)

**Negative: random positions** (same coarse channels and drift mix, ≥ 2 kHz from any
turboSETI hit). These show how often "present" fires on empty sky:

| where | in own scan | other ON | OFF (with free-drift look) |
|---|---|---|---|
| quiet HI coarse channel (npz, 100 positions) | 0 % | 0 % | 0 % |
| filter-2 coarse channels (300) | 4.7 % | 9.5 % | 16.8 % |
| filter-3 coarse channels: GPS L1, Inmarsat (100) | 11 % | 30 % | 52 % |

**Positive: plants that should PASS or be RFI.** A plant in the three ONs must come
back PASSES; a plant in all six scans must come back RFI. The hit is given with a
0.4-quantum drift error. True SNR 40 unless stated.

| where the plants went | injector | ON-only plants | every-scan plants |
|---|---|---|---|
| quiet HI channel, real npz slices (drift −0.12/−0.5/+1, SNR 40/80) | setigen | 18/18 PASSES | 18/18 RFI |
| random positions over the whole 1.1–1.9 GHz h5 files (−0.12) | box | 15/15 PASSES | 14 RFI, 1 OFF-OCCUPIED |
| **3–20 kHz from a filter-2 hit** (crowded, −0.12/−0.5, 60 + 30) | box | **16 PASSES, 6 RFI, 32 OFF-OCCUPIED, 6 NOT-RECOVERED** | 27 RFI, 2 OFF-OCCUPIED, 1 NOT-RECOVERED |
| **3–20 kHz from a filter-3 hit** (GPS/Inmarsat, 30) | box | **6 PASSES, 2 RFI, 16 OFF-OCCUPIED, 6 NOT-RECOVERED** | – |

In crowded bands an OFF usually holds some unrelated carrier near the track. A real
sky signal there would usually fail the cadence (RFI or OFF-OCCUPIED 38 times in 60
near filter-2 hits). Before the RFI/OFF-OCCUPIED split, all 38 would have been called
RFI (an independent review measured 19 of 30). The split moves most of them to the
honest "ambiguous" class, but 6 in 60 (10 %) near filter-2 hits and 2 in 30 (7 %) near
filter-3 hits still meet the same-track test: a neighbour carrier at about the same
drift. That is why RFI calls in crowded bands are E1. The same-track test also costs
some real same-track RFI: 2 of 30 every-scan plants came out OFF-OCCUPIED (still a
cadence failure, never PASSES).

**Sensitivity in the quiet channel** (setigen plants into the six npz slices, 4 per
cell, ON-only and every-scan combined, 96 plants, re-run after the split):

| drift (Hz/s) | true SNR 15 | 20 | 30 | 50 |
|---|---|---|---|---|
| −0.12 | 8/8 | 8/8 | 8/8 | 8/8 |
| +1.0 | 2/8 | 6/8 | 8/8 | 8/8 |
| +3.0 | 0/8 | 1/8 | 2/8 | 8/8 |

All 29 misses came out NOT-RECOVERED (26) or PARTIAL (3). In the quiet channel no plant
was called the wrong way round: no ON-only plant was RFI or OFF-OCCUPIED, and no
every-scan plant was PASSES. The true SNR read back matches the planted value to about
10 %.

The PC's venv has no setigen, so the h5 plants used the setigen-free box injector
(`cadence_forced.py control --injector box`). The selftest checks it against
`inject.plant`: same power per spectrum; centroids differ by up to 0.2 channel because
setigen samples the 1-channel box at channel centres.

## Not verified

- HTTP range reads were checked on one scan and one window: the bytes match the local
  slice exactly, and the read took 30 s. A whole cadence over HTTP was not run (download limit).
- The drift-error model (`drift_error()`) and the same-track drift tolerance were fitted
  to synthetic plants. They are generous on purpose, but no test proves they cover every
  real signal.
- The crowded-band controls used only drifts −0.12 and −0.5 Hz/s at SNR 40.
- One cadence, one night. The "absent" bar is z ≥ 6, about true SNR 8–30 depending on
  drift and band (per-scan `bar_*` columns).

Files: `results_forced/forced_GJ699_f2.csv`, `forced_GJ699_f3.csv` (per scan: z, true
SNR, bar, present, `how_*` = track / crossing / nearby, the free-drift look's
`nearz_*`, `neard_*`, `nearoff_*`, and the reason), `forced_GJ699_summary.json`
(counts, controls, options).

    python cadence_forced.py selftest
    python cadence_forced.py run --events events_f2.csv --files h5s.lst --dats dat/*.dat --out forced_f2.csv
    python cadence_forced.py control --files data/hi_*GJ699*.npz
    python cadence_forced.py control --files h5s.lst --near-events events_f2.csv --where ON --injector box
    python cadence_verdict.py events_f3.csv --forced forced_f3.csv
