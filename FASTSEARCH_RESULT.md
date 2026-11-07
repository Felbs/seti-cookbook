# fastsearch: a one-pass drift search, measured against turboSETI — 2026-11-06

> **Evidence level E2 for the sensitivity comparison** (plants in real data, one coarse
> channel, both searches on the same planted arrays, six seeds). The full-file comparison
> (b) is E1-grade bookkeeping. The cadence search (c) was **not done**. See [EVIDENCE.md](EVIDENCE.md).
> All numbers below are from the current code (row-calibrated statistic); the first
> version's numbers, which were biased for fast drifts, are no longer quoted.

## What it is

`fastsearch.py` reads each HDF5 file in one pass, one coarse channel per hyperslab
(`data[:, 0, c0:c1]`, never the whole file, MemAvailable-guarded). For each coarse channel it:

1. flattens the bandpass with a 1024-channel **block median across frequency** of
   the time-mean spectrum, then divides each spectrum by it and by its own gain, and
   removes the small positive offset a median gain leaves on skewed spectra.
   A per-channel time median is never used: it would erase steady carriers.
2. blanks the DC spike (channel ±1 at the coarse centre);
3. de-Dopplers by brute force over **every drift k × 0.010510 Hz/s, k = −381…381,
   including exactly 0** (±4.0043 Hz/s searched when ±4 is asked);
4. sums each track over a **boxcar of w = round(|drift|·tsamp/df) channels**, the
   width a drifting tone is smeared over inside one spectrum. Rows with w ≥ 16 run
   on data pre-summed over s = 2^j ≤ w/8 channels;
5. gives two strengths per hit:
   - `snr` = **detection significance**. Each boxcar width's rows have their own noise
     median subtracted and are divided by their own noise spread, both measured on one
     real drift row of that width. On noise it has mean 0 and spread 1 at every drift,
     so one threshold means one false-alarm rate. This is the `.dat` SNR column. It
     equals the true SNR only at drift 0.
   - `snr_power` = **true radiometer SNR estimate**: the same boxcar power over the
     zero-drift channel-to-channel noise (`inject.radiometer_noise`'s definition, matched
     to 0.3 % on flat noise). It reads planted zero-drift tones back at 0.97–0.99 of
     their true SNR in real data. For a fast drifter it is 2–5× `snr` (selftest: a
     3.5 Hz/s tone of true SNR 180 has `snr` 34–36, `snr_power` 171–178).
6. keeps one best-per-channel map per boxcar-width octave, then suppresses
   neighbours **by track**, not by a fixed ±406-channel window (`select_hits`
   docstring);
7. keeps at most `--max-hits` hits per coarse channel (default 3000). A capped channel
   is **warned about** on stderr, flagged in the CSV (`capped` column) and listed in
   `<file>.summary.json`. (The `.dat` cannot carry the flag: turboSETI's `read_dat`
   parses every line after the 9th as a hit.)
8. writes a turboSETI-format `.dat` (turboSETI's own `read_dat` parses it; a review run
   of `find_event_pipeline` on six synthetic files gave exactly the planted ON-only event)
   and a CSV.

Backends: numpy (CI), numba if installed (same results to 1e-3, checked in the selftest),
cupy only with `SETI_COOKBOOK_GPU=1` (written, **never run**: no GPU here).

**Coarse channels in these files are 1,033,216 channels, not 2^20** (KNOWN_ISSUES #8).
fastsearch takes the width the way blimpy does. HDF5 chunks are 2^20 channels, so most
coarse channels share a chunk with the next. Each worker therefore takes contiguous
coarse channels with a chunk cache that holds those shared chunks. Measured on a synthetic
file with the same chunk geometry: 1.00× the file's bytes read with the cache, 1.76×
with h5py's default 1 MB cache. A chunk is still read twice at the boundary between two
workers' runs.

### Is the statistic calibrated? (selftest control + real data)

Mean and spread of `snr` on rows **not** used for the calibration, for every
(scrunch × width) in the plan:

| data | worst \|mean\| | spread |
|---|---|---|
| synthetic chi-square noise + bandpass + gain (selftest, 131k ch) | 0.008 | within 2.1 % of 1 |
| real GJ699 0003, coarse 155 (749,568 ch) | 0.012 | 1.000–1.018 |
| same real data, **first version** (review measurement) | 0.71 (3.78 Hz/s) | 1.02–1.08 |

The selftest fails if any width has |mean| ≥ 0.05 or spread off by ≥ 3 %. Checked:
it fails with the old noise estimate (0.062 / 6.1 %) and with both old pieces (0.93).
On synthetic noise the size of the real segment (749,568 channels × 763 drifts, two
seeds), there are 100 hits at `snr` ≥ 5, 6 at ≥ 5.5 and **0 at ≥ 6**. The 3 hits at
threshold 7 on the unplanted real channel are therefore real lines, not noise. Their
count (5/3/3/3 at thresholds 6/7/8/10) does not by itself fix a false-alarm rate.

## (a) Sensitivity: injection-recovery on real data

`fastsearch_recovery.py` uses the GJ699 ON scan 0003 (2016-03-10), coarse channel 155:
the 749,568 channels of the 1418.6–1421.6 MHz slice that lie inside that one coarse
channel. The DC spike is repaired identically for both searches. Doppler-smeared setigen
tones were planted (`inject.plant`), 3,366 per comparison over six seeds (11 drifts, ±,
16–18 SNR levels). The same array was searched by turboSETI 2.3.2 (max_drift 4,
min_drift 0) and by fastsearch at threshold 7. Both were scored one-to-one along the
track (`inject.match`, `drift_ok`). Null hits on the unplanted channel: fastsearch 3;
turboSETI 4 at its SNR 10 (the GJ699 cadence setting) and 7 at SNR 5. turboSETI at SNR 5
therefore runs at a **looser** false-alarm level than fastsearch.

90 %-recovery TRUE SNR per scan, pooled over ± drift and six seeds, log-interpolated.
Brackets: 95 % bootstrap range (plants resampled within each cell). The plant layout
avoids each search's null hits, so the two turboSETI runs used different layouts.

| \|drift\| (Hz/s) | fastsearch | turboSETI SNR 10 | gain | fastsearch | turboSETI SNR 5 | gain |
|---|---|---|---|---|---|---|
| 0 | 9 [7–9] | 91 [85–94] | 10.7× | 9 [8–10] | 114 [42–124] | 12.4× |
| 0.1 | 12 [11–13] | 65 [57–81] | 5.5× | 11 [10–11] | 53 [39–73] | 4.9× |
| 0.5 | 19 [18–21] | 182 [137–318] | 9.5× | 19 [19–20] | 208 [91–250] | 10.7× |
| 1 | 26 [24–28] | 510 [294–833] | 19.9× | 25 [24–26] | 269 [136–519] | 11.0× |
| 2 | 38 [36–39] | 636 [561–2000] | 16.8× | 39 [37–41] | 833 [371–1000] | 21.6× |
| 3.5 | 62 [53–67] | 941 [777–1063] | 15.2× | 59 [50–62] | 1882 [703–2075] | 32.0× |

So fastsearch reaches 90 % at **5–20× lower true SNR than turboSETI at its default SNR
10, and 5–32× lower when turboSETI runs looser at SNR 5**. The gain is smallest (~5×) for
slow drifts and largest for fast ones. fastsearch's own points are tight (±10 %).
turboSETI's are poorly determined: its curve flattens near 90 %, so its bootstrap ranges
span up to 3×, and the gain at a given drift is uncertain by the same amount. The first
version claimed 40 at 3.5 Hz/s. Part of that came from its biased statistic; the honest
figure is 59–62. Raw tables: `results/fastsearch_recovery_gj699_0003_ts{10,5}.json`
(`python fastsearch_recovery.py summary results/...json`).

Where the gain comes from (synthetic, `fastsearch.py` modes): true radiometer noise
instead of a bandpass-inflated one (~4×, KNOWN_ISSUES #1); drift 0 searched; and the
boxcar, which gives ~3.3× in `snr` at 3.5 Hz/s over one channel per spectrum (selftest:
34–36 vs 10.5–10.7). The scrunching of fast rows costs 0 % at ≤ 1 Hz/s and 6–15 % at
1.5–3.5 Hz/s against full resolution (`--widths full`, first-version measurement, not
repeated). In return it has fewer trials, so fewer false alarms at the same threshold.

Per scan, at the radiometer conversion of GJ699_CADENCE_RESULT.md (~0.7 Jy per unit of
true SNR here), 90 % recovery at true SNR ~12 for slow drifts would be ≈ 8 Jy. **This is
NOT a new limit on Barnard's Star.** That needs the full cadence searched with fastsearch
and every event explained, which (c) did not do.

## (b) Full ON file 0003 on the PC (current code)

Threshold 7, 3 workers, Nice 19, CPUQuota 300 %, PC otherwise lightly loaded.

- **Wall time: 11.1 min for all 256 coarse channels** (0.54 worker-minutes of reading).
  The first version took 21.4 min on a busier PC with the default chunk cache
  (1.5 worker-minutes of reading). turboSETI's log for the same file: 65.6 min, but it
  ran as one process alongside three other turboSETI runs.
- **Like-for-like CPU** (same 8 coarse channels 100–107, one process each, CPUQuota 100 %,
  PC load ~2–3 from other users): **fastsearch 51.8 s (numba) / 74.4 s (pure numpy),
  turboSETI 111.1 s.** That is 2.1× / 1.5× less CPU. Peak memory 715 / 620 MB vs
  1,778 MB. On the Mac (review measurement, one 749,568-channel array): 5.4 s numba,
  7.9 s numpy, turboSETI 13.1 s.
- **Hits: fastsearch 63,553, turboSETI 2,644.** fastsearch is far more sensitive, and it
  splits broadband interference into many tracks. **15 coarse channels hit the 3,000-hit
  cap** (86, 94, 102, 108–119: 1,622 MHz and the GPS L1 / Inmarsat 1,535–1,576 MHz
  band). All 15 were warned about and flagged.
- **turboSETI hits also found by fastsearch: 2,443 of 2,644 (92.4 %)**, 2,442 with
  matching drift (`fastsearch_compare.py`: track match within 3 channels + smear at
  both ends). The first version found 2,432 (92.0 %): the suppression change in between
  did **not** recover most of its misses. Each miss was re-probed by evaluating
  fastsearch's statistic along turboSETI's own track:
  - **123 in capped coarse channels** (cc 86, 102, 108, 109, 112–116).
  - **75 above-threshold-unmatched**: fastsearch's statistic along turboSETI's track is
    above threshold, but no kept fastsearch hit matches it. They sit at 1694.4–1694.7
    and 1691.1 MHz (cc 62–63, 24), 1618.2–1626.1 MHz (cc 85, 87, 88; 46), 1646.7 MHz (1),
    1250.0 MHz (2) and 1183.34 MHz (2). Most are drift tracks through broadband
    interference, where fastsearch keeps different (often wider, stronger) tracks than
    turboSETI. The two at 1183.34 MHz were found 4–5 channels away (fastsearch 639 and
    6,033 vs turboSETI 166 and 1,476, review check). The tool does not name the kept
    neighbour.
  - **3 below threshold**: all on the Galactic hydrogen emission at 1420.59–1420.61 MHz
    (cc 155). The time-mean power there is 2.5× the band level. turboSETI responds to
    the broad natural line's shape, which fastsearch's flattening removes.
- **With the cap raised** (`--max-hits 50000`, the 15 capped channels re-run, 3.5 min,
  122,323 hits in those 15 channels): **2,506 of 2,644 found (94.8 %)**, 2,505 with
  matching drift. Remaining misses: 134 above-threshold-unmatched (the 75 above plus 59
  in the former capped channels) and 4 below threshold (3 on HI, 1 at 1556.9 MHz). So
  in broadband-interference channels, fastsearch's hit list is long and picks different
  tracks than turboSETI's. A per-coarse-channel broadband flag would serve better than
  tens of thousands of hits (not built).
  Tables: `results/fastsearch_compare_gj699_0003_cap{3000,50000}.json`.

## (c) Six-file cadence + find_event_pipeline — NOT DONE

All six files were started with an earlier version, then stopped early to protect another
job on the compute machine. No event list exists from fastsearch yet. The 9
GEO-satellite events were not re-tested.

## Not verified

- GPU (cupy) path: never run.
- The cadence (c) and `find_event_pipeline` on real fastsearch `.dat` files (only on six
  synthetic files, in review).
- Completeness outside coarse channel 155; plants at |drift| > 3.5 Hz/s; a second
  observing session; full-file runs of the other five scans with the current code.
- The scrunch cost (6–15 % at 1.5–3.5 Hz/s) was measured with the first version only.
- The chunk-cache gain on the real 12.8 GB files was not isolated (reading fell from
  1.5 to 0.54 worker-minutes, but the PC load differed between the runs).

    python fastsearch.py selftest
    python fastsearch_recovery.py selftest
    python fastsearch_compare.py selftest
    python fastsearch.py FILE.h5 --max-drift 4 --snr 7 --workers 3 --out-dir out/ [--max-hits 3000]
