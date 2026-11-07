# fastsearch: a single-read drift search, measured against turboSETI — 2026-11-06

> **Evidence level E2 for the sensitivity comparison** (plants in real data, one coarse
> channel, both searches on the same planted arrays, six seeds). The full-file comparison
> (b) is E1-grade bookkeeping. The cadence search (c) was **not done**. See [EVIDENCE.md](EVIDENCE.md).

## What it is

`fastsearch.py` reads each HDF5 file once, one coarse channel per hyperslab
(`data[:, 0, c0:c1]`, never the whole file, MemAvailable-guarded). For each coarse channel it:

1. flattens the bandpass with a 1024-channel **block median across frequency** of
   the time-mean spectrum, then divides each spectrum by it and by its own gain.
   A per-channel time median is never used: it would erase steady carriers.
2. blanks the DC spike (channel ±1 at the coarse centre);
3. de-Dopplers by brute force over **every drift k × 0.010510 Hz/s, k = −381…381,
   including exactly 0** (±4.0043 Hz/s searched when ±4 is asked);
4. sums each track over a **boxcar of w = round(|drift|·tsamp/df) channels**, the
   width a drifting tone is smeared over inside one spectrum. Rows with w ≥ 16 run
   on data pre-summed over s = 2^j ≤ w/8 channels;
5. quotes SNR in **true radiometer units**. The noise matches
   `inject.radiometer_noise` to 0.3 % on flat noise (selftest), and the power estimate
   `snr_power` reads planted zero-drift tones back at 0.96–0.99 of their true SNR in
   real data;
6. keeps one best-per-channel map per boxcar-width octave, then suppresses
   neighbours **by track**, not by a fixed ±406-channel window (`select_hits`
   docstring);
7. writes a turboSETI-format `.dat` (turboSETI's own `read_dat` parses it) and a CSV
   (true SNR, power SNR, drift, frequency at t = 0, coarse channel, local noise).

Backends: numpy (CI), numba if installed (same results to 1e-3, checked in the selftest),
cupy only with `SETI_COOKBOOK_GPU=1` (written, **never run**: no GPU here).

**Coarse channels in these files are 1,033,216 channels, not 2^20** (KNOWN_ISSUES #8).
fastsearch takes the width the way blimpy does.

## (a) Sensitivity: injection-recovery on real data, same false-alarm level

`fastsearch_recovery.py` uses the GJ699 ON scan 0003 (2016-03-10), coarse channel 155:
the 749,568 channels of the 1418.6–1421.6 MHz slice that lie inside that one coarse
channel. The DC spike is repaired identically for both searches. Doppler-smeared setigen
tones were planted (`inject.plant`), 3,366 over six seeds (11 drifts, ±, and 16–18 SNR
levels). The same array was searched by turboSETI 2.3.2 (max_drift 4, min_drift 0, snr 10)
and by fastsearch, and both were scored one-to-one along the track
(`inject.match`, `drift_ok`).

**Same false-alarm level:** on the unplanted channel turboSETI reports 4 hits at its
SNR 10. fastsearch reports 3 at true SNR 7, the threshold used for the full files below,
and 4 at 6.5.

90 %-recovery TRUE SNR per scan (pooled over ± drift and six seeds, log-interpolated):

| \|drift\| (Hz/s) | fastsearch (thr 7) | fastsearch (thr 6.5) | turboSETI (snr 10) | gain at thr 7 |
|---|---|---|---|---|
| 0 | 9 | 8 | 91 | 10.7× |
| 0.1 | 11 | 10 | 62 | 5.8× |
| 0.5 | 19 | 15 | 186 | 10.1× |
| 1 | 22 | 20 | 294 | 13.4× |
| 2 | 33 | 31 | 741 | 22.8× |
| 3.5 | 40 | 39 | 1406 | 35.5× |

turboSETI reproduces the published calibration (GJ699_CADENCE_RESULT.md: ~70 at
≤0.15, ~300 at 1, ~1200 at 3.5 Hz/s). Raw tables: `results/fastsearch_recovery_gj699_0003_thr{7,6.5}.json`
(`python fastsearch_recovery.py summary results/...json`). Every fastsearch miss was
below threshold. None was lost to a wrong drift or a stronger neighbour.

Where the gain comes from (synthetic, `fastsearch.py` modes): true radiometer noise
instead of a bandpass-inflated one (~4×, KNOWN_ISSUES #1); drift 0 searched; and the
boxcar, which gives 2.1× SNR at 1 Hz/s, 2.7× at 2 and 3.2× at 3.5 over one channel per
spectrum. The scrunching of fast rows costs 0 % at ≤ 1 Hz/s and 6–15 % at 1.5–3.5 Hz/s
against full resolution (`--widths full`), for ~3× less CPU. Partly in return it has
fewer trials, so fewer false alarms at the same threshold.

Per scan, at the radiometer conversion of GJ699_CADENCE_RESULT.md (~0.7 Jy per unit of
true SNR here), 90 % recovery at true SNR ~11 for slow drifts would be ≈ 8 Jy. **This is
NOT a new limit on Barnard's Star.** That needs the full cadence searched with fastsearch
and every event explained, which (c) did not complete.

## (b) Full ON file 0003 on the PC

Run with the first version of the code (commit 9678b59: brute force + boxcar + scrunch,
simpler suppression), threshold true SNR 7, 3 workers, Nice 19, CPUQuota 300 %.

- **Wall time: 21.4 min for all 256 coarse channels** (each byte read once, 1.5
  worker-minutes of reading). turboSETI's log for the same file: 65.6 min, but it ran
  as one process alongside three other turboSETI runs.
- **Like-for-like CPU** (same 8 coarse channels, one process each, same limits, busy PC):
  turboSETI 94 s, fastsearch 80 s. Per coarse channel, CPU is about the same. The wall-time
  gain comes from reading once and running parallel workers, not from cheaper arithmetic.
  The final code was not benchmarked like this. In the interrupted final runs it took a
  median 18.6 s per coarse channel per worker, against 13.3 s for the first version.
  Those runs had two jobs at once and a heavy compile from another agent on the PC, so
  the final code is **probably slower** than the first.
- **Hits: fastsearch 58,141, turboSETI 2,644.** fastsearch is far more sensitive. It also
  splits broadband interference into many tracks; 81 coarse channels' worth of
  turboSETI hits fell in channels where fastsearch hit its 3,000-hits cap.
- **turboSETI hits also found by fastsearch: 2,432 of 2,644 (92.0 %)**, of which 2,431
  agree in drift (`fastsearch_compare.py`, track match within 3 channels + smear at
  both ends). Each of the 212 misses was re-probed by evaluating fastsearch's statistic
  along turboSETI's own track:
  - **128 "suppressed"** (above threshold, not kept): narrow carriers inside broadband
    interference (Inmarsat 1545 MHz, cc 113; 1618–1626 MHz, cc 85–88; 1183 MHz, cc 236).
    Wide 24-channel boxcar tracks laid over the broadband signal out-scored them and
    removed them. **This was fixed after the run**: per-width-band candidate maps plus a
    measured exception in the suppression. Selftest control: a true-SNR-40 carrier on a
    3.5 Hz/s, SNR-3000 drifter's track is lost under the old rule and kept under the new
    one. **The full file was not re-run with the fix** (see below), so the improved hit
    fraction is not measured.
  - **81 "capped"**: coarse channels 102, 108, 109, 112 (GPS L1, Inmarsat), where
    fastsearch's 3,000-hit cap was reached.
  - **3 below threshold**: all on the Galactic hydrogen emission at 1420.59–1420.61 MHz
    (cc 155). The time-mean power there is 2.5× the band level. fastsearch's statistic
    along those tracks is 2.5–5.5. turboSETI's SNR 12–18 there comes from the broad
    natural line's shape, which fastsearch's flattening removes. This is the intended
    behaviour.

## (c) Six-file cadence + find_event_pipeline — NOT DONE

All six files were started with the final code (two 3-worker jobs). They were **stopped
at 188 of 256 coarse channels** because the D-layer radio test on the PC
(`dlayer-run.service`) stopped logging after 2026-11-06T22:41:30Z, while the service stayed
"active". It normally writes every 30 min. No `overflows=[1-9]` line appeared, and it
was still silent 10 min after these jobs stopped. Another agent's C++ build was loading
the PC at normal priority at the same time, so the cause is unknown. No event list
exists from fastsearch yet. The 9 GEO-satellite events were not re-tested.

## Not verified

- GPU (cupy) path: never run.
- The final code on a full file: only the first version ran end-to-end; the final
  version's sensitivity is measured only on the Mac channel (a).
- Completeness outside coarse channel 155; plants at |drift| > 3.5 Hz/s; a second
  observing session.
- find_event_pipeline on fastsearch `.dat` files: only `read_dat` parsing is tested.

    python fastsearch.py selftest
    python fastsearch_recovery.py selftest
    python fastsearch_compare.py selftest
    python fastsearch.py FILE.h5 --max-drift 4 --snr 7 --workers 3 --out-dir out/
