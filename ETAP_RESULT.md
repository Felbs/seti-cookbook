# Signal recovery rate (eta_P) of our turboSETI search: UCLA protocol, GJ699 data, 2026-11-06

> **Evidence level E2** (signals planted into real data and searched with the
> exact search settings), for ONE scan (GJ699 ON 0003). Not E3: the run has not
> been repeated on a second machine, and the plants were not blind. See [EVIDENCE.md](EVIDENCE.md).
> Plants were added **after channelisation** (setigen, Doppler-smeared), not
> into voltages as UCLA did.

**Headline.** At the UCLA strength (turboSETI should report SNR ~20), our
R9 turboSETI search (max_drift 4, min_drift 0, snr 10) recovers
**eta_P = 3.5 % (95 % CI 3.0-4.1 %, 4,906 plants)** of signals whose drift is
uniform over +-8.86 Hz/s. It finds every plant slower than 0.15 Hz/s, ~39 % at
0.15-0.5 Hz/s, and almost none faster. A simple smear-aware search on the same
data at the same strength recovers 78-90 %. So the low number comes from
turboSETI's settings (it sums one channel per spectrum, and searches only up to
4 Hz/s), not from the data.

## Protocol (`etap.py`)

| item | value |
|---|---|
| start frequency | uniform over each coarse channel; the 1200-1340 MHz GBT notch is excluded |
| drift | uniform over +-8.86 Hz/s; subsets +-4 and +-2 Hz/s reported too |
| strength | expected SNR 20 **in turboSETI's units**: true SNR = 20 x (turboSETI noise / radiometer noise) = 20 x 3.35-4.62 (median 3.99), i.e. **true radiometer SNR 67-92, median 80** |
| density | 75 plants per MHz (223 per 2.93 MHz coarse channel); tracks kept 820+ channels apart, wider than turboSETI's own one-hit-per-window radius (2 x 406 channels) |
| match (UCLA) | one-to-one; start frequency within +-6 Hz **and** drift within +-0.05 Hz/s |
| match (ours) | `inject.match` (one-to-one along the whole track), and that plus `drift_ok` |
| search | turboSETI 2.3.2, max_drift 4, min_drift 0, snr 10, n_coarse_chan 1, on a one-coarse-channel HDF5 slice (exact R9 settings) |
| intervals | exact (Clopper-Pearson) 95 % binomial |

## Results: turboSETI backend, one scan (GJ699 ON 0003, 2016-03-10)

| data | plants | eta_P +-8.86 Hz/s | +-4 Hz/s | +-2 Hz/s |
|---|---|---|---|---|
| Mac: coarse channel 1418.6-1421.6 MHz (HI), 7 seeds | 1,561 | 3.8 % [3.0, 4.9] | 8.7 % [6.7, 11.1] | 17.0 % [13.2, 21.3] |
| PC: 15 coarse channels, 1621-1876 MHz | 3,345 | 3.3 % [2.8, 4.0] | 7.1 % [5.9, 8.5] | 14.0 % [11.7, 16.6] |
| **all together** | **4,906** | **3.5 % [3.0, 4.1]** | **7.6 % [6.5, 8.8]** | **14.9 % [12.9, 17.1]** |
| same, outside `recipe_api.RFI_BANDS` | 4,460 | 3.4 % [2.9, 4.0] | 7.3 % [6.3, 8.6] | 14.5 % [12.4, 16.8] |

All three ways of matching agree within 0.2 points (UCLA 3.51 %, path 3.53 %,
path + correct drift 3.53 %). Of the 4,906 plants, 174 were found, 4,729 were
too faint for the search, 2 were found on the wrong path and 1 was hidden by a
stronger neighbour.

Recovery against |drift|, all 4,906 plants (UCLA match):

| \|drift\| Hz/s | 0-0.15 | 0.15-0.5 | 0.5-1 | 1-2 | 2-3.78 | 3.78-4 | 4-8.86 |
|---|---|---|---|---|---|---|---|
| plants | 85 | 220 | 304 | 545 | 979 | 134 | 2,639 |
| found | **100 %** | **39 %** | 0.3 % | 0 | 0 | 0 | 0 |

The largest drift turboSETI reported for any recovered or stray hit was
2.31 Hz/s. The per-drift curve gives the headline number directly:
(0.15 x 1.00 + 0.35 x 0.39) / 8.86 = 3.2 %.

Only 9 % of the plants landed inside listed RFI bands (the sampled channels
include Iridium 1610-1626 and GOES 1675-1700 MHz), and they were recovered about
as often as the rest. So "outside RFI" changes eta_P very little here. That
could change in the dense bands we did not sample (1525-1610 MHz).

## The same plants, searched by the reference backend (numpy brute-force, smear-aware)

`etap.py`'s own small search (`--backend reference`, +-8.86 Hz/s, box sum over
each spectrum's smear) on the HI coarse channel, 3 seeds:

| strength | plants | eta_P UCLA match | path + drift |
|---|---|---|---|
| true SNR 85 (same as turboSETI's plants) | 447 | **78 % [74, 82]** | **90 % [86, 92]** |
| true SNR 20 | 149 | 4.0 % [1.5, 8.6] | 4.0 % |

So turboSETI loses ~96 % of these signals because of how it searches. The data
do not hide them. The reference backend also shows a limit of the UCLA match
itself: above 6 Hz/s, a search that adds up a wide smear box pins the start
frequency and drift less exactly. 73 % of those plants were found on the right
track, but only 51 % fell inside +-6 Hz and +-0.05 Hz/s. A search that adds up
the smear will score lower on the UCLA match than on the track match.

## Against UCLA's BL-like control (eta_P ~5.7 % for +-8.86 Hz/s)

**Close, and caused by the same thing, but not the same number.** Our 3.5 %
[3.0, 4.1] is about 2 percentage points below 5.7 %, and our 95 % interval does
not include it. In both cases nearly everything faster than ~0.5 Hz/s is lost,
and 0.5 / 8.86 = 5.6 % on its own. We have NOT verified why ours is lower.
Possible reasons:
- **How "SNR 20" was set.** We set it from turboSETI's noise estimate on each
  coarse channel. That estimate includes the bandpass, so true SNR is 3.4-4.6
  times larger. A setup that set "SNR 20" a little brighter would recover more
  of the 0.15-0.5 Hz/s plants, and those are the only drifts where recovery is
  partial.
- **Voltage versus channel injection.** Planting into voltages adds losses
  between fine channels (scalloping, leakage). That should make UCLA's number
  *lower* than ours, so this does not explain the gap.
- **Different data and resolution.** Our scan uses 2.84 Hz x 18 s spectra,
  16 per scan. A different number of spectra moves the 0.15 Hz/s point where
  smearing begins.
- The 5.7 % figure is quoted from our literature notes and was **not**
  recomputed here.

## Transmitter-fraction bound (`reach.py --etap`)

f_TX < [1 - (1 - CL)^(1/n)] / (eta_P x D), CL 95 %, D = fraction of the
1126-1876 MHz file band searched. The 1200-1340 MHz notch is removed
(D = 0.813), and RFI_BANDS are also removed for the outside-RFI figure
(D = 0.365). Drift adds no extra factor, because eta_P was measured over the
whole +-8.86 Hz/s range.

```
python reach.py --etap ETAP_SUMMARY_GJ699.json --n-stars 1
  all band     eta_P 0.0351  D 0.813  ->  bound 33.3 >= 1: NO constraint
  outside RFI  eta_P 0.0341  D 0.365  ->  bound 76.3 >= 1: NO constraint
  for scale, same eta_P and D over 1000 stars: f_TX < 0.105 (all band), 0.24 (outside RFI)
```

One star, one scan, eta_P 3.5 %: **this null puts no limit on the fraction of
stars with transmitters.** The bound only means something across thousands of
stars, and only with a search that keeps fast drifts.

## What was NOT done or verified

- **Band coverage.** Only 16 of 252 coarse channels were sampled: 1419-1422 MHz,
  plus 15 channels between 1621 and 1876 MHz. The PC run planned 35 channels
  across the band, but **it was stopped after 15**. The PC's D-layer experiment
  had missed its 23:01, 23:11 and 23:31 UTC slots. Its last log line was
  22:41:30Z, with overflows=0, and its process was busy on the CPU without
  logging. The D-layer process kept stalling after our jobs stopped (23:16Z),
  while other users' compile jobs held the load at 13-18. We did not touch it,
  and someone needs to look at it.
- No repeat run on a second machine, so this is not E3. The planned PC re-run of
  the Mac's coarse channel with the same seed did not happen.
- One scan only. eta_P here is for a single hit. Needing the signal in all three
  ON scans (find_event) would lower it further, and that was not measured.
- No voltage injection, so scalloping and leakage losses are not included.
- The +-6 Hz start-frequency convention was checked only on recovered slow
  plants (|offset| <= 2.6 Hz). turboSETI found no fast plant to check it on.
- fastsearch is not plugged in yet. `CallableBackend(fn)` is the slot for it,
  and its SNR scale is measured automatically when not given.
- The setigen planting on the PC came from a private copy (`pip --target`,
  setigen 2.7.0), because the PC's venv lacks setigen. The repo venv there was
  left untouched.

## Reproduce

```
python etap.py selftest
python etap.py run data/hi_..._GJ699_0003.gpuspec.0000.npz --seeds 1,2,3,4,5,6 --json data/etap_mac.json   # ~2.5 min/seed
python etap.py h5 FULL_FILE.h5 --coarse 3,9,15,... --seeds 1 --json data/etap_pc.json        # per-channel seeds
python etap.py run FILE.npz --backend reference --max-drift 8.86 --per-mhz 50 --target-snr 85 --seeds 2,3
python etap.py summary data/etap_*.json --json data/etap_summary.json
python reach.py --etap ETAP_SUMMARY_GJ699.json --n-stars 1      # committed summary
```
