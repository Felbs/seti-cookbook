# Signal recovery rate (eta_P) of our turboSETI search: UCLA protocol, GJ699 data, 2026-11-06 (fixed 11-07)

> **Evidence level E2** (signals planted into real data and searched with the
> exact search settings), for ONE scan (GJ699 ON 0003). Not E3: the run has not
> been repeated on a second machine, and the plants were not blind. See [EVIDENCE.md](EVIDENCE.md).
> Plants were added **after channelisation** (setigen, Doppler-smeared), not
> into voltages as UCLA did.
> **Slicing (2026-11-07):** every number here used 2^20-channel slices that are not aligned
> with the 1,033,216-channel coarse channels (KNOWN_ISSUES #8). `etap.py` now cuts true
> coarse channels by default; `--nfine 1048576` (h5) or `--whole-slice` (npz) reproduces
> these results. The aligned re-run was not done.

**Headline.** At the UCLA strength (turboSETI should report SNR ~20), our
R9 turboSETI search (max_drift 4, min_drift 0, snr 10) recovers
**eta_P = 3.5 % (95 % CI 3.0-4.1 %; 173 of 4,906 plants)** of signals whose drift is
uniform over +-8.86 Hz/s. It finds every plant slower than 0.15 Hz/s, ~39 % at
0.15-0.5 Hz/s, and almost none faster. A simple smear-aware search on the same
data at the same strength recovers 78-90 %. So the low number comes from
turboSETI's settings (it sums one channel per spectrum, and searches only up to
4 Hz/s), not from the data.

## Protocol (`etap.py`)

| item | value |
|---|---|
| start frequency | uniform over each slice; the 1200-1340 MHz GBT notch is excluded |
| drift | uniform over +-8.86 Hz/s; subsets +-4 and +-2 Hz/s reported too |
| strength | expected SNR 20 **in turboSETI's units** (UCLA set S/N ~20 "upon recovery" in their own pipeline, which is not the same unit): true SNR = 20 x (turboSETI noise / radiometer noise) = 20 x 3.35-4.62 (median 3.99), i.e. **true radiometer SNR 67-92, median 80** |
| density | 75 plants per MHz (223 per 2.973 MHz slice); tracks kept 820+ channels apart, wider than turboSETI's own one-hit-per-window radius (2 x 406 channels) |
| match (UCLA) | one-to-one; frequency within +-6 Hz **and** drift within +-0.05 Hz/s. The plant's frequency is taken at the moment turboSETI reports its hit (the middle of the first spectrum, t = tsamp/2; see below) |
| match (ours) | `inject.match` (one-to-one along the whole track), and that plus `drift_ok` |
| search | turboSETI 2.3.2, max_drift 4, min_drift 0, snr 10, n_coarse_chan 1 (the R9 settings), on a 2^20-channel (2.973 MHz) HDF5 slice. **The slice is not an aligned PFB coarse channel**: these files have 256 coarse channels of 1,033,216 fine channels, so every slice contains one coarse-channel edge roll-off (down to 0.56 of the median level in the HI slice). R9 ran on the whole file, so its turboSETI noise estimate per coarse channel may differ from ours per slice; not measured |
| intervals | exact (Clopper-Pearson) 95 % binomial |

## Results: turboSETI backend, one scan (GJ699 ON 0003, 2016-03-10)

| data | plants | eta_P +-8.86 Hz/s | +-4 Hz/s | +-2 Hz/s |
|---|---|---|---|---|
| Mac: slice 1418.6-1421.6 MHz (HI), 7 seeds, re-run after the fix | 1,561 | 3.9 % [3.0, 5.0] | 8.9 % [6.8, 11.2] | 17.3 % [13.5, 21.6] |
| PC: 15 slices, 1621-1876 MHz (re-scored, see below) | 3,345 | 3.3 % [2.8, 4.0] | 7.1 % [5.9, 8.5] | 14.0 % [11.7, 16.6] |
| **all together** | **4,906** | **3.5 % [3.0, 4.1]** | **7.6 % [6.6, 8.8]** | **15.0 % [13.0, 17.2]** |
| same, outside `recipe_api.RFI_BANDS` | 4,460 | 3.4 % [2.9, 4.0] | 7.4 % [6.3, 8.6] | 14.6 % [12.5, 16.9] |

All three ways of matching agree within 0.1 points (UCLA 3.53 %, path 3.55 %,
path + correct drift 3.55 %). Of the 4,906 plants, 174 were found, 4,729 were
too faint for the search, 2 were found on the wrong path and 1 was hidden by a
stronger neighbour.

Recovery against |drift|, all 4,906 plants (UCLA match):

| \|drift\| Hz/s | 0-0.15 | 0.15-0.5 | 0.5-1 | 1-2 | 2-3.78 | 3.78-4 | 4-8.86 |
|---|---|---|---|---|---|---|---|
| plants | 85 | 220 | 304 | 545 | 979 | 134 | 2,639 |
| found | **100 %** | **39 %** | 0.7 % | 0 | 0 | 0 | 0 |

The largest drift turboSETI reported for any recovered or stray hit was
2.31 Hz/s. The per-drift curve gives the headline number directly:
(0.15 x 1.00 + 0.35 x 0.39) / 8.86 = 3.2 %.

**The all-band number is measured mostly outside RFI.** Only 9 % of the plants
(446) landed inside listed RFI bands (the sampled slices include Iridium
1610-1626 and GOES 1675-1700 MHz), but RFI_BANDS are 55 % of the searched band.
The RFI-band plants were recovered at 4.5 % (20/446), so weighting the two
strata by their share of the band gives 4.0 %, inside the all-band interval.
That weighting rests on Iridium and GOES only: the dense 1525-1610 MHz bands
and the radar/GNSS bands below the notch were never sampled, and there eta_P
could be much lower. `etap.py summary` prints both numbers, and `reach.py
--etap` prints this caveat under the all-band bound.

### Time reference of turboSETI's frequency (bug found in review, fixed 2026-11-07)

turboSETI reports a hit's frequency at the **middle of the first spectrum**
(t = tsamp/2 = 9 s), while a plant's start frequency is defined at t = 0. Over
the 174 recovered plants the offset fits 0.495 x drift x tsamp (through the
origin). The first version compared them at different times, so a plant found
exactly at 8.86 Hz/s would have sat 80 Hz "off" and failed the +-6 Hz test.
Every backend now carries `t_ref_frac` (turboSETI 0.5, reference 0, a callable
one is measured), and the UCLA test compares like with like. After the fix the
largest |offset| among recovered turboSETI plants is 5.1 Hz (4.3 Hz among
UCLA matches). The headline moved by one plant (172 -> 173 of 4,906; the
0.706 Hz/s plant at 10.6 Hz now counts). The -1.389 Hz/s plant still fails
UCLA, correctly: its offset is 5.1 Hz but its drift is 0.054 Hz/s off.
The Mac channel was re-run with the fixed code; the PC records were re-scored
from their stored path-matched hits (`etap.py rescore`). Re-scoring the old Mac
records the same way gives exactly the re-run's 61 matches, with every record
equal except 0.01 Hz rounding in 11 offsets.

## The same plants, searched by the reference backend (numpy brute-force, smear-aware)

`etap.py`'s own small search (`--backend reference`, +-8.86 Hz/s, box sum over
each spectrum's smear) on the HI slice, 3 seeds:

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

## Against UCLA's BL-like control (eta_P 5.7 % for +-8.88 Hz/s)

Source: Margot et al. 2023, AJ 166, 206 (arXiv:2308.02712), section 4.2 and
Table 2: a "BL-like process" recovered **570 of 10,000** injections (5.7 %,
exact 95 % CI 5.3-6.2 %). Their setup: chirps injected into 8-bit raw voltages,
start frequencies uniform over 1.15-1.73 GHz minus the 1.20-1.34 GHz notch,
drift uniform over +-8.86 Hz/s, strength set to S/N ~20 "upon recovery"; then
2^20-point FFTs (2.98 Hz), 51 spectra summed incoherently (17.11 s), and
turboSETI `-M 8.881784197 -s 10`; credited within +-6 Hz and +-0.05 Hz/s.
99.1 % of what it recovered drifted within +-1 Hz/s. We quote this; we did
not recompute it. (The "51" is the number of FFTs summed into each 17 s
spectrum, the same as in BL's own HSR files and in our data, so it is not a
difference between the two setups.)

**Close, and caused by the same thing, but not the same number.** Our 3.5 %
[3.0, 4.1] is about 2 percentage points below 5.7 %, and the two intervals do
not overlap. In both cases nearly everything faster than ~0.5 Hz/s is lost,
and 0.5 / 8.86 = 5.6 % on its own. We have NOT verified why ours is lower.
Possible reasons:
- **How "SNR 20" was set.** UCLA set S/N ~20 as measured on recovery (by their
  pipeline, as we read it); we set it from turboSETI's noise estimate on each
  slice. That estimate includes the bandpass, so true SNR is 3.4-4.6
  times larger. A setup that set "SNR 20" a little brighter would recover more
  of the 0.15-0.5 Hz/s plants, and those are the only drifts where recovery is
  partial.
- **Voltage versus channel injection.** Planting into voltages adds losses
  between fine channels (scalloping, leakage). That should make UCLA's number
  *lower* than ours, so this does not explain the gap.
- **Different search depth and band.** R9 ran turboSETI with max_drift 4, not
  8.88. That should not matter much (UCLA's recoveries were 99.1 % within
  +-1 Hz/s), but we did not test max_drift 8.88. Their band runs 1.15-1.73 GHz,
  ours was sampled at 1419-1422 and 1621-1876 MHz. Their resolution is
  2.98 Hz x 17.11 s, ours 2.84 Hz x 17.99 s, so smearing starts at a slightly
  lower drift for us. The number of spectra per scan in their files is not
  stated in the passage we read.
- The 5.7 % figure is quoted from the paper (checked against its text on
  2026-11-07) and was **not** recomputed here.

## Transmitter-fraction bound (`reach.py --etap`)

f_TX < [1 - (1 - CL)^(1/n)] / (eta_P x D), CL 95 %, D = fraction of the
1126-1876 MHz file band searched. The 1200-1340 MHz notch is removed
(D = 0.813), and RFI_BANDS are also removed for the outside-RFI figure
(D = 0.365). Drift adds no extra factor, because eta_P was measured over the
whole +-8.86 Hz/s range.

```
python reach.py --etap ETAP_SUMMARY_GJ699.json --n-stars 1
  all band     eta_P 0.0353  D 0.813  ->  bound 33.1 >= 1: NO constraint
               note: measured on a sample only 9% inside RFI_BANDS, which are 55% of the band
  outside RFI  eta_P 0.0343  D 0.365  ->  bound 75.8 >= 1: NO constraint
  for scale, same eta_P and D over 1000 stars: f_TX < 0.104 (all band), 0.239 (outside RFI)
```

One star, one scan, eta_P 3.5 %: **this null puts no limit on the fraction of
stars with transmitters.** The bound only means something across thousands of
stars, and only with a search that keeps fast drifts.

## What was NOT done or verified

- **Band coverage.** Only 16 of the file's 252 2^20-channel slices were sampled
  (the file holds 256 true PFB coarse channels of 1,033,216 fine channels;
  our slices are not aligned with them): 1419-1422 MHz, plus 15 slices between
  1621 and 1876 MHz. `etap.py h5 --nfine 1033216` cuts aligned coarse channels
  for future runs; that was not run. A planned 35-channel run across the band
  was stopped after 15 to free a shared machine.
- No repeat run on a second machine, so this is not E3. The planned PC re-run of
  the Mac's slice with the same seed did not happen.
- One scan only. eta_P here is for a single hit. Needing the signal in all three
  ON scans (find_event) would lower it further, and that was not measured.
- No voltage injection, so scalloping and leakage losses are not included.
- turboSETI's time reference (t = tsamp/2) is fitted on recovered plants up to
  1.39 Hz/s only; turboSETI found nothing faster to check it on.
- fastsearch is not plugged in yet. `CallableBackend(fn)` is the slot for it.
  When not given, its SNR scale is measured with tones drifting by 1 and 2
  steps of its drift grid (never drift 0, which turboSETI-style searches skip:
  zero-drift tones measured 8.3 instead of 4.26 on turboSETI), and its time
  reference is measured with bright 0.5-1 Hz/s tones. Checked on real data
  through the generic path with turboSETI as the callable: scale 4.44 (vs its
  analytic 4.26, +4 %); time reference raw fit 0.32, rounded to the 0.5 that the
  recovered plants fit. The rough raw fit means a new backend's time
  reference should be checked by hand.
- The setigen planting on the PC came from a private copy (`pip --target`,
  setigen 2.7.0), because the PC's venv lacks setigen. The repo venv there was
  left untouched.

## Reproduce

```
python etap.py selftest
python etap.py run data/hi_..._GJ699_0003.gpuspec.0000.npz --seeds 1,2,3,4,5,6 --json data/etap_mac.json   # ~2.5 min/seed
python etap.py h5 FULL_FILE.h5 --coarse 3,9,15,... --seeds 1 --json data/etap_pc.json        # per-channel seeds
python etap.py run FILE.npz --backend reference --max-drift 8.86 --per-mhz 50 --target-snr 85 --seeds 2,3
python etap.py rescore OLD_PRE_FIX.json --json RESCORED.json      # pre-11-07 runs only
python etap.py summary data/etap_mac_on0003_turboseti_tref.json data/pc_on0003_turboseti_35cc_rescored.json --json ETAP_SUMMARY_GJ699.json
python reach.py --etap ETAP_SUMMARY_GJ699.json --n-stars 1      # committed summary
```
