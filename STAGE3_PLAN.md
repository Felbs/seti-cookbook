# Stage 3 plan — raise the odds

Stage 1 stopped us fooling ourselves (exact archive selection, smeared plants,
evidence levels). Stage 2 made the search faster and honest about what it
misses (fastsearch, eta_P, targets, cross-star RFI veto, forced cadence check).
Stage 3 looks where a single scan cannot: **the same star on another date,
many epochs added together, and signals that repeat in time.**

Every item below is written down **before** it is run on real data. A pass/fail
gate that is changed after seeing the data drops the result one evidence level.

## 0. Before any real-data search (done on a laptop)

| step | tool | gate (must pass first) |
|---|---|---|
| 0.1 | multi-epoch census | lists every target with ≥ 2 fine-res sessions on different dates **in the same receiver band**, with real ON/OFF cadences (`bl_archive.py` rules). GJ699 checked explicitly. |
| 0.2 | one-command cadence run (fastsearch → find_event → forced check → verdict) | `--dry-run` prints the plan, time and disk; selftest on synthetic h5 |
| 0.3 | broadband-RFI flag in fastsearch | a synthetic broadband shelf is flagged and a narrow carrier inside it is still kept |
| 0.4 | coarse-channel-aligned slicing | recovery re-measured inside one true coarse channel; the GJ699 slow-drift limit (~56 Jy / 22 GW) is re-checked and corrected if it moves |

## 1. Cross-epoch matching in the star frame (`xepoch.py`)

A transmitter on a planet keeps its frequency **in the star's frame** from one
year to the next; interference is fixed to the telescope.

- Hits from each epoch → barycentric frequency at **mid-scan**, plus the star's
  own motion (proper motion, parallax, radial velocity and its secular change —
  for Barnard's Star about 21 Hz per year at 1.42 GHz).
- Matching by sorted merge under a stated tolerance (no fixed buckets), a
  telescope-frame OFF veto and a star-frame OFF veto, a drift-consistency test
  after removing Earth's predicted drift.
- False-alarm rate from per-epoch hit density over the searched bandwidth.
- **Gate:** the same star-frame tone planted into 2–3 real epochs is recovered;
  a tone fixed at the telescope is vetoed; selftest offline.
- **Evidence:** E2 if the gate passes on real epochs; E1 otherwise.

## 2. Stacking epochs in the star frame (`hi_starframe --stack` / stack)

Adding N epochs can push below one epoch's limit — **if** the signal really
lines up. We measure the gain, never assume √N.

- Per epoch: remove Earth's predicted drift, shift-and-add over a small grid of
  intrinsic drifts, move to the star frame, mask known interference, average.
- Epochs that do not share frequency coverage are refused up front.
- **Gate:** plants at low SNR in every epoch give a measured gain; an
  OFF-as-ON control gives the false-alarm rate. Window, epochs and threshold
  are pre-registered in the result file before the real stack is looked at.

## 3. Repeating signals — the pulsar bounty (`periodicity.py`)

Narrowband drift searches ignore pulses. A beacon could blink.

- Fold or FFA each channel group over a **pre-registered** period range that
  the chosen data product can actually reach (5-min fine-res scans have only
  16 spectra of 18 s; shorter periods need the mid- or high-time-resolution
  products).
- Count trials (channels × periods × duty cycles), label harmonics, give each
  candidate an ON/OFF detection code (e.g. `101010`) with stricter ON than OFF
  cutoffs.
- **Gate:** periodic pulse trains planted into real data at measured noise are
  recovered (fraction vs SNR); unplanted data and OFF-as-ON give the null rate.
- riptide-ffa optional; a pure-numpy path always runs.

## 4. Heavy runs (need a bigger machine; queued until it is free)

In this order, each with the gates above already passed:

1. Full GJ699 cadence through fastsearch → find_event → forced check →
   verdict (77 GB on disk already; ~1–2 h).
2. Cross-epoch and stacking runs on the best census targets (one coarse channel
   per epoch by HTTP range read where possible).
3. Three-class injection benchmark on real background (true / single-true /
   false) and the "where does a published event drop out?" probe on
   Enriquez 2017 / Price 2020 events.
4. Periodicity search on a real mid/high-time-resolution product.

## Rules carried over

True radiometer SNR, smeared plants, one-to-one matching, an evidence level
on every result, missing data never counts as a pass, nothing is ever called a
detection below rung 5.
