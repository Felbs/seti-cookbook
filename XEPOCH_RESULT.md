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
L band, sessions MJD 57587, 57667 and 57804 (2016-07-18, 2016-10-06,
2017-02-20 [calendar dates corrected after registration; the MJDs are what was registered]), the coarse channel that holds the star-frame hydrogen line, all
three ON scans of each session; same window rule, threshold rule and drift
grid. GJ725 A and B share one archive name and one pointing; the star model
uses the first targets.json row (`Star.from_targets("GJ725")`).

---

## 1. Multi-epoch census (`tools/epoch_census.py`)

> Archive metadata only: no signal search, no E-level.

**Question:** which stars of the 1,709-star BL list (targets.json) have ≥ 2
valid fine-resolution (.0000) cadences on DIFFERENT dates in the SAME band?
An epoch is one session (project, MJD day, centre frequency) with a real
cadence: OFFs named `<name>_OFF` (`bl_archive.cadences`) or ABACAD with other
stars as OFFs (`targets.abacad`). Light metadata queries only (one file list
per star, one time-window query per session, ≤ 2 requests/s, raw replies
cached under `data/`), plus one HTTP header read per receiver configuration
for the frequency coverage. Output: `EPOCH_CENSUS.csv` (one row per star and
band, ranked by distance then total GB) and `EPOCH_CENSUS.json` (sessions,
MJDs, files, sizes, pointing offsets, coverage).

**Barnard's Star (GJ699): no second L-band epoch exists in the open archive.**
Its only fine-resolution cadences are L band 2016-03-10 (MJD 57457,
AGBT16A_999_81, 1501.5 MHz, OFFs named GJ699_OFF) and S band 2016-10-28
(MJD 57689, AGBT16B_999_55, ABACAD). Its HIP name (HIP87937) has no archive
rows. So every GJ699 cross-epoch or stacking result here uses PSEUDO-epochs
(the scans of one session), and says so.

## 2. Cross-epoch matching in the star frame (`xepoch.py`)

> **Evidence level: tool validation (E2-style injection gate)** — no real
> multi-date search of GJ699 is possible (one L-band date). Real-date results
> for GJ725 are in section 4.

**What it does.** Every hit of every scan (turboSETI or fastsearch `.dat`, or
fastsearch run directly on fetched coarse channels) is moved to the star's
rest frame at MID-scan time, `f_star = f_obs · γ_s(1 + rv(t)/c) / (γ_o(1 +
v_los(t)/c))`, with the star's position propagated by proper motion and
parallax and its RV including the secular term; the Earth drift predicted at
mid-scan is removed before drifts are compared. Hits that survive a
topocentric OFF veto (±500 Hz, static or carried by the hit's drift) and a
star-frame OFF veto (±30 Hz, OFF hits of every date) are clustered per epoch
(same star-frame frequency within ±10 Hz in ≥ `min_ons` ON scans) and matched
across epochs by a sorted merge under ±10 Hz — no buckets. An epoch is a
session (scans < 6 h apart), never a scan index. Groups must agree in
Earth-removed drift (±0.05 Hz/s); groups whose dates are too alike in Earth
velocity to separate a star from RFI are marked UNTESTABLE.

**Tolerance, stated.** ±10 Hz in the star frame = 3.6 channels of 2.79 Hz:
channel quantisation (±1.4 Hz), one drift step over half a scan (0.0102 Hz/s ×
137 s = ±1.4 Hz), the half-spectrum time convention (≤ 1 Hz), unmodelled
gravitational redshift (< 0.3 Hz between dates) and star-model errors (< 1 Hz
for Barnard's Star between dates years apart).

**Barnard's Star secular term, verified.** μ²d = 4.539 m/s/yr; astropy's
propagated RV changes by 4.541 m/s/yr; at 1420.4 MHz that is **21.5 Hz/yr**
(the review's "~21 Hz/yr" is right). With the J2000 direction and a constant
RV instead, a star-frame beacon lands 465 / 368 / 494 Hz off on 2016-03-10,
+200 d and +700 d: ~100 Hz scatter between dates, ten times the tolerance
(KNOWN_ISSUES #12). The selftest finds a three-date beacon with the full model
and loses it with the naive one.

**False alarms, stated correctly.** With n_i surviving clusters per epoch over
the common star-frame bandwidth B, p_i = n_i·2tol/B and the expected chance
groups spanning ≥ k epochs are (B/2tol)·Σ_{m≥k} P(exactly m)·g(m), where g(m)
is the anchor rule's acceptance factor (g(2) = 1, g(3) = 1.5 exactly). Monte
Carlo through the same matcher (selftest, 3,000/2,500/4,000 random clusters
in 2 MHz, 60 trials): ≥ 2 epochs 292.7 vs 290.5 predicted; ≥ 3 epochs 4.80 vs
4.48. A second, empirical estimate (shift test: each epoch's clusters shifted
by its own random offset, 40 rounds) keeps the real clustering of RFI combs.

**Gate 1 — selftest (offline, ~20 s, CPU only).** 11 checks pass: secular
term; v_los and Earth drift identical to `earth_term` (Voyager-gated); a
three-date beacon matched (full model) and lost (naive frame); RFI in ON and
OFF vetoed; two tones 6 Hz apart across a 10 Hz bucket edge matched; three
scans of one session = one epoch; the false-alarm formula vs Monte Carlo; and
a setigen end-to-end run on synthetic noise (three dates, 1 ON + 1 OFF each,
fastsearch): the star-frame tone is matched on all 3 dates, the topocentric
ON+OFF tone is vetoed in 3/3 ONs, the ON-only topocentric tone is never
matched across dates, and no other group passes.

**Gate 2 — injection on real GJ699 noise, three PSEUDO-epochs** (`python
xepoch.py inject`, ~4 min on 4 cores). The three ON/OFF pairs of the
2016-03-10 cadence, clean part of coarse channel 155 (749,568 channels,
1419.43-1421.56 MHz), each pair relabelled +0 / +200 / +400 days so the frame
physics acts over a year (Earth's velocity toward the star differs by up to
52 km/s between the labels). This is one real session with date LABELS - it
tests the chain on real noise and real RFI, not a second date. fastsearch
(threshold 10, ±0.5 Hz/s) on every scan; 20 star-frame tones (all three ONs,
Doppler-smeared setigen, observed drift = predicted Earth drift), 10
topocentric tones (all ONs and OFFs, −0.02 Hz/s), 10 ON-only topocentric
tones; min 1 ON per epoch (one ON per pseudo-epoch).

| true SNR per scan | matched on all 3 dates | on ≥ 2 dates | found per ON scan | topocentric groups | ON-only topocentric groups | other PASS |
|---|---|---|---|---|---|---|
| 10 | 0/20 | 1/20 | 2/60 | 0 | 0 | 0 |
| 12 | 1/20 | 3/20 | 13/60 | 0 | 0 | 0 |
| 15 | 10/20 | 17/20 | 47/60 | 0 | 0 | 0 |
| 20 | 18/20 | 19/20 | 57/60 | 0 | 0 | 0 |
| 30 | 20/20 | 20/20 | 60/60 | 0 | 0 | 0 |

The matching loses nothing beyond detection: at SNR 15, (47/60)³ × 20 = 9.6
expected on all three dates, 10 observed. Every topocentric ON+OFF tone found
was vetoed (30 of the 40 topocentric vetoes at SNR ≥ 20; the other 10 are the
data's own lines); every ON-only topocentric tone survived the vetoes, as it
must, and was never matched across dates. The unplanted data give 3-4 hits per
scan, all vetoed: 0 groups, 0 expected. Chance veto of a unique tone: 0.4-0.8 %
per ON hit. **90 % matched on all three dates at true SNR 20 per scan** (single-scan
fastsearch detection is the limit: 95 % found per scan at 20).

## 3. Star-frame epoch stack (`stack.py`)

> **Evidence level E2** for the pseudo-epoch run (one session, planted
> Doppler-smeared tones through the same pipeline, OFF-as-ON control). It is
> a single-session result: the three "epochs" are the three ON scans of
> 2016-03-10.

**Gate — selftest (offline, ~15 s):** a drifting track de-drifted into one
channel; 4 synthetic dates × 2 scans with a tone at true SNR 4 per scan:
measured gain 1.98 (sqrt 4 = 2); a noise-only stack has robust spread 1.000
and no z > 4.5; a topocentric tone with a sky-like drift does not stack (epoch
peak 6.6, stack peak 4.4: it lands kHz apart in the star frame, and even the
two scans of one date do not line up); an epoch from another band is refused.

**Run 1 — GJ699 pseudo-epochs, as pre-registered.** Grid 749,516 bins of
2.834 Hz (star frame), 0.0 % masked by the OFF rule (a quiet channel: 3-4
fastsearch hits per scan). 351 empty windows: 99th percentile 5.05 → threshold
**6.00** (the floor). In the pre-registered window (1420.405752 MHz ± 3 kHz,
star frame): stack max z **3.65**; single pseudo-epochs 3.75, 4.34, 4.06;
OFF-as-ON 4.66. **NULL.** OFF-as-ON control: 1 of 351 empty windows above 6.0
= **0.28 % false alarms per window** (the threshold was set at the ON stack's
99th percentile; the OFF stack agrees).

Measured from 40 plants per level (setigen, smeared, true radiometer SNR per
scan, star-frame-fixed, in all three ONs):

| true SNR per scan | 3-bin SNR, stack | 3-bin SNR, one pseudo-epoch | measured gain | detected, stack | detected, one epoch |
|---|---|---|---|---|---|
| 2 | 1.98 | 1.18 | 1.68 | 5 % | 1 % |
| 3 | 2.84 | 1.70 | 1.68 | 5 % | 2 % |
| 4 | 3.71 | 2.21 | 1.67 | 30 % | 4 % |
| 5 | 4.57 | 2.73 | 1.67 | 60 % | 7 % |
| 6 | 5.44 | 3.25 | 1.67 | **92 %** | 15 % |
| 8 | 7.16 | 4.29 | 1.67 | 98 % | 44 % |
| 10 | 8.89 | 5.33 | 1.67 | 98 % | 84 % |
| 12 | 10.62 | 6.37 | 1.67 | 100 % | **96 %** |

**Measured gain 1.67 against an ideal sqrt(3) = 1.73** (97 %): the loss is
the star-frame interpolation between scans whose grids are offset by Earth's
rotation (~37 Hz between scans). 90 % detection moves from true SNR 12 per
scan (one pseudo-epoch) to **6 per scan** (stack). With BL's SEFD 10 Jy, 2.84 Hz
channels and 288 s scans (`reach.py`, 1 Hz signal), that is

> no hydrogen-line beacon held at 1420.405752 MHz in Barnard's Star's frame,
> steady through the three 5-minute ON scans of 2016-03-10, brighter than
> ~4.2 Jy (~1.7 GW EIRP at 5.96 ly), at 90 % completeness.

Compare hi_starframe's single-scan E2 null (< ~12 Jy ≈ 5 GW, which divides by
a 0.68 planting efficiency that the smeared setigen plants here do not need;
the same per-scan SNR 12 here would be 8.4 Jy). Same night, same 15 minutes:
this is a deeper limit on the same data for a beacon that is on in all three
scans, not a new date.

**Changed during the work, stated openly** (none changes the window, the
epochs, the threshold rule or the detection rule):
1. The first null-window picker (random, ≥ 6 kHz apart) gave only 255 windows,
   short of the pre-registered ≥ 300. Replaced by non-overlapping tiling (351).
   The threshold stayed at the 6.0 floor either way (99th pct 4.98 → 5.05);
   the window statistic (3.65) did not change.
2. The first gain estimator (max over ±1 bin, minus its noise floor) gave
   2.0-2.1 at low SNR, i.e. MORE than sqrt(3): a max is not linear in the
   signal. Replaced by the 3-bin sum, whose SNR is mean at the plants over the
   spread at 400 empty places; it gives 1.67-1.68 at every level.
