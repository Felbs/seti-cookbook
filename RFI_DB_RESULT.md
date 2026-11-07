# Cross-star RFI database (`rfi_db.py`) - 2026-11-07

> **Evidence level: tool validation, not a search.** The veto is measured against
> real hit tables (E1-style bookkeeping) and against smeared plants in real noise
> (false-veto rate, E2-style). It explains interference; it can never promote a signal.

**The idea.** A signal from a star appears only when the telescope points at that
star. If the same frequency *with the same drift* is in the hit list of a
DIFFERENT star - that night or another night - it is made on Earth. The database
holds those hit lists; `rfi_db.py query` asks it about one hit, `explain` about a
whole find_event table, and `cadence_verdict.py --rfi-db` adds one line per event.

## What is in it

| source | hits | how obtained | SNR units |
|---|---|---|---|
| Enriquez et al. 2017 (ApJ 849, 104), all hits of 692 stars, GBT L-band 2016-17, 133 nights, 2,457 target names | 28,858,709 | `http://blpd0.ssl.berkeley.edu/lband2017/All_hits_turbo_seti.csv` (10.8 GB, CC BY 4.0), streamed on the PC through `tools/rfidb_reduce_enriquez.py` (raw CSV never stored); gitignored `data/rfi_db_enriquez.npz` | turboSETI (SNR >= 25, \|drift\| <= 2 Hz/s) |
| Price et al. 2020 (AJ 159, 86) L-band events | 9,507 | `https://seti.berkeley.edu/listen2019/lband2019_events.csv` (909 kB, sha256 86ffab24...) | turboSETI |
| Price et al. 2020 Parkes 10-cm events | 77 | `https://seti.berkeley.edu/listen2019/pks10cm2019_events.csv` | turboSETI |
| our GJ699 cadence, 6 scans incl. 3 OFFs, 2016-03-10 | 15,237 | PC `~/seti-work/gj699_57457/dat/*.dat` | turboSETI |
| our GJ699 mid-res file | 87 | `data/star_GJ699.dat` - drift unresolved, flagged `drift_valid=0`, never vetoes | turboSETI |
| blc1 control, Parkes UWL 2019-04-25 and 04-29, 366 public files | 625 | `tools/rfidb_blc1_fetch.py`: HTTP Range reads of +-5 kHz windows from the public .fil files, own de-drift search | **true** (radiometer) |

Everything except Enriquez is in the checked-in **`rfi_seed.csv.gz` (309 kB, 25,533 hits)**.
Rebuild the full DB with
`python rfi_db.py build data/rfi_db_enriquez.npz --enriquez-dir DIR_OF_CHUNKS` (needs ~2 GB RAM; done on the PC).

Not used: (c) Green Bank's public RFI monitoring is published as JPG plots per
band (`greenbankobservatory.org/gbt-rfi-scans/`); the machine-readable archive
(`rfi.gb.nrao.edu`) needs a login, and NRAO's newer DRIFT tool did not offer a
download we could use. (d) `recipe_api.RFI_BANDS` is merged unchanged as an
allocation table (`rfi_db.allocations()`), with review notes on four entries
(AMT upper edge 1525 vs 1535 MHz; the 1610-1626.5 "Iridium uplinks" label - that
range is the MSS uplink band, Iridium is TDD 1616-1626.5; the 1215-1400 radar
entry vs the GBT notch; an overlap). Nothing deleted.

## The rule (full text in the `rfi_db.py` docstring)

Match = the frequency tracks of the two hits overlap (each `[f, f + drift x T_obs]`,
widened by 2 channels + f x 10 m/s / c, i.e. 33 Hz at 1 GHz, 47 Hz at 1.42 GHz; all
tables are topocentric, so Earth-made RFI keeps its frequency) **and** the drifts
agree within 0.02 Hz/s + 20 % + one drift quantum each. The query's own target
family (GJ699, GJ699_OFF, ...) is always left out: a repeat from the SAME star is
evidence for it. Count n_t = other target families; compare with the number
expected by chance from the density of drift-compatible hits in the surrounding
+-200 kHz; score = -log10 Poisson p; **veto when score >= 3**.

## Measured

**GJ699 filter-2 events (1,029, PC `events_f2.csv`), leave-the-target-out:**

| | explained by other targets |
|---|---|
| all other targets, any session | **419 / 1,029 = 40.7 %** |
| other targets, the GJ699 night (MJD 57457) excluded too | 418 / 1,029 = 40.6 % |
| ... only events inside Enriquez's search space (SNR >= 25, \|drift\| <= 2) | 59.0 % of 217 |
| ... events with SNR < 25 (fainter than anything Enriquez recorded) | 36.5 % of 787 |
| matched by >= 1 other target but held back by the chance test | 425 more (crowded satellite bands) |
| the 9 events outside every RFI_BANDS allocation | 0 / 9 (5 are on the HI line, see GJ699_CADENCE_RESULT.md) |

By band: GOES 81 %, AMT 64 %, Inmarsat 54 %, GNSS L1 40 %, Iridium 5 %. Zero-drift
events 66 %, drifting ones 7 % - Iridium's and other LEO satellites' Doppler
changes from scan to scan, so a frequency+drift key cannot follow them (by design;
`RFI_BANDS` still covers them). Filter-3 (the 9 geostationary events, 27 rows):
23 / 27 rows vetoed; 1538.888 MHz (Inmarsat) matched 21 other stars but in a band
so crowded that 11 were expected by chance, so it is not vetoed by the DB.

**False vetoes (unique signals killed).** 234 Doppler-smeared setigen plants
(drifts 0 to +-3.5 Hz/s, true SNR 100-1000) in the real ON 0003 coarse channel,
turboSETI (max_drift 4, snr 10), 164 recovered and matched one-to-one:
**0 / 164 vetoed in place** (1418.6-1421.6 MHz); moved to 20 random frequencies
each across 1.1-1.9 GHz (measured drift kept): **11 / 3,280 = 0.34 %** vetoed,
**0 / 1,251 outside RFI_BANDS**.

**How much band the veto removes** (3,000 random probes per drift, full DB):
drift 0: 0.90 %, -0.12 Hz/s: 0.33 %, +0.5: 0.07 %, -1: 0.23 %, +2: 0.13 %.
Frequency-only keying (drift ignored) would remove 1.37 % at +0.5 Hz/s, ~20x more;
`RFI_BANDS` remove 62.8 % of 1.1-1.9 GHz.

**blc1 control - end to end on the PUBLIC data.** The Parkes UWL files of
Smith et al. 2021 / Sheikh et al. 2021 are in the BL archive (`paperName=Sheikh21a`,
14 GB per 30-min scan). Range reads of a few kHz per spectrum (~1 MB per scan and
window) and our own de-drift search found blc1 in the same five ProxCen scans the
paper lists (Supp. Table 1): 982.002396 / .002472 / .002571 / .002644 / .002808 MHz,
drift +0.027 / +0.029 / +0.033 / +0.025 / +0.012 Hz/s, true SNR 11 / 30 / 18 / 15 / 11
(paper: 982.00241-.00281, 0.031-0.011 Hz/s, SNR 12.6 / 28.1 / 21.2 / 17.7 / 10.2).
- **Recurrence alone does NOT flag blc1**: one other target matches (HIP98959_R,
  2019-04-25, 982.002552 MHz, +0.015 Hz/s), score ~1.9 < 3. (Other 04-25 targets
  show 982.0022 MHz at near-zero drift and J1456-6843 shows 982.002132 at
  +0.031 Hz/s - 340 Hz away, outside our 40 Hz slack.)
- **Intermodulation relatives do**: in the SAME scans, hits at 982.0026 - 135 x
  2.000004 (~712.0019), - 63 x (~856.0021) and + 40 x (~1062.0026) MHz with similar
  drift, tracking blc1's own frequency wander. Chance expectation ~0.005;
  **blc1 flagged in 5 / 5 scans** (`python rfi_db.py blc1`).
- Plainly: this is a **known-answer** test. We fetched only the frequencies the
  paper names (blc1 and three lookalikes, all from Sheikh et al.), and the clock
  spacing 2.000004 MHz is theirs; our chance rate comes from hits in those
  windows. It shows the machinery works on real data, not that it would have
  found the clock blind.

## Not verified

- One target star's events (GJ699) for the leave-target-out number; one coarse
  channel for the plants. Enriquez hits are SNR >= 25 and |drift| <= 2 only, so
  faint or fast RFI is under-represented; Enriquez SNRs are from the 2017 turboSETI.
- Enriquez rows carry a nominal 2.79 Hz / 300 s (not read per file).
- Session = telescope + UTC MJD day; a night crossing 0 UT splits in two.
- The intermod rule is off by default and only knows the Parkes 2019 clock.
- Satellites with large Doppler (GNSS, Iridium) are not followed across sessions.
