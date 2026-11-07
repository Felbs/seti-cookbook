# Evidence levels — how much each result is worth

Every result in this repo carries two labels.

**Rung** (from the definition ladder, for anything we *found*): 0 hit · 1 not
explained by any veto · 2 only when pointed at the target, in every ON ·
3 behaves like a distant source · 4 seen again on another date · 5 seen by
another instrument or team · 6 engineered information. Below rung 5 the honest
word is "not yet explained", never "aliens".

**Evidence level** (for how a search was *done*, which decides what a null or
a limit is worth):

| level | name | what it takes | what you may say |
|---|---|---|---|
| E0 | screen | one file, or ON scans only; no OFF comparison | "nothing obvious"; never a null |
| E1 | cadence | full ON/OFF cadence, every survivor explained; sensitivity from the radiometer equation only | "null", with an *estimated* limit |
| E2 | calibrated | E1 + signals planted into the real data and searched with the same settings; limit quoted where ≥ 90 % are recovered | "null down to X Jy / W" |
| E3 | reproduced | E2 + rerun independently: a second search code, a second machine, or blind plants the analyst did not know about | the strongest claim this repo makes |
| VOID | void | the method was shown unable to support the claim | nothing |

## Every result, labelled

| result | file | level | rung | note |
|---|---|---|---|---|
| Voyager 1 carrier + sidebands | `sideband_pairs.py`, `earth_term.py gate` | positive control | — | drift predicted to 0.2 % |
| Own-rig RF detector panel (GPS, ADS-B, AIS, FT8 …) | `agent.py selftest`, CYCLO/COMB/ENTROPY_RESULT.md | detector control | — | 7/7 flagged, controls silent |
| Epsilon Eridani mid-res | EPSILON_ERIDANI_RESULT.md | VOID | — | drift unresolved on .0002 files (drift_audit.py) |
| 10-star sweep | SWEEP_RESULTS.md | VOID / E0 | — | same reason; single ON files |
| Barnard's Star 1400–1500 MHz (July) | cadence_GJ699_*.json | E1, incomplete | killed | strict re-run after the fix not recorded |
| Barnard's Star hydrogen-line beacon | HI_STARFRAME_RESULT.md | E2 | — | null < ~12 Jy ≈ 5 GW; own local-noise detector; plants not smeared (small effect at −0.11 Hz/s) |
| Barnard's Star full cadence 1.1–1.9 GHz | GJ699_CADENCE_RESULT.md | **E2** (slow drifts), E1 above 0.5 Hz/s | — | null < ~56 Jy ≈ 22 GW for \|drift\| ≤ 0.15 Hz/s; completeness measured in one coarse channel, noise ratio band-wide |
| Barnard's Star keyed carriers | `keying.py` | E1 | — | selftest on real noise; no completeness curve yet |
| SDSS laser lines, 199 stars | LASER_LINES_RESULT.md | E2 | 1 (3 survivors) | 2 of 3 in telluric bands (post hoc); follow-up needs new spectra |
| Cross-star RFI database (veto) | RFI_DB_RESULT.md, `rfi_db.py` | tool validation | — | full DB (PC, 1.8 GB RAM): 40.7 % of 1,029 GJ699 filter-2 rows explained by other stars (seed-only default: 0.9 %); false veto 0/164 in place, 0.34 % band-wide; blc1 relatives present (known-answer) but not significant by its own chance test, 0/5 flagged |

## Rules

1. A result file states its level in its first lines, and README quotes it.
2. A limit from the radiometer equation alone is E1, however careful.
3. Plants must be smeared like real signals (law L10) and scored one-to-one
   along the track (`inject.py`), or the level drops back to E1.
4. Strength is quoted in TRUE radiometer SNR. turboSETI's own SNR is
   typically ~4× smaller on raw Green Bank coarse channels (law L13).
5. Moving a result up a level needs the new evidence committed beside it;
   moving it down needs only the reason.
