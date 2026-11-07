# Known issues

Status: **Open** · **Mitigated** (worked around here) · **Won't fix** (upstream, documented) · **Closed**.

### 1. turboSETI's SNR is ~4× smaller than radiometer SNR on raw Green Bank data — Mitigated
- **Symptom:** planted tones come back at ~0.23× their radiometer SNR (`recovery.py`).
- **Cause:** turboSETI 2.3.2 `helper_functions.comp_stats` divides by the 5–95 %-trimmed std of the whole time-summed coarse channel, which includes the polyphase bandpass shape (~57 % peak-to-peak). Measured ratio 3.2–6.1× (typ. 4.1×) over 62 coarse channels, 1.15–1.86 GHz, GJ699 2016-03-10.
- **Impact:** a radiometer limit at "SNR 10" is ~4× (with drift-tree losses ~7–8×) too optimistic.
- **Workaround:** quote limits from injection-recovery in TRUE SNR (`inject.radiometer_noise`), never from turboSETI's threshold.

### 2. turboSETI never searches drift exactly zero — Won't fix
- **Cause:** `find_doppler.py` tests `abs(drift_rate) > min_drift`, strictly; `min_drift` may not be negative.
- **Impact:** a perfectly steady tone is found only at ±1 drift bin, at ~half strength. Harmless for sky sources (Earth's term is never zero) and it suppresses Earth-locked RFI, but explains the +0.0105 Hz/s "zero-drift" satellites.

### 3. turboSETI searches less drift than asked — Won't fix
- **Symptom:** `max_drift=4` searched ±3.78 Hz/s on GJ699 fine-res files (largest drift in the hit list).
- **Workaround:** quote the searched range from the output (law L11).

### 4. turboSETI reads each coarse channel many times — Mitigated
- **Symptom:** ~18× re-reads of each file; the 77 GB GJ699 cadence took hours.
- **Workaround:** `fastsearch.py` reads each file in one pass, one coarse channel per hyperslab, with a chunk cache so chunks shared by neighbouring coarse channels are not re-read (synthetic file with GBT chunk geometry: 1.00× the file's bytes vs 1.76× with h5py's default cache). It searches drift 0 and the whole asked range, and quotes true radiometer SNR (FASTSEARCH_RESULT.md). Measured like-for-like on 8 coarse channels of GJ699 0003: 51.8 s CPU (numba) or 74.4 s (numpy) vs turboSETI 111.1 s; a full 12.8 GB file took 11.1 min with 3 workers.

### 5. Archive names match by prefix — Mitigated
- **Symptom:** `?target=HIP2` returns 1,982 files of other stars.
- **Workaround:** `bl_archive.py` (exact names, real OFFs, one session and band).

### 6. Some cadences' OFF scans carry other stars' names — Open
- **Symptom:** Epsilon Eridani (GJ144) has only ON scans under its name; its OFFs are filed under other targets, which the archive cannot be searched for without a name.
- **Impact:** `bl_archive.py` reports no valid cadence rather than building a wrong one.
- **Workaround:** `targets.py` `abacad()` rebuilds these cadences from a time-window archive query (`target=%`), with the same session/band/product/gap/separation rules; not yet merged into `bl_archive.py`.

### 7. Mid-resolution (.0002) drift vetoes are void — Mitigated
- **Cause:** drift quantum ~9.8 Hz/s, wider than any planet's drift (`drift_audit.py`). `star_sweep.py` prints UNRESOLVED.

### 8. 2016 GBT fine-res coarse channels are 1,033,216 channels, not 2^20 — Mitigated
- **Symptom:** the GJ699 2016-03-10 files have 264,503,296 channels = 256 coarse channels of 1,033,216 (2.9296875 MHz at 2.8355 Hz); the DC spike sits at k × 1,033,216 + 516,608 (checked in the 0003 slice: channel 232,960 of the slice = coarse 155's centre).
- **Impact:** the `hi_starframe` .npz slices (1,048,576 channels from ch0 = 160,432,128 = 153 × 2^20, i.e. starting 27 % into coarse 155) straddle two coarse channels; `recovery.py` ran turboSETI on the whole slice, with its DC blanking at the slice centre instead of at the real spike and a coarse boundary inside. Re-measured on the clean 749,568-channel part of coarse 155 with the spike repaired (`fastsearch_recovery.py`), turboSETI's 90 % points are close to the published ones (62 / 294 / 1406 at 0.1 / 1 / 3.5 Hz/s in the first six-seed run vs ~70 / ~300 / ~1200; a later six-seed run with a different plant layout gave 65 / 510 / 941, so turboSETI's 90 % points are uncertain by up to ~2×, FASTSEARCH_RESULT.md).
- **Workaround:** `fastsearch.py` takes the coarse width the way blimpy does (2^20 if it divides, else 2.9296875 MHz for GBT); `fastsearch_recovery.py` uses only the part of a slice inside one coarse channel.

### 9. find_event's "absent from OFF" means "no hit logged", not "measured absent" — Mitigated
- **Symptom:** all 9 GJ699 filter-3 events (GEO satellites) are present in all three OFF scans at true SNR 29–673; the 1674.9977 MHz track is in both later OFFs.
- **Cause:** an OFF counts as clean when turboSETI logged no hit near it; carrier combs (one hit per neighbourhood) and RFI that changes drift between scans slip through.
- **Workaround:** `cadence_forced.py` measures every scan along the predicted track (plus a free-drift look in OFFs); `cadence_verdict.py --forced`.
- **Limit:** in crowded bands an OFF often holds an unrelated carrier near the track. A sky-like ON-only plant 3–20 kHz from a filter-2 hit came back RFI 6/60 and OFF-OCCUPIED 32/60, so a forced RFI call there is E1, not proof the hit is local. NOT-RECOVERED/INCOMPLETE rows were not tested and never count as explained.

## Step 3: cross-epoch work (epoch_census, xepoch, stack)

Numbers 10-14 below may be renumbered when merged with other Step 3 branches.

### 10. Some archive "ON" scans point far from the star they are named after — Open
- **Symptom:** `tools/epoch_census.py` compares every ON scan's archive pointing (`ra`/`decl`, which agree with the file header's `src_raj`/`src_dej` where checked) with the star's position propagated to the scan date. GJ65A's two L-band sessions (MJD 57573 and 57579) point at Dec −16.730° while UV Ceti is at −17.951°: 73' away, ~9 beams at 1.5 GHz. Same pattern in the same projects (AGBT16A_999_232/233/236/237) for GJ54.1 (79'), HIP103039 (76'), HIP101027 (57'), HIP1242 (24'), HIP104139 (12.7') and GJ866 (5.8'); some single S-band sessions too (HIP42333 56', HIP14954 17.1', HIP12444 13.7', HIP11783 11.2').
- **Cause:** unknown; the archive metadata and the header agree, so either the telescope pointed there or both copied the same wrong coordinates.
- **Impact:** a search of those scans is not a search of that star; a multi-epoch census that ignores it ranks GJ65A first.
- **Workaround:** the census flags `pointing_ok = false` when the offset exceeds half the GBT beam (12.45'/f_GHz) and counts `n_days_pointing_ok`; `xepoch.py fetch` skips such sessions.

### 11. SIMBAD's rvz_radvel is not always a stellar radial velocity — Mitigated
- **Symptom:** targets.json has GJ866 at RV = 6,824.7 km/s.
- **Impact:** a star-frame frequency off by ~32 MHz at 1.42 GHz; a secular term computed with it is meaningless.
- **Workaround:** the census flags `rv_ok = false` beyond ±600 km/s; `xepoch.Star.from_targets` refuses such a star unless an RV is given by hand.

### 12. The star frame between dates needs the star's motion, not its J2000 entry — Mitigated
- **Symptom:** for Barnard's Star a beacon fixed in the star frame, converted with the J2000 direction and a constant RV (what a pipeline using one SIMBAD position does), lands 465, 368 and 494 Hz from where it should on 2016-03-10, +200 d and +700 d - ~100 Hz scatter between dates, 10× a 10 Hz matching tolerance.
- **Cause:** the direction moves 10.4"/yr (Earth's velocity projected on it changes by ±11 m/s by 2016) and the RV changes by the secular (perspective) acceleration mu^2 d = +4.54 m/s/yr, 21.5 Hz/yr at 1420 MHz (checked: astropy's propagated RV equals mu^2 d).
- **Workaround:** `xepoch.frame` propagates position and RV to every scan's mid-time (astropy `apply_space_motion`); `xepoch.py selftest` shows the naive frame loses a three-date beacon that the full model finds. The absolute star-frame frequency still carries the RV's own error (±0.5 km/s = ±2.4 kHz) and the RV's unknown reference epoch (up to ~80 m/s for Barnard's Star); both are constant between dates and do not affect matching.

### 13. A topocentric OFF veto of ±500 Hz kills unique tones by chance in crowded bands — Mitigated
- **Symptom:** a tone that IS only in the ON is vetoed whenever any OFF hit of the same session lies within ±500 Hz of it (where it sits, or where its drift carries it).
- **Impact:** measured per ON hit on real data: 0.7 % over the whole 1.1-1.9 GHz GJ699 turboSETI hit lists, 0.2-0.8 % in the quiet hydrogen-line coarse channel; much worse in dense RFI bands.
- **Workaround:** `xepoch.analyse` reports the chance-veto probability per ON hit (random frequencies through the same vetoes); quote completeness from injections through the full chain (XEPOCH_RESULT.md), never from detection alone.

### 14. Two dates can be blind to each other — Mitigated
- **Symptom:** dates almost exactly a year apart can have nearly the same Earth velocity toward the star: for Barnard's Star 2016-03-10, +365.25 d and +730.5 d differ by only 0.2-0.5 km/s, so a topocentric tone and a star-frame tone move by < 1 kHz between those dates.
- **Impact:** such a pair cannot tell local interference from a star-fixed tone.
- **Workaround:** `xepoch.analyse` marks a group UNTESTABLE when two of its dates differ by less than 2 × the topocentric veto width in star-frame offset; it never counts as a PASS.
