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
