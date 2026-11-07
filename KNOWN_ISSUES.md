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

### 4. turboSETI reads each coarse channel many times — Open
- **Symptom:** ~18× re-reads of each file; the 77 GB GJ699 cadence took hours.
- **Plan:** single-read coarse-channel search (step 2).

### 5. Archive names match by prefix — Mitigated
- **Symptom:** `?target=HIP2` returns 1,982 files of other stars.
- **Workaround:** `bl_archive.py` (exact names, real OFFs, one session and band).

### 6. Some cadences' OFF scans carry other stars' names — Open
- **Symptom:** Epsilon Eridani (GJ144) has only ON scans under its name; its OFFs are filed under other targets, which the archive cannot be searched for without a name.
- **Impact:** `bl_archive.py` reports no valid cadence rather than building a wrong one.

### 7. Mid-resolution (.0002) drift vetoes are void — Mitigated
- **Cause:** drift quantum ~9.8 Hz/s, wider than any planet's drift (`drift_audit.py`). `star_sweep.py` prints UNRESOLVED.

### 8. find_event's "absent from OFF" means "no hit logged", not "measured absent" — Mitigated
- **Symptom:** all 9 GJ699 filter-3 events (GEO satellites) are present in all three OFF scans at true SNR 29–673; the 1674.9977 MHz track is in both later OFFs.
- **Cause:** an OFF counts as clean when turboSETI logged no hit near it; carrier combs (one hit per neighbourhood) and RFI that changes drift between scans slip through.
- **Workaround:** `cadence_forced.py` measures every scan along the predicted track (plus a free-drift look in OFFs); `cadence_verdict.py --forced`.
