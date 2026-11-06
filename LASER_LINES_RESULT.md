# Hunt result: laser lines in 199 Sloan star spectra — 2026-11-06

**Idea.** A laser aimed at us shows up in a star's spectrum as a needle: a
line as narrow as the spectrograph allows, at a wavelength where nothing
natural belongs, steady from one exposure to the next (Townes 1983). The
Sloan Digital Sky Survey (DR17) keeps every separate exposure of every star —
often from several nights — plus the sky it subtracted. `laser_lines.py`.

**Sample.** 200 K and M dwarfs (SDSS DR17 legacy, signal-to-noise > 15, on
plates with many exposures), ~150 MB. 199 searched; one had no comparison
stars of its type.

**What a candidate had to survive** (pre-registered thresholds in the script):
unresolved; in ≥ 80 % of exposures; steady flux; not a sky line; not a known
stellar line; not in the comparison spectrum of the star's own type; not in a
stellar absorption line; the star model fits there; does not recur in ≥ 3 %
of the sample.

| stage | peaks |
|---|---|
| coadd matched-filter SNR ≥ 7 | 431 |
| star model does not fit there | 369 vetoes |
| wider than the instrument (resolved) | 315 |
| at a stellar absorption line (shallower than comparison stars) | 273 |
| known line / sky line / not every exposure / unsteady / recurring | 46 / 35 / 36 / 15 / 18 |
| **survivors** | **3** |

(A peak can collect several vetoes.)

| star | wavelength | SNR | exposures | note |
|---|---|---|---|---|
| spec-0271-51883-0440 | 7181.25 Å | 19.6 | 4/4 | inside Earth's 7160–7340 Å water band |
| spec-0278-51900-0030 | 8113.35 Å | 8.8 | 3/3 | inside Earth's 8100–8350 Å water band |
| spec-0296-51984-0232 | 8090.96 Å | 8.8 | 4/4 | just blueward of that band |

Telluric (Earth-atmosphere) bands were **not** in the pre-registered veto
list; that note was added after seeing the survivors. Residuals of telluric
correction are the leading explanation for two, and probably the third. All
three are rung 1, "not yet explained". None is a laser until another
spectrum of the same star shows the same line.

**Sensitivity, measured by planting lines:** inside the searchable 60 % of
each spectrum, planted unresolved lines are recovered ~89 % of the time at
coadd SNR ≥ 12 (≈ 12 × 10⁻¹⁷ erg s⁻¹ cm⁻²), 65 % at SNR 8, ~0 % at SNR 5.
Recovery plateaus just under 90 %: sky lines and exposure-to-exposure flux
scatter take the rest. The other 40 % of a K/M dwarf spectrum sits in its own
absorption lines or where the comparison spectrum fails, and this method
cannot search there.

**Bugs found on the way (fixed, kept as lessons):** Sloan stores blue-camera
frames in *descending* wavelength order — before sorting, the "in every
exposure" test looked at the wrong pixels and lost half of all planted lines.
M dwarfs' molecular bands and shallower-than-usual metal lines (Mg b, Ca I)
both mimic narrow emission until compared with same-type stars.
