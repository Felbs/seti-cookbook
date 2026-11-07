# Hunt result: the hydrogen line in Barnard's Star's own frame — 2026-11-06

> **Evidence level E2** — see [EVIDENCE.md](EVIDENCE.md).

**Idea.** If someone at Barnard's Star (GJ699, 5.96 ly) runs a beacon at the
hydrogen line — the frequency every radio astronomer knows — tuned in *their*
star's rest frame, then the star's radial velocity (−110.5 km/s) plus Earth's
motion predict exactly where it lands in our data. That shrinks the haystack
from 750 MHz to ±2.5 kHz.

**Data.** The public Green Bank fine-resolution (2.84 Hz × 18 s) ON/OFF
cadence of 2016-03-10, three 288 s scans on the star and three off it. We read
only the one 2.93 MHz coarse channel that holds the prediction, straight from
the archive over HTTP (64 MB of each 12.8 GB file). `python hi_starframe.py
fetch && python hi_starframe.py run`.

**Where it should land.** 1421.052 MHz — not 1420.93: that March, Earth's
orbit was carrying Green Bank toward the star at ~26 km/s on top of the star's
own approach. A second, pre-declared window uses the Solar-System barycentre
frame (1420.528 MHz). Earth's drift on the track (−0.114 Hz/s) comes from
`earth_term.py`; the search allows ±0.3 Hz/s around it for the star's planet.

**Pre-registered** before looking: threshold = max(SNR 10, 99th percentile of
the same search at ~900 empty windows) → 10.0. Detection = above threshold in
all three ONs at the same frequency and in no OFF.

| frame | predicted | best SNR in the 6 scans | verdict |
|---|---|---|---|
| star rest frame | 1421.052026 MHz | 4.1, 4.8, 4.8, 4.3, 5.0, 4.7 | **null** |
| barycentre | 1420.528481 MHz | 4.8, 4.6, 5.2, 4.3, 4.3, 5.2 | **null** |

**How sensitive was it — measured, not assumed.** Tones planted at 30 empty
spots in each ON scan come back at 68 % of their nominal strength (they smear
across neighbouring channels), and are recovered 90 % of the time at SNR 12.
With BL's published SEFD (10 Jy) that is:

> **no steady, narrowband (~1 Hz) hydrogen-line beacon in Barnard's Star's
> frame brighter than ~12 Jy — about 5 GW aimed at us — was on during those
> three 5-minute scans in March 2016.**
>
> The 90 % is per scan; the detection rule needs all three ONs (≈ 78 % at
> measured SNR 12), so the limit for that rule is ~12-16 Jy (~5-6 GW).

That is 1/2000 of the Arecibo planetary radar. Caveat (2026-11-06 peer review): the planted
tones were not smeared within each 18 s spectrum. At the injected drift (Earth's −0.114 Hz/s,
< 1 channel per spectrum) that barely matters, but for a beacon whose OWN drift adds more than
~±0.04 Hz/s the real sensitivity is worse than stated. It is not a clearance of
Barnard's Star: one night, one frequency convention, a 15-minute window.

**Amendment, stated openly.** The plan said "an SNR-8 tone must be recovered".
SNR 8 is below the threshold, so that test could not show detection works; it
was replaced by the completeness curve above. That changes the sensitivity
claim, not what counted as a detection.
