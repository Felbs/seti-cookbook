# Nearby-star SETI sweep

> **Evidence level VOID / E0** (single ON files, mid-resolution) — see [EVIDENCE.md](EVIDENCE.md).

> **⚠️ CORRECTION (2026-11-06): the "RFI" verdicts below are VOID.**
> They rest on "all hits share one or two drift rates, so they are instrumental."
> On these mid-resolution (.0002) files one drift step is ~9.8 Hz/s (turboSETI
> reports multiples of ~5.2 Hz/s), wider than the entire range a planet can
> produce (~0.05–4 Hz/s). So **every** signal, real or not, reports the same
> drift. `python drift_audit.py data/star_GJ144.h5` plants ten signals with
> drifts from −3 to +15 Hz/s into the real file: all ten are found, all ten come
> back at ±5.2 Hz/s, and the old rule would have called all ten RFI. The same
> test on the fine-resolution Voyager file gets 7/9 drifts right. These targets
> are **not cleared**; they need a fine-resolution (.0000) ON/OFF re-run.


Swept 10 targets (GBT open data, turboSETI + forensics).

| target | hits | drift rates | verdict |
|---|---|---|---|
| GJ699 | 87 | 2 | RFI (87 hits all at 2 drift rate(s) = instrumental) |
| GJ411 | 55 | 2 | RFI (55 hits all at 2 drift rate(s) = instrumental) |
| GJ887 | - | - | no fine-frequency GBT HDF5 (only HTR) |
| GJ273 | 89 | 2 | RFI (89 hits all at 2 drift rate(s) = instrumental) |
| GJ15A | 14 | 2 | RFI (14 hits all at 2 drift rate(s) = instrumental) |
| GJ71 | 11 | 2 | RFI (11 hits all at 2 drift rate(s) = instrumental) |
| HIP54035 | - | - | no fine-frequency GBT HDF5 (only HTR) |
| HIP57548 | - | - | no fine-frequency GBT HDF5 (only HTR) |
| GJ876 | 18 | 2 | RFI (18 hits all at 2 drift rate(s) = instrumental) |
| GJ581 | - | - | no fine-frequency GBT HDF5 (only HTR) |

**Candidates needing review: none**
## Honest notes
- 6 stars screened (Barnard's/GJ699, GJ411, GJ273, GJ15A, Tau Ceti/GJ71, GJ876);
  4 skipped (archive only had high-time-resolution files, not turboSETI-usable).
- ALL 6 = RFI: every hit pinned to ~2 drift rates (a real planet's signal drifts
  at VARIED rates; a single drift bin = instrumental artifact). Zero candidates.
- Caveat: these are MEDIUM-resolution (.0002) products with coarse drift
  resolution (~9.8 Hz/s) - a SCREENING pass, not a deep search. The hi-spectral
  (.0000, ~3 Hz) files are 15-20 GB each; a deep dive would pull those + run
  find_event over the full ON/OFF cadence. Infra is ready for it.
- The frontier isn't more narrowband screening (all clear so far) - it's the
  NOVEL detectors (cyclostationarity etc, NOVEL_DETECTORS.md) that catch what a
  drift search can't see.
