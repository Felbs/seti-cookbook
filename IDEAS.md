# seti-cookbook — novel search ideas (2026-07-25, user + Fable brainstorm)
1. SIDEBAND-PAIR SEARCH (the star idea): turboSETI finds isolated lines; nobody
   post-searches the hit list for SYMMETRIC PAIRS/TRIPLETS sharing a common
   drift rate (f0 +/- delta, same drift) = a MODULATED transmitter — exactly
   Voyager's carrier+2-sidebands signature, generalized. Buildable today:
   cross-correlate .dat hits within each file for drift-matched symmetry.
   "An alien sideband humans wouldn't know to look for" — found by structure,
   not by strength.
2. DRIFT-CURVATURE PRIOR: RFI drifts ~linearly or not at all; a rotating/
   orbiting planet's Doppler is SINUSOIDAL -> fit curvature (2nd derivative)
   across long/multi observations; curvature consistent with plausible orbital
   dynamics = high-value candidate class no linear search ranks.
3. ANTI-CADENCE BEACON FILTER: find_event kills anything present in OFF
   pointings — but an isotropic mega-beacon WOULD appear in OFFs via sidelobes.
   Search instead for ON/OFF strength RATIO matching the beam-gain ratio.
4. COMPLEXITY SCAN: sliding-tile compression ratio (zlib on quantized
   spectrogram tiles): pure noise = incompressible, pure carrier = trivial;
   INTERMEDIATE complexity = modulation/structure. Modulation-agnostic detector
   for signals that aren't narrowband at all.
