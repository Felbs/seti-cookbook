# Data notice

The MIT licence in [LICENSE](LICENSE) covers this repository's **code only**.
Data tables shipped here are derived from other people's published data and keep
their owners' terms. If you reuse them, cite the sources below.

| file | contents | source | terms |
|---|---|---|---|
| `rfi_seed.csv.gz` | 9,584 rows (frequency, drift, SNR, source, file id) from the Price et al. 2020 event tables | Price, D. C. et al. 2020, AJ 159, 86, [doi:10.3847/1538-3881/ab65f1](https://doi.org/10.3847/1538-3881/ab65f1); `lband2019_events.csv`, `pks10cm2019_events.csv` from [seti.berkeley.edu/listen2019/opendata.html](https://seti.berkeley.edu/listen2019/opendata.html) | No explicit licence stated for these files (checked 2026-11-07). Released by Breakthrough Listen as public data for reanalysis; this paper's Zenodo deposit ([10.5281/zenodo.7600403](https://doi.org/10.5281/zenodo.7600403)) and the 2017 release are CC BY 4.0. Redistributed with attribution; if the owners object, `python rfi_db.py build` regenerates the seed without them. |
| `rfi_seed.csv.gz` | our own turboSETI hits on Breakthrough Listen open data (GJ699 GBT; Parkes UWL blc1 windows) | Breakthrough Listen open data archive; Lebofsky, M. et al. 2019, PASP 131, 124505, [doi:10.1088/1538-3873/ab3e82](https://doi.org/10.1088/1538-3873/ab3e82) | our measurements of public data; cite the archive |
| `targets.json` | Breakthrough Listen nearby-star target list with SIMBAD values | Isaacson, H. et al. 2017, PASP 129, 054501 (VizieR J/PASP/129/E4501); SIMBAD (CDS, Strasbourg) | catalogue facts with attribution |
| not shipped | Enriquez et al. 2017 hit tables (built locally by `rfi_db.py`) | Enriquez, J. E. et al. 2017, ApJ 849, 104, [seti.berkeley.edu/lband2017](http://seti.berkeley.edu/lband2017/) | CC BY 4.0 |

Breakthrough Listen is managed by the Breakthrough Initiatives, sponsored by the
Breakthrough Prize Foundation.
