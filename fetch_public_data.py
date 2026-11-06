#!/usr/bin/env python3
"""fetch_public_data.py - get REAL telescope data showing the natural phenomena
in SETI_HISTORY.md. Fetch scripts only: this repo never ships anyone else's data.

Everything below is free, public, and produced by other people's telescopes. Each
entry says what PHENOMENON it shows and roughly how big it is, so you can decide
before downloading. Files land in data/ , which is gitignored.

  python fetch_public_data.py list                 # what is available and why
  python fetch_public_data.py callisto             # solar radio bursts (small!)
  python fetch_public_data.py voyager              # the standing calibration
  python fetch_public_data.py bl --target GJ699    # a nearby-star SETI pointing

Credit where it is due - if you publish anything based on these, cite the source:
  * Breakthrough Listen open data (Lebofsky et al. 2019, PASP 131:124505;
    Price et al. 2020, AJ 159:86). Berkeley SETI Research Center.
  * e-CALLISTO, the international network of solar radio spectrometers
    (Benz, Monstein & Meyer 2005, Solar Physics 226:143), FHNW Switzerland.
  * CHIME/FRB Catalog 1 (CHIME/FRB Collaboration 2021, ApJS 257:59).
  * turboSETI (Enriquez & Price 2019) and blimpy (Price et al. 2019) for reading
    and searching the BL products.
"""
import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

CATALOG = {
    "voyager": dict(
        shows="a REAL confirmed technosignature at interstellar distance: "
              "Voyager 1's X-band carrier + both telemetry sidebands, drifting "
              "-0.38 Hz/s. seti-cookbook's standing calibration.",
        size="~50 MB",
        how="direct URL from the BL open-data archive",
        url="http://blpd0.ssl.berkeley.edu/Voyager_data/"
            "Voyager1.single_coarse.fine_res.h5",
        out="Voyager1.single_coarse.fine_res.h5"),
    "callisto": dict(
        shows="SOLAR RADIO BURSTS in a real waterfall - Type III bursts sweep "
              "hundreds of MHz downward in seconds, Type II drift slowly. The "
              "cheapest real dynamic spectrum of a violent natural radio source, "
              "and the closest thing to 'aliens' shape a hobbyist can see daily.",
        size="~1-3 MB per 15-minute FITS file",
        how="e-CALLISTO archive, browse by date: "
            "http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/",
        url=None,
        out=None,
        note="Pick a date with a known flare (e.g. 2011-08-09, an X6.9) and a "
             "station in daylight; file names look like "
             "BLEN7M_20110809_080000_59.fit.gz. Needs astropy to read FITS."),
    "bl": dict(
        shows="a real SETI pointing at a nearby star: GBT L- or S-band, ~3 Hz "
              "channels, and (if the band covers 1420 MHz) the galaxy's own "
              "21 cm hydrogen line - run natural_signals.py on it.",
        size="~200-300 MB per file",
        how="queries the BL archive API by target name, downloads the smallest "
            "fine-frequency HDF5. This is exactly what star_sweep.py does.",
        url="http://seti.berkeley.edu/opendata/api/query-files",
        out=None),
    "frb": dict(
        shows="FAST RADIO BURSTS with their dispersion sweep. Breakthrough Listen "
              "recorded 21 bursts from the repeater FRB 121102 at 4-8 GHz with the "
              "GBT (Gajjar et al. 2018, ApJ 863:2); the data is public. The "
              "CHIME/FRB Catalog 1 (ApJS 257:59) publishes waterfalls for 536 "
              "bursts and is the friendlier starting point.",
        size="BL FRB 121102 raw products are LARGE (many GB); CHIME catalogue "
             "waterfalls are small",
        how="browse https://seti.berkeley.edu/frb-machine/ (BL FRB data) or "
            "https://www.chime-frb.ca/catalog (CHIME/FRB Catalog 1)",
        url=None, out=None),
    "pulsar": dict(
        shows="PULSARS - periodic dispersed pulses. The classic bright targets are "
              "PSR B0329+54 (P=0.7145 s, DM 26.8) and Vela/B0833-45 (P=0.089 s, "
              "DM 68). Fold with waterfall.py fold --period.",
        size="varies",
        how="the BL archive has pulsar pointings (query the API by pulsar name, "
            "e.g. target=B0329+54); the EPN Database of Pulsar Profiles "
            "(https://psrweb.jb.man.ac.uk/epndb/) publishes profiles",
        url=None, out=None),
}


def cmd_list():
    print(__doc__.split("Credit")[0])
    for k, v in CATALOG.items():
        print("=" * 74)
        print(f"{k}   [{v['size']}]")
        print(f"  shows: {v['shows']}")
        print(f"  how:   {v['how']}")
        if v.get("note"):
            print(f"  note:  {v['note']}")
    print("=" * 74)
    print("data/ is gitignored - nothing downloaded here is redistributed by "
          "this repo.")


def _curl(url, out):
    DATA.mkdir(exist_ok=True)
    dest = DATA / out
    if dest.exists() and dest.stat().st_size > 1e6:
        print(f"already have {dest} ({dest.stat().st_size/1e6:.0f} MB)")
        return dest
    print(f"fetching {url}\n     -> {dest}")
    subprocess.run(["curl", "-fL", "--progress-bar", "-o", str(dest), url],
                   check=False)
    if dest.exists():
        print(f"got {dest.stat().st_size/1e6:.1f} MB")
        return dest
    print("download failed")
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", nargs="?", default="list", choices=list(CATALOG) + ["list"])
    ap.add_argument("--target", default="GJ699", help="BL target name")
    args = ap.parse_args()
    if args.what == "list":
        cmd_list()
        return 0
    entry = CATALOG[args.what]
    if args.what == "bl":
        print("delegating to star_sweep.py, which queries the archive, downloads "
              "the smallest fine-frequency HDF5 and runs the search:")
        print(f"    python star_sweep.py {args.target}")
        return 0
    if not entry.get("url"):
        print(f"{args.what}: no single direct URL - by design.\n")
        print(f"  shows: {entry['shows']}\n  how:   {entry['how']}")
        if entry.get("note"):
            print(f"  note:  {entry['note']}")
        return 0
    _curl(entry["url"], entry["out"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
