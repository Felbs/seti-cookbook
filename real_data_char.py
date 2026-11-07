#!/usr/bin/env python3
"""real_data_char.py - do the novel detectors behave sanely on REAL RF?

Before trusting cyclo/comb/entropy on sky data, characterize them on real,
KNOWN-type terrestrial captures from the rig: they should FLAG engineered/digital
signals (the signal class turboSETI misses) and stay honest on plain analog
carriers. This is the real-world true-positive / false-alarm study for the
novel-detector suite (task #31), no telescope needed.

  python real_data_char.py            # runs the panel over the known captures
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cyclo
import comb
import entropy

SRC = HERE.parent      # parent folder holding grid-atlas/

# (path, label, kind, is_digital) — kind is the ground-truth modulation
CAPTURES = [
    ("grid-atlas/captures/gps_fix2_20260725.cs16", "GPS L1",   "DSSS spread",   True),
    ("grid-atlas/captures/adsb1090.cs16",          "ADS-B",    "pulsed PPM",    True),
    ("grid-atlas/captures/ais_162.cs16",           "AIS",      "GMSK",          True),
    ("grid-atlas/captures/flex_timed.cs16",        "FLEX pager","4-FSK",        True),
    ("grid-atlas/captures/ft8_7074.cs16",          "FT8",      "8-MFSK",        True),
    ("gr-radiotuna/lab/am_frontier/am_wwfd820_001123Z.cs16", "WWFD 820", "HD-AM OFDM", True),
    ("grid-atlas/captures/chu_7850.cs16",          "CHU",      "AM voice+FSK TC", True),
    ("grid-atlas/captures/dlayer_ref_wwv10.cs16",  "WWV 10",   "AM carrier+tones", False),
]


def load_fs(path):
    for ext in (".json", ".cs16.json"):
        j = Path(str(path) + ext) if not str(path).endswith(".cs16") else Path(str(path)[:-5] + ext)
    # try both <name>.json and <name>.cs16.json
    for cand in (Path(str(path)[:-5] + ".json"), Path(str(path) + ".json")):
        if cand.exists():
            try:
                d = json.load(open(cand))
            except Exception:
                continue
            for k in ("fs_hz", "sample_rate", "fs", "samp_rate"):
                if k in d:
                    return float(d[k])
    return None


def load_iq(path, nmax=6_000_000):
    raw = np.fromfile(path, np.int16, count=2 * nmax).astype(np.float32) / 32768.0
    return (raw[0::2] + 1j * raw[1::2]).astype(np.complex64)


def run():
    print("=" * 92)
    print("NOVEL-DETECTOR real-data characterization (task #31): behavior on KNOWN terrestrial RF")
    print("=" * 92)
    print(f"{'capture':11s} {'modulation':17s} {'truth':7s} | "
          f"{'cyclo(sig)':13s} {'comb(sig)':12s} {'entropy(SFM)':14s} | flags")
    print("-" * 92)
    rows = []
    for rel, label, kind, digital in CAPTURES:
        p = SRC / rel
        if not p.exists():
            print(f"{label:11s} {kind:17s} {'--':7s} | (capture not on disk)")
            continue
        fs = load_fs(p) or 250_000.0
        x = load_iq(str(p))
        if len(x) < 50_000:
            print(f"{label:11s} {kind:17s} | (too short)")
            continue
        ca, cs, _ = cyclo.detect(x, fs)
        ma, ms, _ = comb.detect(x, fs)
        sfm, _ = entropy.detect(x, fs)
        cyc_f = cs >= 8.0
        comb_f = ms >= 10.0
        ent_f = entropy.LO <= sfm <= entropy.HI
        flags = [n for n, f in (("cyclo", cyc_f), ("comb", comb_f), ("entropy", ent_f)) if f]
        rows.append((label, digital, bool(flags), flags, cs, ms, sfm))
        print(f"{label:11s} {kind:17s} {'DIGITAL' if digital else 'analog':7s} | "
              f"{ca:6.0f}Hz {cs:4.1f}  {ma:6.0f}Hz {ms:4.1f}  {sfm:0.3f}        | "
              f"{','.join(flags) if flags else '(silent)'}")
    # negative controls: truly featureless input MUST stay silent (false-alarm check)
    rng = np.random.default_rng(0)
    N = 2_000_000
    noise = ((rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)).astype(np.complex64)
    n = np.arange(N)
    carrier = (5.0 * np.exp(2j * np.pi * 12_345.0 / 250_000.0 * n)
               + 0.1 * noise).astype(np.complex64)
    for lbl, x in (("~noise~", noise), ("~carrier~", carrier)):
        fs = 250_000.0
        ca, cs, _ = cyclo.detect(x, fs); ma, ms, _ = comb.detect(x, fs); sfm, _ = entropy.detect(x, fs)
        flags = [nm for nm, f in (("cyclo", cs >= 8.0), ("comb", ms >= 10.0),
                                  ("entropy", entropy.LO <= sfm <= entropy.HI)) if f]
        print(f"{lbl:11s} {'(control)':17s} {'none':7s} | "
              f"{ca:6.0f}Hz {cs:4.1f}  {ma:6.0f}Hz {ms:4.1f}  {sfm:0.3f}        | "
              f"{','.join(flags) if flags else '(silent - correct)'}")
    print("-" * 92)
    # scorecard: digital signals SHOULD trip >=1 detector; analog is the honesty check
    dig = [r for r in rows if r[1]]
    ana = [r for r in rows if not r[1]]
    dig_hit = sum(1 for r in dig if r[2])
    print(f"DIGITAL/engineered flagged by >=1 detector: {dig_hit}/{len(dig)} "
          f"(these are the class turboSETI misses)")
    print(f"analog references: " + "; ".join(f"{r[0]} -> {','.join(r[3]) or 'silent'}" for r in ana))
    print("=" * 92)
    return rows


if __name__ == "__main__":
    run()
