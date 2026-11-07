#!/usr/bin/env python3
"""Stream Enriquez et al. 2017's All_hits_turbo_seti.csv (10.8 GB, CC BY 4.0)
from stdin into compact npz chunks, never storing the raw CSV.

    curl -s http://blpd0.ssl.berkeley.edu/lband2017/All_hits_turbo_seti.csv \
      | python tools/rfidb_reduce_enriquez.py OUTDIR

Keeps per hit: Freq (MHz, turboSETI start frequency = t=0), DriftRate (Hz/s),
SNR (turboSETI units of the 2017 pipeline), MJD, Source, obs_run.
"""
import os
import sys

import numpy as np
import pandas as pd

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
cols = ["DriftRate", "Freq", "MJD", "SNR", "Source", "obs_run"]
n = 0
for i, ch in enumerate(pd.read_csv(sys.stdin, usecols=cols, chunksize=2_000_000,
                                   dtype={"Source": str, "obs_run": str})):
    ch = ch.dropna(subset=["Freq", "DriftRate", "MJD"])
    np.savez_compressed(os.path.join(out, f"chunk{i:03d}.npz"),
                        freq=ch["Freq"].to_numpy(np.float64),
                        drift=ch["DriftRate"].to_numpy(np.float32),
                        snr=ch["SNR"].to_numpy(np.float32),
                        mjd=ch["MJD"].to_numpy(np.float64),
                        source=ch["Source"].fillna("").to_numpy(str),
                        run=ch["obs_run"].fillna("").to_numpy(str))
    n += len(ch)
    print(f"chunk {i} rows {len(ch)} total {n}", flush=True)
open(os.path.join(out, "DONE"), "w").write(f"{n}\n")
print("done", n)
