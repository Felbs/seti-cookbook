#!/usr/bin/env python3
"""tools/fetch_coarse.py - read whole TRUE coarse channels from a BL HDF5 file
(local path or URL, HTTP range reads via fsspec) and save each as an npz slice.

The npz has the same keys as hi_starframe.py's (data, ch0, fch1, foff, tsamp,
tstart, ra_h, dec_deg, kind) plus nchans and coarse, so recovery.py,
fastsearch_recovery.py and cadence_forced.py read it. The coarse width comes
from coarse_align.nfine_per_coarse (1,033,216 for the 2016 GJ699 fine-res files).
A remote read transfers the compressed HDF5 chunks the channel touches (2016 GBT
fine-res: 2^20-channel chunks, ~0.75 compression, so ~100 MB per coarse channel).

    python tools/fetch_coarse.py URL_OR_PATH --coarse 112,113 --out-dir data/coarse --kind ON
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import coarse_align  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("src")
    ap.add_argument("--coarse", required=True, help="comma list of coarse-channel indices")
    ap.add_argument("--out-dir", default="data/coarse")
    ap.add_argument("--kind", default="ON", choices=["ON", "OFF"])
    a = ap.parse_args()
    import h5py
    import hdf5plugin  # noqa: F401
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    name = a.src.rstrip("/").split("/")[-1]
    if a.src.startswith("http"):
        import fsspec
        fh = fsspec.open(a.src, "rb", block_size=8 * 2 ** 20).open()
    else:
        fh = open(a.src, "rb")
    with fh, h5py.File(fh, "r") as h:
        d = h["data"]
        at = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in d.attrs.items()}
        nchans = int(at["nchans"])
        nfpc = coarse_align.nfine_per_coarse(float(at["foff"]), nchans, int(at.get("telescope_id", 6)))
        for cc in (int(c) for c in a.coarse.split(",")):
            fn = out / f"cc{cc:03d}_{name.replace('.h5', '.npz')}"
            if fn.exists():
                print(f"have {fn}")
                continue
            c0 = cc * nfpc
            x = d[:, 0, c0:c0 + nfpc].astype(np.float32)
            np.savez(fn, data=x, ch0=c0, fch1=at["fch1"], foff=at["foff"], tsamp=at["tsamp"],
                     tstart=at["tstart"], ra_h=at["src_raj"], dec_deg=at["src_dej"], kind=a.kind,
                     nchans=nchans, coarse=cc)
            f0 = float(at["fch1"]) + float(at["foff"]) * c0
            print(f"coarse {cc}: ch {c0}..{c0 + nfpc} = {f0:.4f} to {f0 + float(at['foff']) * nfpc:.4f} MHz, "
                  f"{x.nbytes / 1e6:.0f} MB -> {fn}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
