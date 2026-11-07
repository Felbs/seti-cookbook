#!/usr/bin/env python3
"""rfidb_blc1_fetch.py - the blc1 control for rfi_db.py, from the PUBLIC Parkes data.

blc1 (Smith et al. 2021, Sheikh et al. 2021, Nature Astronomy 5) was a 982.002 MHz
narrowband signal seen toward Proxima Centauri on 2019-04-29 that passed ON/OFF
and was later explained as an intermodulation product of local interference.
The Parkes UWL filterbank files used by those papers are public (BL open-data
archive, paperName=Sheikh21a; 14 GB per 30-min scan). We do NOT download them:
filterbank is header + raw spectra, so HTTP Range requests fetch only a few kHz
around each frequency of interest from every spectrum (~1 MB per scan per window).

For each file and window we run a small brute-force de-drift search (drift grid
+-DRIFT_MAX Hz/s, step one channel per scan) and keep hits with TRUE radiometer
SNR >= SNR_MIN (time-summed, de-drifted power over the channel-to-channel noise
with a 65-channel running median removed - the same definition as
inject.radiometer_noise). Hits are written in rfi_db's schema.

    python tools/rfidb_blc1_fetch.py list   [--days 58598,58602]
    python tools/rfidb_blc1_fetch.py fetch  OUT.csv [--days ...] [--windows ...]

Windows (MHz, +-HALF_KHZ): blc1 at 982.0026 and the three lookalikes Sheikh et
al. relate to it by multiples of a 2.000004 MHz clock (~712, ~856, ~1062 MHz).
"""
import argparse
import csv
import json
import re
import struct
import sys
import time
import urllib.parse
import urllib.request

import numpy as np

API = "http://seti.berkeley.edu/opendata/api/query-files"
CLOCK = 2.000004                      # MHz, Sheikh et al. 2021 section 4.4.2
BLC1 = 982.0026                       # MHz, Sheikh et al. 2021 Supp. Table 1 (982.00241-982.00281)
WINDOWS = {"blc1": BLC1,
           "lk712": BLC1 - 135 * CLOCK,     # 712.0019
           "lk856": BLC1 - 63 * CLOCK,      # 856.0021
           "lk1062": BLC1 + 40 * CLOCK}     # 1062.0026
HALF_KHZ = 5.0
DRIFT_MAX = 0.5
SNR_MIN = 10.0

INT_KEYS = {"telescope_id", "machine_id", "data_type", "barycentric", "pulsarcentric",
            "nbits", "nchans", "nifs", "nbeams", "ibeam", "nsamples"}
DBL_KEYS = {"tstart", "tsamp", "fch1", "foff", "refdm", "az_start", "za_start",
            "src_raj", "src_dej", "period"}
STR_KEYS = {"source_name", "rawdatafile"}


def query_paper(paper="Sheikh21a"):
    p = urllib.parse.urlencode({"paperName": paper, "target": ""})
    with urllib.request.urlopen(API + "?" + p, timeout=120) as r:
        return json.loads(r.read())["data"]


def get_range(url, a, b, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status != 206:
                    raise IOError(f"server ignored Range (status {r.status})")
                return r.read()
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


def parse_header(buf):
    def rd_str(i):
        n = struct.unpack("<i", buf[i:i + 4])[0]
        return buf[i + 4:i + 4 + n].decode("latin1"), i + 4 + n
    s, i = rd_str(0)
    if s != "HEADER_START":
        raise ValueError("not a sigproc filterbank")
    h = {}
    while True:
        key, i = rd_str(i)
        if key == "HEADER_END":
            h["hdr_len"] = i
            return h
        if key in INT_KEYS:
            h[key] = struct.unpack("<i", buf[i:i + 4])[0]; i += 4
        elif key in DBL_KEYS:
            h[key] = struct.unpack("<d", buf[i:i + 8])[0]; i += 8
        elif key in STR_KEYS:
            h[key], i = rd_str(i)
        else:
            raise ValueError(f"unknown header key {key}")


def read_window(url, size, f_mhz, half_khz=HALF_KHZ):
    h = parse_header(get_range(url, 0, 4095))
    if h.get("nbits") != 32 or h.get("nifs", 1) != 1:
        raise ValueError(f"unsupported nbits/nifs {h.get('nbits')}/{h.get('nifs')}")
    nc, foff, fch1 = h["nchans"], h["foff"], h["fch1"]
    c = (f_mhz - fch1) / foff
    hw = half_khz * 1e-3 / abs(foff)
    i0, i1 = int(np.floor(c - hw)), int(np.ceil(c + hw))
    if i0 < 0 or i1 >= nc:
        return None, h
    row = nc * 4
    nt = (size - h["hdr_len"]) // row
    out = np.empty((nt, i1 - i0), np.float32)
    for t in range(nt):
        a = h["hdr_len"] + t * row + i0 * 4
        out[t] = np.frombuffer(get_range(url, a, a + (i1 - i0) * 4 - 1), "<f4")
    h["i0"], h["nt"] = i0, nt
    h["f0_mhz"] = fch1 + foff * i0          # frequency of window channel 0
    return out, h


def robust_std(x):
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))


def dedrift_search(data, foff_hz, tsamp, drift_max=DRIFT_MAX, snr_min=SNR_MIN, guard=8):
    """Brute-force de-drift. Returns list of (chan_at_t0, drift_hz_s, true_snr)."""
    from scipy.ndimage import median_filter
    nt, nw = data.shape
    T = nt * tsamp
    step = abs(foff_hz) / T
    drifts = np.arange(-drift_max, drift_max + step / 2, step)
    norm = data / np.median(data, axis=1, keepdims=True)
    t = np.arange(nt) * tsamp
    best = np.full(nw, -np.inf)
    bestd = np.zeros(nw)
    pad = int(np.ceil(drift_max * T / abs(foff_hz))) + 1
    for d in drifts:
        sh = np.round(d * t / foff_hz).astype(int)        # channel shift per spectrum
        s = np.zeros(nw)
        for k in range(nt):
            s += np.roll(norm[k], -sh[k])
        s[:pad] = np.nan
        s[-pad:] = np.nan
        ok = np.isfinite(s)
        r = np.full(nw, np.nan)
        base = median_filter(np.where(ok, s, np.nanmedian(s)), size=65)
        r[ok] = s[ok] - base[ok]
        z = r / robust_std(r[ok])
        better = np.nan_to_num(z, nan=-np.inf) > best
        best[better] = z[better]
        bestd[better] = d
    hits = []
    order = np.argsort(best)[::-1]
    taken = np.zeros(nw, bool)
    for i in order:
        if best[i] < snr_min:
            break
        if taken[max(0, i - guard):i + guard + 1].any():
            continue
        taken[i] = True
        hits.append((int(i), float(bestd[i]), float(best[i])))
    return hits


def select(rows, days, centres):
    out = []
    for x in rows:
        if int(x["mjd"]) in days and x["file_type"] == "filterbank" \
                and any(abs(x["center_freq"] - c) < 1 for c in centres):
            out.append(x)
    return sorted(out, key=lambda x: x["mjd"])


def subband_centre(f_mhz):
    # UWL BL subbands: 128 MHz wide, centres 768, 896, 1024, ... (704-4032 MHz)
    return 768 + 128 * round((f_mhz - 768) / 128)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["list", "fetch"])
    ap.add_argument("out", nargs="?")
    ap.add_argument("--days", default="58598,58602")
    ap.add_argument("--windows", default=",".join(WINDOWS))
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--min-size-mb", type=float, default=500,
                    help="skip short (<~3 min) scans")
    a = ap.parse_args()
    days = {int(d) for d in a.days.split(",")}
    wins = {k: WINDOWS[k] for k in a.windows.split(",")}
    rows = query_paper()
    sel = select(rows, days, {subband_centre(f) for f in wins.values()})
    sel = [x for x in sel if int(x["size"]) >= a.min_size_mb * 1e6]
    if a.cmd == "list":
        for x in sel:
            print(f"{x['mjd']:.4f} {x['target']:16s} {x['center_freq']:7.0f} "
                  f"{int(x['size'])//10**6:6d} MB {x['url']}")
        print(len(sel), "files")
        return 0
    w = csv.writer(open(a.out, "w", newline=""))
    w.writerow(["freq_mhz", "drift_hz_s", "snr", "snr_units", "target", "session", "mjd",
                "telescope", "product", "res_hz", "tobs_s", "drift_valid", "source", "window"])
    jobs = [(x, name, f) for x in sel for name, f in wins.items()
            if abs(x["center_freq"] - subband_centre(f)) <= 1]

    def one(job):
        x, name, f = job
        try:
            data, h = read_window(x["url"], int(x["size"]), f)
        except Exception as e:                           # noqa: BLE001
            return job, None, str(e)
        if data is None or data.shape[0] < 4:
            return job, None, "window outside file"
        return job, h, dedrift_search(data, h["foff"] * 1e6, h["tsamp"])

    from concurrent.futures import ThreadPoolExecutor
    n = 0
    with ThreadPoolExecutor(a.threads) as ex:
        for (x, name, f), h, hits in ex.map(one, jobs):
            if h is None:
                print("  skip", x["url"].rsplit("/", 1)[-1], name, hits, file=sys.stderr, flush=True)
                continue
            foff_hz = h["foff"] * 1e6
            T = h["nt"] * h["tsamp"]
            for ch, d, snr in hits:
                fm = h["f0_mhz"] + ch * h["foff"]
                w.writerow([f"{fm:.6f}", f"{d:+.5f}", f"{snr:.1f}", "true", x["target"],
                            f"Parkes:{int(x['mjd'])}", f"{x['mjd']:.6f}", "Parkes",
                            f"UWL fine {abs(foff_hz):.2f}Hz", f"{abs(foff_hz):.3f}", f"{T:.1f}", 1,
                            "BL open data " + x["url"].rsplit("/", 1)[-1], name])
                n += 1
            print(f"{x['mjd']:.4f} {x['target']:14s} {name:7s} nt={h['nt']:4d} hits={len(hits)} "
                  + " ".join(f"{h['f0_mhz'] + c * h['foff']:.6f}@{d:+.3f}/{s:.0f}" for c, d, s in hits[:4]),
                  flush=True)
    print("rows", n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
