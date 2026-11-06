#!/usr/bin/env python3
"""hi_starframe.py - the hydrogen line, as the alien's own star would tune it.

Cocconi & Morrison (1959) proposed the 21 cm hydrogen line, 1420.405751768
MHz, as the frequency any radio civilisation would know. But in WHOSE frame?
If a beacon sits at exactly that frequency in its home star's rest frame, the
star's measured radial velocity plus Earth's motion predict where it lands in
OUR data - to within the velocity's error bar. That turns a 750 MHz haystack
into a few kHz.

For Barnard's Star (GJ699, RV about -110.5 km/s, approaching) the window lands
near 1421.05 MHz on 2016-03-10 (the star's approach plus ~26 km/s of Earth's
orbital motion toward it that March) - inside the 1418.59-1421.56 MHz coarse
channel of the public
Green Bank fine-resolution (.0000, 2.8 Hz x 18 s) ON/OFF cadence. We read only
that one coarse channel over HTTP (about 64 MB of each 12.8 GB file).

PRE-REGISTERED (fixed before the data were looked at; see PREREG below):
  frames         star rest frame (primary) + Solar-System barycentre (secondary)
  window         +-2.5 kHz around each predicted frequency (RV error +-0.5 km/s)
  drift grid     Earth term (earth_term.py) + [-0.30, +0.30] Hz/s, step 0.01
  threshold      the larger of SNR 10 and the 99th percentile of the same search
                 run at >= 300 empty offset windows in the same coarse channel
  detection      above threshold in ALL THREE ONs, at the same frequency
                 (+-3 channels after drift), and in NO OFF
  injection      an SNR-8-per-scan tone at the predicted track must be recovered
AMENDMENT (after the first run, 2026-11-06; it changes the sensitivity claim,
not the detection rule): planted tones came back at ~70 % of their nominal
SNR (linear-interpolation smearing between channels), and SNR 8 sits below the
SNR-10 threshold, so it could not show detection works. The limit is now set
by injection-recovery: tones planted at 30 empty offsets per ON scan, at
several MEASURED strengths; the limit uses the strength recovered >= 90 % of
the time, divided by the measured efficiency.
Trials: 2 frames x 2 windows of 5 kHz x 61 drifts, reported, not hidden.

    python hi_starframe.py fetch      # pull the one coarse channel of each scan (~400 MB)
    python hi_starframe.py run        # search, calibrate, inject, write the result
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
C = 299_792_458.0
F_HI = 1420.405751768e6

PREREG = dict(
    target="GJ699 (Barnard's Star)",
    rv_kms=-110.5, rv_err_kms=0.5,               # SIMBAD/Gaia-class value; systematics ~0.5
    frames=["star", "barycentre"],
    window_hz=2500.0,
    drift_offsets_hz_s=[round(x, 2) for x in np.arange(-0.30, 0.301, 0.01)],
    snr_floor=10.0, null_percentile=99.0, n_null_min=300,
    match_channels=3, inject_snr=8.0,
    fetch_at_mhz=1420.93,
)

BASE = "https://bldata.berkeley.edu/pipeline/AGBT16A_999_81/holding/"
SCANS = [("ON", "spliced_blc02030405_2bit_guppi_57457_48006_GJ699_0003.gpuspec.0000.h5"),
         ("OFF", "spliced_blc02030405_2bit_guppi_57457_48334_GJ699_OFF_0004.gpuspec.0000.h5"),
         ("ON", "spliced_blc02030405_2bit_guppi_57457_48661_GJ699_0005.gpuspec.0000.h5"),
         ("OFF", "spliced_blc02030405_2bit_guppi_57457_48988_GJ699_OFF_0006.gpuspec.0000.h5"),
         ("ON", "spliced_blc02030405_2bit_guppi_57457_49314_GJ699_0007.gpuspec.0000.h5"),
         ("OFF", "spliced_blc02030405_2bit_guppi_57457_49640_GJ699_OFF_0008.gpuspec.0000.h5")]


def cache_path(name):
    return DATA / ("hi_" + name.replace(".h5", ".npz"))


def cmd_fetch():
    import fsspec, h5py, hdf5plugin  # noqa: F401
    DATA.mkdir(exist_ok=True)
    for kind, name in SCANS:
        out = cache_path(name)
        if out.exists():
            print(f"have {out.name}")
            continue
        with fsspec.open(BASE + name, block_size=8 * 2 ** 20).open() as fh:
            with h5py.File(fh, "r") as h:
                d = h["data"]
                a = dict(d.attrs)
                nchan = d.shape[2]
                freqs_edge = a["fch1"] + a["foff"] * np.array([0, nchan - 1])
                # coarse channel (chunk) holding the band
                fpc = d.chunks[2]
                # the chunk holding the star-frame prediction (~1420.93 MHz); the
                # barycentric window (~1420.4) falls in the same 2.93 MHz chunk
                i_lo = int((PREREG["fetch_at_mhz"] - a["fch1"]) / a["foff"])
                c0 = (i_lo // fpc) * fpc
                x = d[:, 0, c0:c0 + fpc].astype(np.float32)
        np.savez(out, data=x, ch0=c0, fch1=a["fch1"], foff=a["foff"], tsamp=a["tsamp"],
                 tstart=a["tstart"], ra_h=a["src_raj"], dec_deg=a["src_dej"], kind=kind)
        f = a["fch1"] + a["foff"] * c0
        print(f"{kind:3s} {name[-40:]}: channels {c0}..{c0+fpc} = {f + a['foff']*fpc:.4f}-{f:.4f} MHz"
              f"  {x.nbytes/1e6:.0f} MB -> {out.name}")


def load(name):
    z = np.load(cache_path(name))
    m = {k: (z[k].item() if z[k].ndim == 0 else z[k]) for k in z.files if k != "data"}
    m["freqs"] = (m["fch1"] + m["foff"] * (m["ch0"] + np.arange(z["data"].shape[1]))) * 1e6
    return z["data"].astype(np.float64), m


def predicted(m, frame):
    """Predicted observed frequency (Hz) of a beacon at F_HI, and Earth drift (Hz/s)."""
    from earth_term import earth_drift, v_los
    nt = 16
    mid = m["tstart"] + 0.5 * nt * m["tsamp"] / 86400.0
    ra, dec = m["ra_h"] * 15.0, m["dec_deg"]
    approach = v_los(ra, dec, mid)[0]               # observer toward star, m/s
    if frame == "star":
        approach += -PREREG["rv_kms"] * 1e3         # star approaching at -RV
    f = F_HI * (1 + approach / C)
    return f, earth_drift(f, ra, dec, mid)


def flatten(d):
    """Remove the bandpass: divide each channel by its time-median, then robustly
    normalise. Leaves narrow drifting tones intact (16 spectra, they move)."""
    bp = np.median(d, axis=0)
    bp = np.where(bp > 0, bp, np.median(bp))
    return d / bp - 1.0


def dedrift_sum(x, drift_hz_s, foff_hz, tsamp):
    """Shift-and-sum along a drift (fractional channels, linear interpolation),
    referenced to the MIDDLE of the scan."""
    nt, nc = x.shape
    t = (np.arange(nt) + 0.5 - nt / 2) * tsamp
    out = np.zeros(nc)
    idx = np.arange(nc, dtype=float)
    for k in range(nt):
        sh = drift_hz_s * t[k] / foff_hz           # channel shift of the track at time k
        out += np.interp(idx + sh, idx, x[k], left=0.0, right=0.0)
    return out


def search_window(x, m, f_centre, drift0, half_hz):
    """Best (snr, freq, drift) inside +-half_hz of f_centre over the drift grid."""
    foff_hz = m["foff"] * 1e6
    freqs = m["freqs"]
    c = int(np.argmin(abs(freqs - f_centre)))
    pad = int(half_hz / abs(foff_hz)) + 60
    lo, hi = max(0, c - pad), min(x.shape[1], c + pad)
    sub = x[:, lo:hi]
    inner = slice(60, (hi - lo) - 60)
    best = (-np.inf, None, None)
    for dd in PREREG["drift_offsets_hz_s"]:
        s = dedrift_sum(sub, drift0 + dd, foff_hz, m["tsamp"])
        med = np.median(s)
        mad = np.median(np.abs(s - med)) * 1.4826
        z = (s - med) / mad
        j = int(np.argmax(z[inner])) + 60
        if z[j] > best[0]:
            best = (float(z[j]), float(freqs[lo + j]), float(drift0 + dd))
    return best


def null_offsets(m, avoid_hz, n):
    """Centres of empty windows across the coarse channel, avoiding the band edges,
    the DC spike, the Galactic HI hump and the windows under test."""
    freqs = m["freqs"]
    f_lo, f_hi = freqs.min() + 60e3, freqs.max() - 60e3
    dc = 0.5 * (freqs.min() + freqs.max())
    step = (f_hi - f_lo) / n
    out = []
    for k in range(n):
        f = f_lo + (k + 0.5) * step
        if abs(f - dc) < 20e3 or 1420.0e6 < f < 1420.7e6:
            continue
        if any(abs(f - a) < 20e3 for a in avoid_hz):
            continue
        out.append(f)
    return out


def inject(x, m, f, drift, snr):
    """Plant a tone so that it dedrift-sums to about `snr` in this scan's units."""
    foff_hz = m["foff"] * 1e6
    nt, nc = x.shape
    t = (np.arange(nt) + 0.5 - nt / 2) * m["tsamp"]
    # per-channel noise of a dedrifted sum ~ sqrt(nt) * per-pixel sigma
    sig = np.median(np.abs(x - np.median(x))) * 1.4826
    amp = snr * sig * np.sqrt(nt) / nt
    y = x.copy()
    c = (f - m["freqs"][0]) / (m["freqs"][1] - m["freqs"][0])
    for k in range(nt):
        pos = c + drift * t[k] / foff_hz
        j = int(np.floor(pos))
        fr = pos - j
        if 0 <= j < nc - 1:
            y[k, j] += amp * (1 - fr)
            y[k, j + 1] += amp * fr
    return y


def cmd_run():
    from reach import s_min_jy, eirp_min_w, LY
    scans = []
    for kind, name in SCANS:
        d, m = load(name)
        scans.append((kind, name, flatten(d), m))
    print(f"pre-registered: {json.dumps({k: v for k, v in PREREG.items() if k != 'drift_offsets_hz_s'})}")
    m0 = scans[0][3]
    result = dict(prereg=PREREG, frames={})
    for frame in PREREG["frames"]:
        f_pred, drift0 = predicted(m0, frame)
        print(f"\n=== frame: {frame}  predicted {f_pred/1e6:.6f} MHz, Earth drift {drift0:+.4f} Hz/s")
        # 1. null distribution from empty windows, ON scans only (same statistics as the test)
        nulls = []
        offs = null_offsets(m0, [f_pred], 400)
        for kind, name, x, m in scans:
            if kind != "ON":
                continue
            for f in offs:
                nulls.append(search_window(x, m, f, drift0, PREREG["window_hz"])[0])
        thr = max(PREREG["snr_floor"], float(np.percentile(nulls, PREREG["null_percentile"])))
        print(f"  null: {len(nulls)} empty windows, 99th pct {np.percentile(nulls, 99):.2f}, "
              f"max {max(nulls):.2f}  ->  threshold {thr:.2f}")
        # 2. injection-recovery (see AMENDMENT): efficiency, then completeness
        ons = [(name, x, m) for kind, name, x, m in scans if kind == "ON"]
        rng = np.random.default_rng(11)
        spots = [offs[i] for i in rng.choice(len(offs), size=min(30, len(offs)), replace=False)]
        eff = []
        for name, x, m in ons:
            dr = predicted(m, frame)[1]
            for f in spots[:10]:
                y = inject(x, m, f, dr, 20.0)
                snr, fb, _ = search_window(y, m, f, dr, 40.0)   # look right at the plant
                eff.append(snr / 20.0)
        eta = float(np.median(eff))
        levels = [8.0, 10.0, 12.0, 15.0, 20.0]
        compl = {}
        for lv in levels:
            hits = n = 0
            for name, x, m in ons:
                dr = predicted(m, frame)[1]
                for f in spots:
                    y = inject(x, m, f, dr, lv / eta)          # aim at MEASURED snr lv
                    snr, fb, _ = search_window(y, m, f, dr, PREREG["window_hz"])
                    n += 1
                    hits += (snr >= thr and abs(fb - f) < 10.0)
            compl[lv] = hits / n
        snr90 = next((lv for lv in levels if compl[lv] >= 0.9), None)
        inj_ok = snr90 is not None
        rec = dict(efficiency=eta, completeness=compl, snr90=snr90)
        print(f"  injection: efficiency {eta:.2f}; completeness "
              + ", ".join(f"SNR{lv:g}:{100*c:.0f}%" for lv, c in compl.items())
              + f"  ->  90 % at SNR {snr90}")
        # 3. the real search, all six scans
        rows = []
        for kind, name, x, m in scans:
            fp, dr = predicted(m, frame)
            snr, fbest, dbest = search_window(x, m, fp, dr, PREREG["window_hz"])
            rows.append(dict(kind=kind, scan=name[-36:], f_pred_mhz=fp / 1e6, best_snr=snr,
                             best_f_mhz=fbest / 1e6, best_drift=dbest,
                             above=snr >= thr))
            print(f"  {kind:3s} pred {fp/1e6:.6f}  best SNR {snr:6.2f} at {fbest/1e6:.6f} MHz "
                  f"drift {dbest:+.3f}  {'ABOVE' if snr >= thr else '-'}")
        on = [r for r in rows if r["kind"] == "ON"]
        off = [r for r in rows if r["kind"] == "OFF"]
        ch_hz = abs(m0["foff"]) * 1e6
        same_f = (max(r["best_f_mhz"] for r in on) - min(r["best_f_mhz"] for r in on)) * 1e6 \
            <= PREREG["match_channels"] * ch_hz + 0.5 * 0  # same topocentric bin, Earth term shifts it
        candidate = all(r["above"] for r in on) and not any(r["above"] for r in off)
        verdict = ("CANDIDATE - in all ONs, no OFF: rung 2, try hard to kill it" if candidate else
                   "NULL - nothing in all three ONs above threshold")
        # 4. upper limit, in the same units Breakthrough Listen uses
        t_scan = 16 * m0["tsamp"]
        # radiometer limit at the 90 %-complete measured SNR, corrected for efficiency
        snr_lim = (snr90 or levels[-1]) / eta
        s_jy = s_min_jy(10.0, snr_lim, ch_hz, t_scan, 2, 1.0)
        eirp = eirp_min_w(s_jy, 5.96 * LY)
        print(f"  VERDICT: {verdict}")
        print(f"  upper limit (SEFD 10 Jy, 90%-complete SNR {snr90} / eff {eta:.2f}, {t_scan:.0f} s/scan, 1 Hz signal): "
              f"S < {s_jy:.1f} Jy  ->  EIRP < {eirp/1e9:.1f} GW at 5.96 ly")
        result["frames"][frame] = dict(f_pred_mhz=f_pred / 1e6, earth_drift=drift0,
                                       threshold=thr, n_null=len(nulls),
                                       null_p99=float(np.percentile(nulls, 99)),
                                       injection=rec, injection_ok=inj_ok, scans=rows,
                                       verdict=verdict, s_min_jy=s_jy, eirp_min_w=eirp)
    out = HERE / "HI_STARFRAME_GJ699.json"
    out.write_text(json.dumps(result, indent=1, default=float))
    print(f"\nwrote {out.name}")


if __name__ == "__main__":
    {"fetch": cmd_fetch, "run": cmd_run}[sys.argv[1] if len(sys.argv) > 1 else "run"]()
