#!/usr/bin/env python3
"""seti_io.py - one loader for everything seti-cookbook looks at.

Every detector, recipe, waterfall and MCP tool in this repo speaks ONE data
type: a `Spectrogram` - a (ntime, nchan) block of power with frequency
ASCENDING and physical axes attached. That single convention is what lets a
stranger write a recipe (see RECIPES.md) without knowing whether the bytes came
from Green Bank, Parkes, an RTL-SDR capture or a synthetic test signal.

Loaders:
  load_bl(path, ...)     Breakthrough Listen / SIGPROC: .h5 (bitshuffle) or .fil
  load_iq(path, fs, ...) raw interleaved int16 IQ (.cs16) -> STFT spectrogram
  synth(kind, ...)       physically-correct SYNTHETIC signals (drift, FRB,
                         pulsar, comb, spread, noise) - so selftests and
                         teaching figures need no downloads at all
  open_any(path, ...)     dispatch on extension

GPU is OPTIONAL everywhere (see gpu_optional_law): set SETI_COOKBOOK_GPU=1 and the
heavy loops (dedispersion, folding) run on cupy if it is installed; with no GPU,
or no cupy, or the variable unset, the identical numpy path runs. Nothing in
this repo ever REQUIRES a GPU.

  python seti_io.py <file>            # header / summary
  python seti_io.py selftest
"""
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# --- physical constants used across the repo -------------------------------
HI_MHZ = 1420.405751768        # neutral hydrogen 21 cm rest frequency
OH_MHZ = (1612.231, 1665.402, 1667.359, 1720.530)   # the OH maser quartet
C_KMS = 299792.458
# Dispersion delay constant. t = K * DM * nu_MHz^-2 seconds, K in MHz^2 pc^-1 cm^3 s
DM_CONST = 4.148808e3          # => t_ms = 4.148808 * DM * nu_GHz^-2


def xp(force_cpu=False):
    """Return the array module to compute with: cupy if the user asked for GPU
    and it is actually importable, else numpy. Honest fallback, never fatal."""
    if force_cpu or os.environ.get("SETI_COOKBOOK_GPU", os.environ.get("SETITUNA_GPU", "0")) not in ("1", "true", "yes"):
        return np
    try:
        import cupy
        cupy.zeros(1)          # prove a device exists, not just the module
        return cupy
    except Exception:
        return np


def _tonp(a):
    return a.get() if hasattr(a, "get") else np.asarray(a)


# --------------------------------------------------------------------------
@dataclass
class Spectrogram:
    """(ntime, nchan) power, frequency ASCENDING. The repo's only data type.

    data    float32/float64 array, data[t, c] = power in channel c at time t
    f0_mhz  centre frequency of channel 0 (the LOWEST frequency)
    df_mhz  channel width in MHz, always POSITIVE
    t0_s    time of sample 0, seconds from the start of the observation
    dt_s    sample cadence in seconds
    meta    dict: source_name, telescope, tstart_mjd, ra_hr, dec_deg, origin...
    """
    data: np.ndarray
    f0_mhz: float
    df_mhz: float
    t0_s: float = 0.0
    dt_s: float = 1.0
    meta: dict = field(default_factory=dict)

    # --- axes -------------------------------------------------------------
    @property
    def ntime(self):
        return self.data.shape[0]

    @property
    def nchan(self):
        return self.data.shape[1]

    def freqs_mhz(self):
        return self.f0_mhz + self.df_mhz * np.arange(self.nchan)

    def times_s(self):
        return self.t0_s + self.dt_s * np.arange(self.ntime)

    @property
    def f_lo(self):
        return self.f0_mhz

    @property
    def f_hi(self):
        return self.f0_mhz + self.df_mhz * (self.nchan - 1)

    @property
    def duration_s(self):
        return self.dt_s * self.ntime

    @property
    def res_hz(self):
        return self.df_mhz * 1e6

    def chan_of(self, f_mhz):
        return int(round((f_mhz - self.f0_mhz) / self.df_mhz))

    # --- shaping ----------------------------------------------------------
    def crop(self, f_start=None, f_stop=None, t_start=None, t_stop=None):
        c0 = 0 if f_start is None else max(0, self.chan_of(f_start))
        c1 = self.nchan if f_stop is None else min(self.nchan, self.chan_of(f_stop) + 1)
        i0 = 0 if t_start is None else max(0, int(round((t_start - self.t0_s) / self.dt_s)))
        i1 = self.ntime if t_stop is None else min(
            self.ntime, int(round((t_stop - self.t0_s) / self.dt_s)) + 1)
        if c1 <= c0 or i1 <= i0:
            raise ValueError("crop selected an empty region")
        return Spectrogram(self.data[i0:i1, c0:c1], self.f0_mhz + self.df_mhz * c0,
                           self.df_mhz, self.t0_s + self.dt_s * i0, self.dt_s,
                           dict(self.meta, cropped=True))

    def integrated(self):
        """Time-averaged spectrum (the 'total power' spectrum)."""
        return np.asarray(self.data, np.float64).mean(0)

    def timeseries(self):
        """Frequency-averaged power vs time (the light curve)."""
        return np.asarray(self.data, np.float64).mean(1)

    def bandpass_normalized(self, width=None):
        """Divide out the instrument bandpass with a running median so that a
        LINE (HI, a maser, a carrier) stands out from the telescope's own shape.
        Returns a new Spectrogram in units of 'times the local continuum'."""
        from scipy.ndimage import median_filter
        w = width or max(31, min(1001, self.nchan // 16 * 2 + 1))
        sp = self.integrated()
        base = median_filter(sp, w, mode="nearest")
        base = np.where(base <= 0, np.nanmedian(base[base > 0]) or 1.0, base)
        out = np.asarray(self.data, np.float64) / base[None, :]
        return Spectrogram(out, self.f0_mhz, self.df_mhz, self.t0_s, self.dt_s,
                           dict(self.meta, normalized=True))

    def decimate(self, tfac=1, ffac=1):
        """Block-average in time and/or frequency (for plotting big files)."""
        d = np.asarray(self.data, np.float64)
        if ffac > 1:
            n = (self.nchan // ffac) * ffac
            d = d[:, :n].reshape(d.shape[0], n // ffac, ffac).mean(2)
        if tfac > 1:
            n = (d.shape[0] // tfac) * tfac
            d = d[:n].reshape(n // tfac, tfac, d.shape[1]).mean(1)
        return Spectrogram(d, self.f0_mhz + self.df_mhz * (ffac - 1) / 2,
                           self.df_mhz * ffac, self.t0_s + self.dt_s * (tfac - 1) / 2,
                           self.dt_s * tfac, dict(self.meta))

    # --- physics ----------------------------------------------------------
    def dedisperse(self, dm, force_cpu=False):
        """Shift each channel back by the cold-plasma dispersion delay for a
        trial DM (pc cm^-3), referenced to the HIGHEST frequency. An FRB or
        pulsar pulse that was a curved sweep becomes a vertical line.
        delay(nu) = 4.148808e3 * DM * (nu_MHz^-2 - nu_hi_MHz^-2) seconds."""
        m = xp(force_cpu)
        f = self.freqs_mhz()
        delay = DM_CONST * dm * (f ** -2 - self.f_hi ** -2)
        shift = np.round(delay / self.dt_s).astype(int)
        d = m.asarray(self.data)
        out = m.zeros_like(d)
        for c in range(self.nchan):
            s = int(shift[c])
            if s == 0:
                out[:, c] = d[:, c]
            elif s < self.ntime:
                out[:self.ntime - s, c] = d[s:, c]
        return Spectrogram(_tonp(out), self.f0_mhz, self.df_mhz, self.t0_s,
                           self.dt_s, dict(self.meta, dedispersed_dm=dm))

    def fold(self, period_s, nbins=64, force_cpu=False):
        """Fold the light curve at a trial period - how a hobbyist actually
        sees a pulsar: single pulses are buried, the FOLDED profile is not.
        Returns (phase, profile)."""
        m = xp(force_cpu)
        ts = m.asarray(self.timeseries())
        t = m.asarray(self.times_s())
        ph = (t / period_s) % 1.0
        idx = m.asarray((ph * nbins).astype(int)) if hasattr(ph, "astype") else None
        idx = (ph * nbins).astype(int)
        prof = m.zeros(nbins)
        cnt = m.zeros(nbins)
        for b in range(nbins):
            sel = idx == b
            n = sel.sum()
            if n:
                prof[b] = ts[sel].sum()
                cnt[b] = n
        prof = _tonp(prof) / np.maximum(_tonp(cnt), 1)
        return np.arange(nbins) / nbins, prof

    def summary(self):
        d = np.asarray(self.data, np.float64)
        return {
            "source_name": self.meta.get("source_name", "?"),
            "telescope": self.meta.get("telescope", "?"),
            "origin": self.meta.get("origin", "?"),
            "ntime": self.ntime, "nchan": self.nchan,
            "f_lo_mhz": round(self.f_lo, 6), "f_hi_mhz": round(self.f_hi, 6),
            "res_hz": round(self.res_hz, 4),
            "dt_s": self.dt_s, "duration_s": round(self.duration_s, 3),
            "tstart_mjd": self.meta.get("tstart_mjd"),
            "power_median": float(np.median(d)), "power_max": float(d.max()),
        }

    def __repr__(self):
        return (f"<Spectrogram {self.ntime}x{self.nchan} "
                f"{self.f_lo:.4f}-{self.f_hi:.4f} MHz @ {self.res_hz:.3f} Hz, "
                f"{self.duration_s:.1f}s src={self.meta.get('source_name','?')}>")


# --------------------------------------------------------------------------
TELESCOPE_IDS = {0: "Fake", 1: "Arecibo", 2: "Ooty", 3: "Nancay", 4: "Parkes",
                 5: "Jodrell", 6: "GBT", 8: "Effelsberg", 10: "SRT",
                 64: "MeerKAT", 65: "KAT-7"}


def _bl_attrs(a):
    tid = int(a.get("telescope_id", -1)) if "telescope_id" in a else -1
    return {
        "source_name": (a.get("source_name", b"?").decode()
                        if isinstance(a.get("source_name"), bytes)
                        else str(a.get("source_name", "?"))),
        "telescope": TELESCOPE_IDS.get(tid, f"id{tid}"),
        "telescope_id": tid,
        "tstart_mjd": float(a["tstart"]) if "tstart" in a else None,
        "ra_hr": float(a["src_raj"]) if "src_raj" in a else None,
        "dec_deg": float(a["src_dej"]) if "src_dej" in a else None,
    }


def load_bl(path, f_start=None, f_stop=None, t_start=None, t_stop=None,
            max_chans=1 << 20, max_time=4096):
    """Load a Breakthrough Listen / SIGPROC data product.

    .h5  read directly with h5py (+hdf5plugin for BL's bitshuffle filter) so a
         window can be pulled out of a multi-GB file without loading it all.
    .fil read with blimpy.

    f_start/f_stop in MHz, t_start/t_stop in seconds. Frequency is flipped to
    ASCENDING on the way out - BL files usually store it descending."""
    path = Path(path)
    if path.suffix.lower() in (".fil", ".filterbank"):
        return _load_fil(path, f_start, f_stop, t_start, t_stop)
    try:
        import hdf5plugin      # noqa: F401  (registers BL's bitshuffle filter)
    except ImportError:
        pass
    try:
        import h5py
    except ImportError:
        raise SystemExit("reading BL .h5 needs h5py + hdf5plugin: "
                         "pip install h5py hdf5plugin   (or blimpy)")
    with h5py.File(path, "r") as h:
        d = h["data"]
        a = dict(d.attrs)
        fch1 = float(a["fch1"])
        foff = float(a["foff"])
        tsamp = float(a["tsamp"])
        nch = d.shape[-1]
        nt = d.shape[0]
        # channel index range for the requested MHz window (foff may be < 0)
        def idx(fm):
            return int(round((fm - fch1) / foff))
        if f_start is None and f_stop is None:
            c0, c1 = 0, nch
        else:
            lo = idx(f_stop if foff < 0 else f_start) if (f_stop if foff < 0 else f_start) is not None else 0
            hi = idx(f_start if foff < 0 else f_stop) if (f_start if foff < 0 else f_stop) is not None else nch - 1
            c0, c1 = max(0, min(lo, hi)), min(nch, max(lo, hi) + 1)
        i0 = 0 if t_start is None else max(0, int(t_start / tsamp))
        i1 = nt if t_stop is None else min(nt, int(np.ceil(t_stop / tsamp)))
        if c1 - c0 > max_chans:
            raise ValueError(f"{c1-c0} channels requested (> max_chans={max_chans}); "
                             "narrow f_start/f_stop or raise max_chans")
        i1 = min(i1, i0 + max_time)
        blk = np.array(d[i0:i1, 0, c0:c1], dtype=np.float32)
        meta = _bl_attrs(a)
    f_c0 = fch1 + foff * c0
    if foff < 0:                       # flip to ascending
        blk = blk[:, ::-1]
        f0 = f_c0 + foff * (blk.shape[1] - 1)
    else:
        f0 = f_c0
    meta.update(origin=str(path.name), path=str(path), fmt="BL-HDF5",
                chan_offset=c0, time_offset=i0)
    return Spectrogram(blk, f0, abs(foff), i0 * tsamp, tsamp, meta)


def _load_fil(path, f_start, f_stop, t_start, t_stop):
    try:
        from blimpy import Waterfall
    except ImportError:
        raise SystemExit("reading .fil needs blimpy: pip install blimpy")
    w = Waterfall(str(path), f_start=f_start, f_stop=f_stop,
                  t_start=None if t_start is None else int(t_start),
                  t_stop=None if t_stop is None else int(t_stop))
    d = np.asarray(w.data[:, 0, :], np.float32)
    foff = float(w.header["foff"])
    fch1 = float(w.header["fch1"])
    meta = _bl_attrs({k: v for k, v in w.header.items()})
    if foff < 0:
        d = d[:, ::-1]
        f0 = fch1 + foff * (d.shape[1] - 1)
    else:
        f0 = fch1
    meta.update(origin=path.name, path=str(path), fmt="SIGPROC-fil")
    return Spectrogram(d, f0, abs(foff), 0.0, float(w.header["tsamp"]), meta)


def load_iq(path, fs, nfft=4096, nmax=8_000_000, centre_mhz=0.0, overlap=0):
    """STFT a raw interleaved int16 (.cs16) capture into a Spectrogram, so an
    RTL-SDR/SDRplay capture from your own rig goes through the same recipes as
    Green Bank data. centre_mhz labels the axis with the tuner's centre."""
    raw = np.fromfile(str(path), np.int16, count=2 * nmax).astype(np.float32) / 32768.0
    x = raw[0::2] + 1j * raw[1::2]
    step = nfft - overlap
    nrow = max(1, (len(x) - nfft) // step + 1)
    win = np.hanning(nfft).astype(np.float32)
    out = np.empty((nrow, nfft), np.float32)
    for i in range(nrow):
        seg = x[i * step:i * step + nfft] * win
        out[i] = np.fft.fftshift(np.abs(np.fft.fft(seg)) ** 2)
    df_mhz = (fs / nfft) / 1e6
    f0 = centre_mhz - (fs / 2) / 1e6
    return Spectrogram(out, f0, df_mhz, 0.0, step / fs,
                       dict(source_name=Path(path).stem, telescope="local SDR",
                            origin=Path(path).name, fmt="cs16-IQ", fs_hz=fs))


# --------------------------------------------------------------------------
def synth(kind="drift", ntime=64, nchan=1024, f0_mhz=1400.0, df_mhz=2.79e-6,
          dt_s=1.0, seed=0, **kw):
    """Physically-correct synthetic spectrograms. Every teaching figure and
    recipe selftest in this repo can run with NO download - and every signature
    below is generated from the same formula quoted in SETI_HISTORY.md.

    kind:
      noise     pure radiometer noise (the null hypothesis)
      drift     narrowband tone with a Doppler DRIFT RATE (Hz/s) - what
                turboSETI hunts; kw: snr, drift_hz_s, f_mhz
      zerodrift the same tone with drift exactly 0 - the classic RFI tell
      frb       dispersed broadband ms pulse; kw: dm, width_s, snr
      pulsar    periodic dispersed pulse train; kw: period_s, dm, snr
      comb      uniformly Hz-spaced tones (our comb detector's target)
      spread    wideband noise-like digital signal (cyclostationary, spectrum-flat)
      maser     bright narrow natural line, no drift, ~few kHz wide
      hi        broad ~100 kHz galactic hydrogen line with a velocity offset
    """
    rng = np.random.default_rng(seed)
    d = rng.chisquare(2, size=(ntime, nchan)).astype(np.float32) / 2.0
    f = f0_mhz + df_mhz * np.arange(nchan)
    t = dt_s * np.arange(ntime)
    hz = df_mhz * 1e6

    def add_tone(fmhz, drift, snr, width_ch=1.0):
        for i, ti in enumerate(t):
            fc = fmhz + drift * ti / 1e6
            c = (fc - f0_mhz) / df_mhz
            lo, hi = int(c - 4 * width_ch), int(c + 4 * width_ch) + 1
            for cc in range(max(0, lo), min(nchan, hi)):
                d[i, cc] += snr * np.exp(-0.5 * ((cc - c) / width_ch) ** 2)

    if kind == "noise":
        pass
    elif kind in ("drift", "zerodrift"):
        add_tone(kw.get("f_mhz", f0_mhz + df_mhz * nchan / 2),
                 0.0 if kind == "zerodrift" else kw.get("drift_hz_s", -0.35),
                 kw.get("snr", 20.0), kw.get("width_ch", 1.0))
    elif kind == "maser":
        add_tone(kw.get("f_mhz", f0_mhz + df_mhz * nchan * 0.4), 0.0,
                 kw.get("snr", 25.0), kw.get("width_ch", max(1.0, 3000.0 / hz)))
    elif kind == "hi":
        v = kw.get("v_kms", 20.0)
        fc = HI_MHZ * (1 - v / C_KMS)
        sig_ch = max(1.0, (kw.get("width_kms", 20.0) / C_KMS * HI_MHZ) / df_mhz)
        c = (fc - f0_mhz) / df_mhz
        prof = kw.get("snr", 1.0) * np.exp(-0.5 * ((np.arange(nchan) - c) / sig_ch) ** 2)
        d += prof[None, :].astype(np.float32)
    elif kind == "frb":
        dm = kw.get("dm", 560.0)
        w = kw.get("width_s", 0.003)
        t_arr = kw.get("t_s", dt_s * ntime * 0.35)
        snr = kw.get("snr", 30.0)
        fhi = f[-1]
        for c in range(nchan):
            tt = t_arr + DM_CONST * dm * (f[c] ** -2 - fhi ** -2)
            i = int(round(tt / dt_s))
            if 0 <= i < ntime:
                d[i, c] += snr * (1.0 + 0.3 * rng.standard_normal())
                if w > dt_s and i + 1 < ntime:
                    d[i + 1, c] += snr * 0.5
    elif kind == "pulsar":
        per = kw.get("period_s", 0.714)     # PSR B0329+54
        dm = kw.get("dm", 26.8)
        snr = kw.get("snr", 3.0)
        fhi = f[-1]
        k = 0
        while k * per < dt_s * ntime:
            for c in range(nchan):
                tt = k * per + DM_CONST * dm * (f[c] ** -2 - fhi ** -2)
                i = int(round(tt / dt_s))
                if 0 <= i < ntime:
                    d[i, c] += snr
            k += 1
    elif kind == "comb":
        n = kw.get("teeth", 12)
        spacing_ch = kw.get("spacing_ch", nchan // (n + 2))
        for j in range(n):
            add_tone(f0_mhz + df_mhz * (10 + j * spacing_ch), 0.0, kw.get("snr", 12.0))
    elif kind == "spread":
        c0 = int(nchan * 0.3)
        c1 = int(nchan * 0.7)
        d[:, c0:c1] += kw.get("snr", 1.2) * rng.chisquare(2, (ntime, c1 - c0)) / 2.0
    else:
        raise ValueError(f"unknown synth kind {kind!r}")
    return Spectrogram(d, f0_mhz, df_mhz, 0.0, dt_s,
                       dict(source_name=f"SYNTH-{kind}", telescope="synthetic",
                            origin=f"synth:{kind}", fmt="synthetic", synth_kind=kind,
                            synth_kw=kw, truth=dict(kind=kind, **kw)))


def open_any(path, **kw):
    """Dispatch by extension.

    'synth:<kind>' loads a synthetic signal; parameters may be appended as
    comma-separated key=value pairs, e.g.
        synth:drift,drift_hz_s=-0.35,snr=25,ntime=64
        synth:frb,dm=560,f0_mhz=4000,df_mhz=4,dt_s=0.001,ntime=400,nchan=512
    """
    s = str(path)
    if s.startswith("synth:"):
        body = s.split(":", 1)[1]
        parts = [p for p in body.split(",") if p]
        kind = parts[0]
        for kv in parts[1:]:
            k, _, v = kv.partition("=")
            try:
                kw[k.strip()] = int(v) if v.strip().lstrip("-").isdigit() else float(v)
            except ValueError:
                kw[k.strip()] = v.strip()
        return synth(kind, **kw)
    ext = Path(s).suffix.lower()
    if ext in (".h5", ".hdf5", ".fil", ".filterbank"):
        return load_bl(s, **kw)
    if ext in (".cs16", ".iq", ".bin", ".dat", ".c16"):
        fs = kw.pop("fs", None)
        if fs is None:
            raise ValueError("raw IQ needs fs= (sample rate in Hz)")
        return load_iq(s, fs, **kw)
    raise ValueError(f"don't know how to open {path}")


def list_data(root=None, pattern=("*.h5", "*.fil", "*.cs16")):
    """Everything openable under data/ - the MCP server's file list."""
    root = Path(root or (HERE / "data"))
    out = []
    if not root.exists():
        return out
    for pat in pattern:
        for p in sorted(root.glob(pat)):
            out.append({"path": str(p), "name": p.name,
                        "size_mb": round(p.stat().st_size / 1e6, 1)})
    return out


# --------------------------------------------------------------------------
def selftest():
    print("=" * 70)
    print("seti_io selftest - the shared data contract")
    print("=" * 70)
    ok = []
    s = synth("drift", drift_hz_s=-0.4, snr=30)
    ok.append(("ascending frequency", s.df_mhz > 0 and s.f_hi > s.f_lo))
    ok.append(("axes length", len(s.freqs_mhz()) == s.nchan and len(s.times_s()) == s.ntime))
    span = s.f_hi - s.f_lo
    c = s.crop(f_start=s.f_lo + 0.25 * span, f_stop=s.f_lo + 0.75 * span)
    ok.append(("crop narrows", c.nchan < s.nchan and c.f_lo >= s.f_lo))

    # dedispersion must actually straighten a synthetic FRB
    fr = synth("frb", ntime=200, nchan=256, f0_mhz=4000, df_mhz=4.0,
               dt_s=0.001, dm=560, snr=40)
    raw_peak = fr.timeseries().max() / np.median(fr.timeseries())
    dd = fr.dedisperse(560.0)
    dd_peak = dd.timeseries().max() / np.median(dd.timeseries())
    ok.append((f"dedisperse concentrates ({raw_peak:.1f}x -> {dd_peak:.1f}x)",
               dd_peak > 2.5 * raw_peak))

    # folding must recover a pulsar profile
    ps = synth("pulsar", ntime=1200, nchan=64, dt_s=0.01, period_s=0.714,
               dm=0.0, snr=1.5)
    ph, prof = ps.fold(0.714, nbins=32)
    wrong = ps.fold(0.714 * 1.37, nbins=32)[1]
    contrast = (prof.max() - prof.mean()) / prof.std()
    wrong_c = (wrong.max() - wrong.mean()) / wrong.std()
    ok.append((f"fold at true period beats wrong period ({contrast:.2f} vs {wrong_c:.2f})",
               contrast > wrong_c))

    # bandpass normalisation exposes a line hidden under instrument shape.
    # A 20 km/s HI line at 1420 MHz is ~95 kHz = ~34 channels wide here, so the
    # fair test smooths to the line width before taking the peak (matched filter).
    hi = synth("hi", ntime=128, nchan=2048, f0_mhz=1419.5, df_mhz=0.0028,
               snr=0.8, v_kms=30)
    hi.data *= (1.0 + 3.0 * np.exp(-0.5 * ((np.arange(2048) - 1024) / 400.0) ** 2))[None, :]
    n = hi.bandpass_normalized()
    sm = np.convolve(n.integrated(), np.ones(34) / 34, "same")
    fpk = n.freqs_mhz()[int(np.argmax(sm[20:-20]) + 20)]
    vpk = -C_KMS * (fpk - HI_MHZ) / HI_MHZ
    ok.append((f"bandpass-normalise recovers HI at {vpk:+.0f} km/s (truth +30)",
               abs(vpk - 30) < 12))

    m = xp()
    ok.append((f"array backend = {m.__name__} (GPU optional, CPU always works)", True))
    for name, good in ok:
        print(f"  {'PASS' if good else 'FAIL'}  {name}")
    good = all(g for _, g in ok)
    print("=" * 70)
    print(f"RESULT: {sum(g for _, g in ok)}/{len(ok)} checks. "
          f"{'PASS' if good else 'FAIL'}")
    return 0 if good else 1


def main():
    if len(sys.argv) < 2 or sys.argv[1] == "selftest":
        sys.exit(selftest())
    if sys.argv[1] == "list":
        for r in list_data():
            print(f"  {r['size_mb']:9.1f} MB  {r['name']}")
        return
    kw = {}
    if len(sys.argv) >= 3:
        kw["fs"] = float(sys.argv[2])
    s = open_any(sys.argv[1], **kw)
    print(s)
    for k, v in s.summary().items():
        print(f"  {k:16s} {v}")


if __name__ == "__main__":
    main()
