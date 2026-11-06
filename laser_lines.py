#!/usr/bin/env python3
"""laser_lines.py - optical SETI on spectra the Sloan survey already took.

A laser aimed at us would show up in a star's spectrum as a needle: a bright
line exactly as narrow as the spectrograph can make it (unresolved), at a
wavelength where no element, molecule or sky glow belongs (Townes 1983). SDSS
keeps every separate exposure inside each star's "full" spec file - often
from several NIGHTS - plus the sky spectrum subtracted from each. So a real
line can be demanded to be:

  1. unresolved          matched filter at the instrument's own line width
  2. in EVERY exposure   cosmic rays hit one exposure; a laser hits all
  3. steady in flux      chi^2 across exposures consistent with one value
  4. not the sky         the subtracted sky has no feature there
  5. not a star's line   not at a known stellar emission/absorption line, and
                         not in the median spectrum of OTHER stars of the same
                         type (M dwarfs are full of sharp molecular bands), and
                         not sitting in an absorption line of those stars (a
                         shallower-than-usual line leaves a laser-like bump)
  6. not the instrument  the same observed wavelength must not recur in many
                         OTHER stars of the sample (that is sky/detector)
Survivors are "not yet explained" (rung 1), never "lasers".

PRE-REGISTERED: coadd matched-filter SNR >= 7 to enter; each exposure that
covers the wavelength must show SNR >= 1.5 at the same pixel (+-1) in >= 80 %
of exposures; flux chi^2/dof < 4; fitted width <= 1.6x instrument width; the
observed wavelength (+-2 px) must recur in < 3 % of sample stars. Limits come
from injection-recovery, not from the noise formula.

    python laser_lines.py fetch --n 200     # SDSS DR17 M/K dwarfs, ~2 MB each
    python laser_lines.py run               # search + injection + summary
"""
import argparse
import json
import subprocess
import sys
import urllib.parse
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
DIR = HERE / "data" / "sdss"
SQL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"

PREREG = dict(z_enter=7.0, z_exp=1.5, exp_frac=0.8, chi2_max=4.0, width_max=1.6,
              recur_frac=0.03, recur_px=2, cont_win=61)

# vacuum wavelengths (SDSS uses vacuum), Angstrom: stellar + strongest sky lines
KNOWN = [3727.1, 3934.8, 3969.6, 3971.2, 4102.9, 4341.7, 4862.7, 4960.3, 5008.2,
         5578.9, 5877.3, 5891.6, 5897.6, 6302.0, 6365.5, 6549.9, 6564.6, 6585.3,
         6718.3, 6732.7, 7667.0, 7701.1, 8500.4, 8544.4, 8664.5]
KNOWN_KMS = 150.0       # within this velocity of a known line = vetoed


def fetch(n):
    DIR.mkdir(parents=True, exist_ok=True)
    q = (f"SELECT TOP {n} s.plate, s.mjd, s.fiberid, s.subclass, s.snMedian "
         "FROM specObj s JOIN plateX p ON s.plateID=p.plateID "
         "WHERE s.class='STAR' AND s.run2d='26' AND p.nExp>=12 AND s.snMedian>15 "
         "AND (s.subclass LIKE 'M%' OR s.subclass LIKE 'K%') ORDER BY s.specObjID")
    url = SQL + "?" + urllib.parse.urlencode(dict(cmd=q, format="json"))
    out = DIR / "sample.json"
    subprocess.run(["curl", "-s", "--max-time", "120", "-o", str(out), url], check=True)
    rows = json.loads(out.read_text())[0]["Rows"]
    got = 0
    for r in rows:
        fn = DIR / f"spec-{r['plate']:04d}-{r['mjd']}-{r['fiberid']:04d}.fits"
        if not fn.exists():
            u = (f"https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/full/"
                 f"{r['plate']:04d}/{fn.name}")
            subprocess.run(["curl", "-s", "--max-time", "120", "-o", str(fn), u])
        got += fn.exists() and fn.stat().st_size > 1e5
    print(f"{got}/{len(rows)} spectra in {DIR}")


def smooth_cont(flux, ivar, win):
    from scipy.ndimage import median_filter
    good = ivar > 0
    f = np.where(good, flux, np.nan)
    return median_filter(np.nan_to_num(f, nan=np.nanmedian(f[good]) if good.any() else 0.0),
                         size=win, mode="nearest")


def matched(flux, ivar, sig_px, cont_win, model=None):
    """Matched-filter SNR at every pixel after subtracting a model of the star.
    model=None: a running median (fine for smooth stars). For M/K dwarfs the
    model is (broad continuum) x (median of OTHER stars of the same type), so
    the star's own molecular bands cancel and only what is unique remains."""
    good = ivar > 0
    cont = model if model is not None else smooth_cont(flux, ivar, cont_win)
    r = np.where(good, flux - cont, 0.0)
    w = np.where(good, ivar, 0.0)
    h = int(np.ceil(4 * sig_px))
    k = np.exp(-0.5 * (np.arange(-h, h + 1) / sig_px) ** 2)
    num = np.convolve(r * w, k, "same")
    den = np.sqrt(np.convolve(w, k ** 2, "same"))
    amp = num / np.where(den > 0, den ** 2, np.inf)     # best-fit line amplitude
    z = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)
    return z, amp, r


def read(fn):
    from astropy.io import fits
    h = fits.open(fn, memmap=False)
    co = h[1].data
    exps = []
    for hdu in h[4:]:
        d = hdu.data
        o = np.argsort(d["loglam"])      # blue-camera frames are stored in DESCENDING
        exps.append(dict(name=hdu.name,  # wavelength order; everything downstream needs ascending
                         lam=10 ** d["loglam"][o], flux=d["flux"][o].astype(float),
                         ivar=d["ivar"][o].astype(float), sky=d["sky"][o].astype(float),
                         wdisp=float(np.median(d["wdisp"]))))
    co = dict(lam=10 ** co["loglam"], flux=co["flux"].astype(float),
              ivar=co["ivar"].astype(float), wdisp=float(np.median(co["wdisp"])))
    return co, exps


def near_known(lam):
    return any(abs(lam - k) / k * 3e5 < KNOWN_KMS for k in KNOWN)


def model_for(lam, flux, ivar, template):
    if template is None:
        return None
    tl, tt = template
    ok = (ivar > 0) & (lam > 1000)
    T = np.interp(lam, tl, tt, left=np.nan, right=np.nan)
    # Each star moves at its own radial velocity; a 1-pixel (69 km/s) misalignment
    # leaves derivative-shaped residuals at every metal line. Find the shift.
    best, best_s = np.inf, 0.0
    n0 = flux / np.where(smooth_cont(flux, ivar, 301) > 0, smooth_cont(flux, ivar, 301), np.nan)
    for sh in np.arange(-4.0, 4.01, 0.25):             # pixels of 1e-4 dex
        Ts = np.interp(np.log10(np.where(lam > 1000, lam, 1.0)) - sh * 1e-4,
                       np.log10(tl), tt, left=np.nan, right=np.nan)
        m = ok & np.isfinite(Ts) & np.isfinite(n0)
        if m.sum() < 500:
            continue
        r = np.nanmedian(np.abs(n0[m] - Ts[m]))
        if r < best:
            best, best_s = r, sh
    T = np.interp(np.log10(np.where(lam > 1000, lam, 1.0)) - best_s * 1e-4,
                  np.log10(tl), tt, left=np.nan, right=np.nan)
    broad = smooth_cont(flux / np.where(np.isfinite(T) & (T > 0), T, 1.0), ivar, 301)
    m = broad * T
    return np.where(np.isfinite(m), m, smooth_cont(flux, ivar, PREREG["cont_win"]))


def examine(co, exps, inject=None, template=None):
    """Return candidate dicts for one star. inject=(lam, amp) plants an unresolved
    line of that peak amplitude (flux units) into every exposure AND the coadd."""
    def plant(lam, flux, sig_px):
        if inject is None:
            return flux
        l0, a0 = inject
        v = np.where(lam > 1000)[0]                    # exposure files pad with empty rows
        pix = np.interp(l0, lam[v], v.astype(float))
        return flux + a0 * np.exp(-0.5 * ((np.arange(len(lam)) - pix) / sig_px) ** 2)
    sig = max(co["wdisp"], 0.8)
    fco = plant(co["lam"], co["flux"], sig)
    mco = model_for(co["lam"], co["flux"], co["ivar"], template)   # model from the UNPLANTED star
    z, amp, r = matched(fco, co["ivar"], sig, PREREG["cont_win"], mco)
    cands = []
    edge = PREREG["cont_win"]
    peaks = np.where((z[1:-1] > z[:-2]) & (z[1:-1] >= z[2:]) & (z[1:-1] >= PREREG["z_enter"]))[0] + 1
    for p in peaks:
        if p < edge or p > len(z) - edge:
            continue
        lam = co["lam"][p]
        c = dict(lam=float(lam), z=float(z[p]), amp=float(amp[p]), vetoes=[])
        # Does the star model fit here at all? Residual chi^2 in +-40 px, peak masked.
        w = np.r_[max(0, p - 40):max(0, p - 4), p + 5:min(len(r), p + 41)]
        iv = co["ivar"][w]
        chi = float(np.sum(r[w] ** 2 * iv) / max(1, (iv > 0).sum()))
        c["model_chi2"] = chi
        if chi > 3.0:
            c["vetoes"].append(f"star model does not fit here (chi2 {chi:.1f})")
        if near_known(lam):
            c["vetoes"].append("known line")
        # A star whose absorption line is SHALLOWER than its comparison stars'
        # leaves a narrow positive residual exactly at that line - the commonest
        # laser impostor once same-type templates are used.
        if mco is not None:
            sm = smooth_cont(mco, co["ivar"], 61)
            q0, q1 = max(0, p - 2), min(len(mco), p + 3)
            depth = float(np.min(mco[q0:q1] / np.where(sm[q0:q1] > 0, sm[q0:q1], np.inf)))
            c["model_depth"] = depth
            if depth < 0.97:
                c["vetoes"].append(f"at a stellar absorption line (comparison stars dip to {depth:.2f})")
        # width: compare matched SNR at 1x vs 2x instrument width
        z2, _, _ = matched(fco, co["ivar"], 2 * sig, PREREG["cont_win"], mco)
        if z2[p] > z[p] * 1.15:
            c["vetoes"].append("resolved (wider than the instrument)")
        # per-exposure persistence + flux consistency + sky
        zs, amps, errs, skyz = [], [], [], []
        for e in exps:
            v = e["lam"][(e["ivar"] > 0) & (e["lam"] > 1000)]   # files pad with empty rows
            if v.size == 0 or not (v.min() + 5 < lam < v.max() - 5):
                continue
            se = max(e["wdisp"], 0.8)
            ok_e = (e["ivar"] > 0) & (e["lam"] > 1000)
            me = model_for(e["lam"], e["flux"], e["ivar"], template) if ok_e.any() else None
            ze, ae, _ = matched(plant(e["lam"], e["flux"], se), e["ivar"], se, PREREG["cont_win"], me)
            q = int(np.argmin(abs(e["lam"] - lam)))
            j = q - 1 + int(np.argmax(ze[q - 1:q + 2]))
            zs.append(ze[j])
            amps.append(ae[j])
            errs.append(ae[j] / ze[j] if ze[j] != 0 else np.inf)
            # sky: is there a sky feature here? (sky spectrum, its own noise ~ sqrt)
            s = e["sky"]
            loc = s[max(0, q - 30):q + 31]
            skyz.append((s[q] - np.median(loc)) / (np.median(np.abs(loc - np.median(loc))) * 1.4826 + 1e-9))
        if zs:
            frac = float(np.mean(np.array(zs) >= PREREG["z_exp"]))
            c["exp_frac"], c["n_exp"] = frac, len(zs)
            if frac < PREREG["exp_frac"]:
                c["vetoes"].append(f"not in every exposure ({frac:.0%} of {len(zs)})")
            a, e_ = np.array(amps), np.array(errs)
            ok = np.isfinite(e_) & (e_ > 0)
            if ok.sum() >= 2:
                wm = np.sum(a[ok] / e_[ok] ** 2) / np.sum(1 / e_[ok] ** 2)
                chi2 = float(np.sum(((a[ok] - wm) / e_[ok]) ** 2) / (ok.sum() - 1))
                c["chi2"] = chi2
                if chi2 > PREREG["chi2_max"]:
                    c["vetoes"].append(f"flux not steady (chi2/dof {chi2:.1f})")
            if np.median(skyz) > 5:
                c["vetoes"].append("sky line here")
        else:
            c["vetoes"].append("no exposures cover it")
        cands.append(c)
    return cands


def build_templates(files):
    """Leave-one-out median spectrum of the OTHER stars of the same spectral
    type (e.g. all M3 but this one), normalised by a broad continuum, on the
    common SDSS log-lambda grid (1e-4 dex steps)."""
    meta = {f"spec-{r['plate']:04d}-{r['mjd']}-{r['fiberid']:04d}": r["subclass"][:2]
            for r in json.loads((DIR / "sample.json").read_text())[0]["Rows"]}
    grid = np.arange(35500, 40200) * 1e-4
    norm = {}
    for fn in files:
        co, _ = read(fn)
        good = co["ivar"] > 0
        n = co["flux"] / np.where(smooth_cont(co["flux"], co["ivar"], 301) > 0,
                                  smooth_cont(co["flux"], co["ivar"], 301), np.nan)
        n[~good] = np.nan
        norm[fn.stem] = np.interp(grid, np.log10(co["lam"]), n, left=np.nan, right=np.nan)
    out = {}
    for stem, typ in meta.items():
        if stem not in norm:
            continue
        def num(t):
            return ("KM".index(t[0]) * 10 + int(t[1])) if t[0] in "KM" and t[1:].isdigit() else None
        k = num(typ)
        for width in (0, 1, 2):
            mates = [norm[o] for o, t in meta.items() if o != stem and o in norm
                     and k is not None and num(t) is not None and abs(num(t) - k) <= width]
            if len(mates) >= 3:
                out[stem] = (10 ** grid, np.nanmedian(np.array(mates), axis=0))
                break
    print(f"templates: {len(out)}/{len(files)} stars have >= 3 same-type mates")
    return out


def run(args):
    files = sorted(DIR.glob("spec-*.fits"))
    if not files:
        sys.exit("no spectra: run `python laser_lines.py fetch` first")
    print(f"pre-registered: {PREREG}")
    templates = build_templates(files)
    allc, skipped = [], []
    for fn in files:
        try:
            co, exps = read(fn)
        except Exception as ex:
            print(f"  skip {fn.name}: {ex}")
            continue
        if fn.stem not in templates:
            skipped.append(fn.stem)
            continue
        for c in examine(co, exps, template=templates.get(fn.stem)):
            c["star"] = fn.stem
            allc.append(c)
    # recurrence across stars: the same observed wavelength in many stars = not a star
    nstar = len(files)
    lams = np.array([c["lam"] for c in allc])
    for c in allc:
        dl = c["lam"] * 69e-6 * PREREG["recur_px"]      # SDSS pixel ~69 km/s
        stars = {d["star"] for d in allc if abs(d["lam"] - c["lam"]) < dl}
        c["recur"] = len(stars)
        if len(stars) / nstar >= PREREG["recur_frac"]:
            c["vetoes"].append(f"recurs in {len(stars)} stars (sky/instrument)")
    survivors = [c for c in allc if not c["vetoes"]]
    from collections import Counter
    vc = Counter(v.split(" (")[0] for c in allc for v in c["vetoes"])
    print(f"\nskipped (no same-type comparison stars): {skipped}")
    print(f"{nstar - len(skipped)} stars searched, {len(allc)} peaks with coadd SNR >= {PREREG['z_enter']}")
    for v, n in vc.most_common():
        print(f"  vetoed: {v:45s} {n}")
    print(f"  SURVIVORS (rung 1, 'not yet explained'): {len(survivors)}")
    for c in survivors[:20]:
        print(f"    {c['star']}  {c['lam']:.2f} A  SNR {c['z']:.1f}  in {c.get('exp_frac',0):.0%} of "
              f"{c.get('n_exp',0)} exposures  chi2 {c.get('chi2',float('nan')):.1f}")
    # injection-recovery: plant unresolved lines in 40 stars at clean wavelengths
    rng = np.random.default_rng(3)
    pick = [files[i] for i in rng.choice(len(files), size=min(40, len(files)), replace=False)]
    levels = [5.0, 8.0, 12.0, 15.0, 20.0]  # planted coadd SNR (by the coadd's own noise)
    comp = {}
    flux_at = {}
    searchable = []
    for lv in levels:
        ok = n = tried = 0
        fl = []
        for fn in pick:
            if fn.stem not in templates:
                continue
            co, exps = read(fn)
            for _ in range(3):
                lam = float(rng.uniform(co["lam"][200], co["lam"][-200]))
                if near_known(lam):
                    continue
                # amplitude giving matched SNR ~ lv on the coadd
                sig = max(co["wdisp"], 0.8)
                q = int(np.argmin(abs(co["lam"] - lam)))
                iv = co["ivar"][max(0, q - 5):q + 6]
                k = np.exp(-0.5 * (np.arange(-5, 6) / sig) ** 2)
                a0 = lv / np.sqrt(np.sum(k ** 2 * iv)) if np.sum(iv) > 0 else 0
                if a0 <= 0:
                    continue
                found = [c for c in examine(co, exps, inject=(lam, a0), template=templates.get(fn.stem))
                         if abs(c["lam"] - lam) < 3 * lam * 69e-6
                         and not [v for v in c["vetoes"] if not v.startswith("recurs")]]
                # Is this spot searchable at all, before anything is planted?
                mco = model_for(co["lam"], co["flux"], co["ivar"], templates[fn.stem])
                sm = smooth_cont(mco, co["ivar"], 61)
                q0 = int(np.argmin(abs(co["lam"] - lam)))
                dep = float(np.min(mco[q0 - 2:q0 + 3] / np.where(sm[q0 - 2:q0 + 3] > 0, sm[q0 - 2:q0 + 3], np.inf)))
                _, _, r0 = matched(co["flux"], co["ivar"], sig, PREREG["cont_win"], mco)
                w = np.r_[max(0, q0 - 40):q0 - 4, q0 + 5:q0 + 41]
                chi0 = float(np.sum(r0[w] ** 2 * co["ivar"][w]) / max(1, (co["ivar"][w] > 0).sum()))
                tried += 1
                if dep < 0.97 or chi0 > 3.0:
                    continue        # in an absorption line, or the model fails: not searchable
                n += 1
                ok += bool(found)
                fl.append(a0 * sig * np.sqrt(2 * np.pi) * lam * 69e-6 / 1.0)   # line flux, 1e-17 erg/s/cm^2
        comp[lv] = ok / max(n, 1)
        searchable.append(n / max(tried, 1))
        flux_at[lv] = float(np.median(fl)) if fl else float("nan")
        print(f"  injection SNR {lv:4.0f}: recovered {ok}/{n} = {comp[lv]:.0%}   "
              f"median line flux {flux_at[lv]:.1f}e-17 erg/s/cm^2")
    sfrac = float(np.mean(searchable))
    print(f"  searchable fraction of the spectrum (not in a stellar absorption line,"
          f" star model fits): {sfrac:.0%}")
    # flagged AFTER the run, not pre-registered: Earth's own water-vapour and O2
    # bands, where telluric-correction residuals can leave narrow bumps
    TELLURIC = [(6865, 6935), (7160, 7340), (7590, 7700), (8100, 8350), (8900, 9200)]
    for c in survivors:
        c["telluric_band"] = any(a <= c["lam"] <= b for a, b in TELLURIC)
    nt = sum(c["telluric_band"] for c in survivors)
    print(f"  post-hoc: {nt}/{len(survivors)} survivors sit in a telluric (Earth-atmosphere) band")
    lim = next((lv for lv in levels if comp[lv] >= 0.9), None)
    print(f"  90 %-complete at coadd SNR {lim} -> line flux ~{flux_at.get(lim, float('nan')):.0f}e-17 erg/s/cm^2")
    (HERE / "LASER_LINES_SDSS.json").write_text(json.dumps(dict(
        prereg=PREREG, n_stars=nstar, searchable_fraction=sfrac, n_peaks=len(allc), veto_counts=dict(vc),
        survivors=survivors, completeness=comp, flux_1e17=flux_at, snr90=lim), indent=1, default=float))
    print("wrote LASER_LINES_SDSS.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--n", type=int, default=200)
    sub.add_parser("run")
    a = ap.parse_args()
    fetch(a.n) if a.cmd == "fetch" else run(a)
