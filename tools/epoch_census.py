#!/usr/bin/env python3
"""epoch_census.py - which stars were observed on TWO OR MORE dates in the same band?

Cross-epoch matching (xepoch.py) and epoch stacking (stack.py) need the same
star observed in the same receiver band on different days, each time as a
real ON/OFF cadence. This census answers that from the Breakthrough Listen
open-data archive, with the exact-name rules of bl_archive.py:

  * an EPOCH is one observing session = one (project, MJD day, centre
    frequency) holding a valid fine-resolution (.0000) cadence - either OFFs
    named <name>_OFF (bl_archive.cadences) or ABACAD with other stars as the
    OFFs (targets.abacad, from a time-window query). Never a scan index.
  * a star qualifies in a band when >= 2 epochs in that band fall on DIFFERENT
    MJD days.
  * for every ON scan the archive's pointing is compared with the star's
    position PROPAGATED TO THAT DATE (proper motion, parallax, RV); an ON more
    than half a beam (GBT FWHM ~ 12.45'/f_GHz) from the star is flagged
    pointing_ok = false. (GJ65A's L-band ONs are 1.2 deg off - KNOWN_ISSUES.)
  * SIMBAD radial velocities outside +-600 km/s are flagged rv_ok = false
    (GJ866's rvz_radvel is 6,824.7 km/s - a redshift field, not a star's RV).
  * frequency coverage per epoch comes from one header read (HTTP range read,
    a few hundred kB) per distinct (centre frequency, node set) configuration.

Only light metadata queries: one file-list query per star and one time-window
query per session, at most 2 requests per second, every raw reply cached under
data/ (gitignored), so a re-run is offline.

    python tools/epoch_census.py                 # stars targets.json already checked (+ GJ699 always)
    python tools/epoch_census.py --unchecked     # also the stars targets.json never file-checked
    python tools/epoch_census.py --offline       # from the cache only
    python tools/epoch_census.py selftest        # offline rules check

Writes EPOCH_CENSUS.csv (one row per star and band, ranked by distance then
total GB) and EPOCH_CENSUS.json (the sessions behind each row: MJDs, files,
sizes, pointing offsets, coverage). Archive metadata only: no signal search,
no E-level (EVIDENCE.md).
"""
import argparse
import csv
import json
import math
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

OUT_CSV = HERE / "EPOCH_CENSUS.csv"
OUT_JSON = HERE / "EPOCH_CENSUS.json"
RAW = HERE / "data" / "epoch_census_raw"
GBT_FWHM_ARCMIN_GHZ = 12.45      # GBT beam FWHM ~ 747"/f_GHz
RV_SANE_KMS = 600.0
COARSE_BYTES = 16 * 2 ** 20 * 4  # one 2^20-channel coarse channel of a 16-spectrum fine-res scan, float32
MIN_INTERVAL_S = 0.5             # <= 2 requests per second to BL's servers

_last = [0.0]


def polite():
    dt = time.time() - _last[0]
    if dt < MIN_INTERVAL_S:
        time.sleep(MIN_INTERVAL_S - dt)
    _last[0] = time.time()


def mjd_to_iso(mjd):
    return datetime.fromtimestamp((mjd - 40587.0) * 86400.0, tz=timezone.utc).strftime("%Y-%m-%d")


def nodes_of(name):
    m = re.search(r"_(blc[0-9]+)", name)
    return m.group(1) if m else "?"


def propagate(star, mjds):
    """(ra, dec) of the star at each MJD from its J2000 catalogue position and
    space motion. Falls back to the catalogue position when motion is missing."""
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from astropy.utils import iers
    iers.conf.auto_download = False
    ra, dec = star["ra"], star["dec"]
    pmra, pmde, plx = star.get("pmra_mas_yr"), star.get("pmdec_mas_yr"), star.get("parallax_mas")
    if pmra is None or pmde is None or not plx or plx <= 0:
        return [(ra, dec)] * len(mjds)
    rv = star.get("rv_km_s")
    if rv is None or abs(rv) > RV_SANE_KMS:
        rv = 0.0
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, distance=(1000.0 / plx) * u.pc,
                 pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmde * u.mas / u.yr,
                 radial_velocity=rv * u.km / u.s, obstime=Time("J2000"))
    n = c.apply_space_motion(new_obstime=Time(np.asarray(mjds, float), format="mjd"))
    return list(zip(np.atleast_1d(n.ra.deg).tolist(), np.atleast_1d(n.dec.deg).tolist()))


def sep_arcmin(ra1, de1, ra2, de2):
    import bl_archive
    return 60.0 * bl_archive.sep_deg(ra1, de1, ra2, de2)


# ------------------------------------------------------------------ archive
def star_sessions(token, offline=False):
    """Every valid fine-res cadence of one exact archive token, as sessions with
    per-scan urls/sizes/pointings. Returns (sessions, notes) or (None, notes)."""
    import bl_archive
    import targets as T
    notes = []
    if not offline:
        polite()
    rows, trunc = T.bl_rows(token, offline=offline)
    if rows is None:
        return None, ["not in cache (offline)"]
    if trunc:
        notes.append("file-list reply truncated")
    tel = {x["url"]: x.get("telescope", "?") for x in rows}
    out = []
    good, _rej = bl_archive.cadences(token, rows, "0000")
    for c in good:
        out.append(dict(kind="named-OFF", project=c["project"], day=c["day"], fc=c["fc"],
                        telescope=tel.get(c["scans"][0]["url"], "?"),
                        scans=[dict(kind=s["kind"], target=token if s["kind"] == "ON" else token + "_OFF",
                                    name=s["name"], url=s["url"], size=s["size"], mjd=s["mjd"],
                                    ra=s["ra"], dec=s["dec"]) for s in c["scans"]]))
    ons = [x for x in bl_archive.scans(token, rows, "0000") if x["kind"] == "ON"]
    for key in sorted({(x["project"], x["day"], x["fc"]) for x in ons}):
        if any((s["project"], s["day"], s["fc"]) == key for s in out):
            continue
        g = [x for x in ons if (x["project"], x["day"], x["fc"]) == key]
        m0, m1 = min(x["mjd"] for x in g) - 2 / 1440, max(x["mjd"] for x in g) + 25 / 1440
        ck = re.sub(r"[^A-Za-z0-9._+-]", "_", f"{token}_{key[0]}_{key[1]}_{key[2]}")
        cached = (T.RAW / "bl_window" / (ck + ".json")).exists()
        if not offline and not cached:
            polite()
        win, wtrunc = T.window_rows(m0, m1, offline, ck)
        if win is None:
            notes.append(f"window {key} not in cache")
            continue
        if wtrunc:
            notes.append(f"window {key} truncated: not used")
            continue
        by_name = {r["url"].split("/")[-1]: r for r in list(win) + list(rows)}
        for c in T.abacad(token, rows, win):
            if (c["project"], c["mjd_day"], c["fc_mhz"]) != key:
                continue
            sc = []
            for nm in c["files"]:
                r = by_name[nm]
                sc.append(dict(kind="ON" if r["target"] == token else "OFF", target=r["target"],
                               name=nm, url=r["url"], size=r.get("size", 0), mjd=float(r["mjd"]),
                               ra=float(r["ra"]), dec=float(r["decl"])))
            out.append(dict(kind="ABACAD", project=key[0], day=key[1], fc=key[2],
                            telescope=c["telescope"], scans=sc))
    return out, notes


def header_coverage(url, fc, cache, offline=False):
    """(f_lo, f_hi, nchans, foff_hz, n_coarse) from one HTTP header read, cached
    per (node set, centre frequency) - the two things that fix the coverage."""
    key = f"{nodes_of(url.split('/')[-1])}|{fc}"
    if key in cache:
        return cache[key]
    if offline:
        return None
    import fsspec
    import h5py
    import hdf5plugin  # noqa: F401
    polite()
    try:
        with fsspec.open(url, block_size=2 ** 18).open() as fh:
            with h5py.File(fh, "r") as h:
                d = h["data"]
                a = dict(d.attrs)
                nch = int(a["nchans"])
                f1, df = float(a["fch1"]), float(a["foff"])
                lo, hi = sorted([f1, f1 + df * (nch - 1)])
                ncoarse = nch // 2 ** 20 if nch % 2 ** 20 == 0 else int(round(abs(df) * nch / 2.9296875))
                cache[key] = dict(f_lo_mhz=round(lo, 4), f_hi_mhz=round(hi, 4), nchans=nch,
                                  foff_hz=round(df * 1e6, 6), n_coarse=ncoarse,
                                  chunk=list(d.chunks), from_file=url.split("/")[-1])
    except Exception as e:                       # a header we cannot read is recorded, not fatal
        cache[key] = dict(error=f"{type(e).__name__}: {e}"[:200])
    return cache[key]


# ------------------------------------------------------------------ census
def band_rows(star, sessions, cov_cache, offline=False, read_headers=True):
    import targets as T
    by_band = {}
    for s in sessions:
        by_band.setdefault((s["telescope"], T.band_of(s["fc"])), []).append(s)
    out = []
    for (tel, band), ss in sorted(by_band.items()):
        days = sorted({s["day"] for s in ss})
        if len(days) < 2:
            continue
        ss = sorted(ss, key=lambda s: s["scans"][0]["mjd"])
        sess = []
        for s in ss:
            ons = [x for x in s["scans"] if x["kind"] == "ON"]
            pos = propagate(star, [x["mjd"] for x in ons])
            offs = [sep_arcmin(r, d, x["ra"], x["dec"]) for (r, d), x in zip(pos, ons)]
            fwhm = GBT_FWHM_ARCMIN_GHZ / (s["fc"] / 1000.0) if tel == "GBT" else None
            cov = header_coverage(s["scans"][0]["url"], s["fc"], cov_cache, offline) if read_headers else None
            url0 = s["scans"][0]["url"]
            sess.append(dict(kind=s["kind"], project=s["project"], mjd_day=s["day"],
                             date=mjd_to_iso(s["scans"][0]["mjd"]), fc_mhz=s["fc"],
                             nodes=nodes_of(s["scans"][0]["name"]),
                             mjd_mid=[round(x["mjd"], 5) for x in s["scans"]],
                             kinds="".join("A" if x["kind"] == "ON" else "-" for x in s["scans"]),
                             off_targets=[x["target"] for x in s["scans"] if x["kind"] == "OFF"],
                             gb=round(sum(x["size"] for x in s["scans"]) / 1e9, 2),
                             pointing_offset_arcmin=round(max(offs), 2),
                             beam_fwhm_arcmin=round(fwhm, 2) if fwhm else None,
                             pointing_ok=(max(offs) <= 0.5 * fwhm) if fwhm else None,
                             coverage=cov,
                             url_dir=url0.rsplit("/", 1)[0],
                             files=[x["name"] for x in s["scans"]]))
        good = [x for x in sess if x["pointing_ok"] is not False]
        gdays = sorted({x["mjd_day"] for x in good})
        lo = max((x["coverage"] or {}).get("f_lo_mhz", -1) for x in sess) if read_headers else None
        hi = min((x["coverage"] or {}).get("f_hi_mhz", 1e9) for x in sess) if read_headers else None
        rv = star.get("rv_km_s")
        out.append(dict(star=star["star"], token=star["bl_token"], common=star.get("common"),
                        distance_ly=star.get("distance_ly"), telescope=tel, band=band,
                        n_sessions=len(sess), n_days=len(days), n_days_pointing_ok=len(gdays),
                        mjd_days=days, first=sess[0]["date"], last=sess[-1]["date"],
                        span_days=days[-1] - days[0], total_gb=round(sum(x["gb"] for x in sess), 1),
                        one_coarse_gb=round(len(sess) * 6 * COARSE_BYTES / 1e9, 2),
                        common_coverage_mhz=[lo, hi] if lo and hi and lo > 0 and hi < 1e9 and hi > lo else None,
                        pointing_ok=all(x["pointing_ok"] is not False for x in sess),
                        rv_km_s=rv, rv_ok=(rv is not None and abs(rv) <= RV_SANE_KMS),
                        pm_mas_yr=(round(math.hypot(star["pmra_mas_yr"], star["pmdec_mas_yr"]), 1)
                                   if star.get("pmra_mas_yr") is not None and star.get("pmdec_mas_yr") is not None
                                   else None),
                        parallax_mas=star.get("parallax_mas"),
                        secular_hz_per_yr_at_1420=secular_hz_per_yr(star),
                        same_simbad_object_as=star.get("same_simbad_object_as"),
                        sessions=sess))
    return out


def secular_hz_per_yr(star, f_hz=1420.405751768e6):
    """Perspective (secular) radial acceleration mu^2 d, as a frequency drift."""
    pmra, pmde, plx = star.get("pmra_mas_yr"), star.get("pmdec_mas_yr"), star.get("parallax_mas")
    if pmra is None or pmde is None or not plx or plx <= 0:
        return None
    mu = math.radians(math.hypot(pmra, pmde) / 3.6e6)          # rad / yr
    d_m = 1000.0 / plx * 3.0856775814913673e16
    yr = 365.25 * 86400.0
    a = (mu / yr) ** 2 * d_m                                    # m / s^2
    return round(f_hz * a * yr / 299792458.0, 3)


def run(unchecked=False, offline=False, read_headers=True, only=None):
    import targets as T
    data = json.loads((HERE / "targets.json").read_text())
    stars = data["stars"]
    RAW.mkdir(parents=True, exist_ok=True)
    cov_path = RAW / "headers.json"
    cov_cache = json.loads(cov_path.read_text()) if cov_path.exists() else {}
    pick, why = [], {}
    seen_tokens = set()
    for s in stars:
        tok = s.get("bl_token")
        if not tok:
            continue
        if only and s["star"] not in only and tok not in only:
            continue
        st = T.cadence_state(s)
        bl = s.get("bl", {})
        cads = (bl.get("cadences") or []) + (bl.get("abacad") or [])
        n_band = {}
        for c in cads:
            n_band.setdefault(c["band"], set()).add(c["mjd_day"])
        multi = any(len(v) >= 2 for v in n_band.values())
        if tok == "GJ699" or only:
            pick.append(s)
            why[s["star"]] = "explicit"
        elif st == "yes" and multi:
            pick.append(s)
            why[s["star"]] = "targets.json shows >= 2 same-band dates"
        elif st == "unchecked" and unchecked:
            pick.append(s)
            why[s["star"]] = "never file-checked by targets.py"
    print(f"census: {len(pick)} stars to (re)check"
          f"{' (offline)' if offline else ' (<= 2 requests/s, replies cached)'}")
    rows, skipped = [], []
    t0 = time.time()
    for i, s in enumerate(pick):
        tok = s["bl_token"]
        try:
            sess, notes = star_sessions(tok, offline)
        except Exception as e:
            skipped.append(dict(star=s["star"], why=f"{type(e).__name__}: {e}"[:160]))
            continue
        if sess is None:
            skipped.append(dict(star=s["star"], why="; ".join(notes)))
            continue
        if tok in seen_tokens:
            # two catalogue rows, one archive name (e.g. GJ725 A/B): keep both rows, flag the twin
            pass
        seen_tokens.add(tok)
        br = band_rows(s, sess, cov_cache, offline, read_headers)
        for r in br:
            r["notes"] = notes
            r["why_checked"] = why[s["star"]]
        rows += br
        if tok == "GJ699":
            gj699 = dict(sessions=[dict(project=x["project"], day=x["day"], fc=x["fc"], kind=x["kind"],
                                        band=T.band_of(x["fc"])) for x in sess])
        if (i + 1) % 25 == 0:
            cov_path.write_text(json.dumps(cov_cache, indent=0))
            print(f"  {i+1}/{len(pick)} stars, {len(rows)} multi-epoch rows, {time.time()-t0:.0f} s")
    cov_path.write_text(json.dumps(cov_cache, indent=0))
    rows.sort(key=lambda r: ((r["distance_ly"] or 1e9), r["total_gb"], r["band"]))
    for k, r in enumerate(rows, 1):
        r["rank"] = k
    return rows, skipped, locals().get("gj699"), pick


CSV_COLS = ["rank", "star", "token", "common", "distance_ly", "telescope", "band", "n_days",
            "n_days_pointing_ok", "n_sessions", "mjd_days", "first", "last", "span_days", "total_gb",
            "one_coarse_gb", "common_coverage_mhz", "pointing_ok", "max_pointing_offset_arcmin",
            "rv_km_s", "rv_ok", "pm_mas_yr", "secular_hz_per_yr_at_1420", "same_simbad_object_as"]


def write(rows, skipped, gj699, pick, unchecked):
    with open(OUT_CSV, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(CSV_COLS)
        for r in rows:
            vals = []
            for c in CSV_COLS:
                if c == "max_pointing_offset_arcmin":
                    v = max(x["pointing_offset_arcmin"] for x in r["sessions"])
                else:
                    v = r.get(c)
                if isinstance(v, list):
                    v = " ".join(str(x) for x in v)
                vals.append("" if v is None else v)
            w.writerow(vals)
    stars_ok = sorted({r["star"] for r in rows if r["n_days_pointing_ok"] >= 2})
    meta = dict(
        what="stars with >= 2 valid fine-res (.0000) cadences on different MJD days in the same band",
        evidence="archive metadata only - no signal search, no E-level",
        rules=dict(epoch="one (project, MJD day, centre freq) session holding a valid cadence "
                         "(bl_archive.cadences named-OFF, or targets.abacad ABACAD)",
                   same_band="telescope + band letter from the centre frequency (L 1-2, S 2-4, C 4-8, X 8-12 GHz)",
                   pointing_ok="max ON-scan pointing offset from the star propagated to the scan date "
                               "<= half the GBT FWHM (12.45'/f_GHz)",
                   rv_ok=f"|SIMBAD RV| <= {RV_SANE_KMS:.0f} km/s"),
        checked=dict(n_stars=len(pick), unchecked_included=unchecked,
                     scope=("targets.json's file-checked stars showing >= 2 same-band dates, GJ699 always"
                            + (", and every star targets.json had never file-checked" if unchecked else ""))),
        n_rows=len(rows), n_stars=len({r["star"] for r in rows}),
        n_stars_pointing_ok=len(stars_ok),
        gj699=gj699, skipped=skipped,
        when=datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        source="Breakthrough Listen open data (Lebofsky et al. 2019, PASP 131, 124505)")
    OUT_JSON.write_text(json.dumps(dict(meta=meta, rows=rows), indent=1))
    print(f"wrote {OUT_CSV.name} ({OUT_CSV.stat().st_size/1e3:.0f} kB), "
          f"{OUT_JSON.name} ({OUT_JSON.stat().st_size/1e3:.0f} kB)")


def show(rows, gj699, n=25):
    print(f"\n{'rk':>3} {'star':<10} {'ly':>6} {'band':<5} {'days':>4} {'span':>5} {'GB':>6} "
          f"{'1-coarse':>8} {'point':>6} {'RV':>8}  MJD days")
    for r in rows[:n]:
        print(f"{r['rank']:>3} {r['star']:<10} {r['distance_ly'] or 0:6.2f} {r['telescope'][:3]}-{r['band']:<1} "
              f"{r['n_days']:>4} {r['span_days']:>5} {r['total_gb']:>6.0f} {r['one_coarse_gb']:>8.2f} "
              f"{'ok' if r['pointing_ok'] else 'OFF!':>6} {r['rv_km_s'] if r['rv_km_s'] is not None else '-':>8}  "
              f"{' '.join(str(d) for d in r['mjd_days'])}")
    if gj699:
        print("\nGJ699 (Barnard's Star), every valid fine-res session in the archive:")
        for s in gj699["sessions"]:
            print(f"  {s['band']}-band  MJD {s['day']}  {s['fc']} MHz  {s['project']}  ({s['kind']})")
        bands = {}
        for s in gj699["sessions"]:
            bands.setdefault(s["band"], set()).add(s["day"])
        two = [b for b, d in bands.items() if len(d) >= 2]
        print("  -> a second same-band epoch: " + (", ".join(two) if two else
              "NONE (one L-band date, one S-band date)"))


def selftest():
    checks = []
    star = dict(star="X", bl_token="X", ra=10.0, dec=5.0, pmra_mas_yr=0.0, pmdec_mas_yr=0.0,
                parallax_mas=100.0, rv_km_s=10.0, distance_ly=32.6)

    def sess(day, fc, kind="ABACAD", dec_on=5.0, proj="P1"):
        sc = []
        for i in range(6):
            on = i % 2 == 0
            sc.append(dict(kind="ON" if on else "OFF", target="X" if on else f"H{i}", name=f"blc00_{day}_{i}.0000.h5",
                           url=f"https://x/pipeline/{proj}/holding/spliced_blc00_guppi_{day}_{i}.0000.h5",
                           size=int(15e9), mjd=day + 0.5 + i * 5 / 1440, ra=10.0, dec=dec_on if on else 7.0))
        return dict(kind=kind, project=proj, day=day, fc=fc, telescope="GBT", scans=sc)
    r = band_rows(star, [sess(57500, 1475.1), sess(57600, 1475.1, proj="P2"), sess(57700, 2300.4)],
                  {}, offline=True, read_headers=False)
    checks.append(("two L dates qualify, a lone S date does not", len(r) == 1 and r[0]["band"] == "L"
                   and r[0]["n_days"] == 2))
    r = band_rows(star, [sess(57500, 1475.1), sess(57500, 1475.1, proj="P2")], {}, True, False)
    checks.append(("two sessions on ONE day are one date, not two epochs", r == []))
    r = band_rows(star, [sess(57500, 1475.1), sess(57600, 1475.1, dec_on=6.2)], {}, True, False)
    checks.append(("an ON 1.2 deg off the star is flagged", len(r) == 1 and not r[0]["pointing_ok"]
                   and r[0]["n_days_pointing_ok"] == 1))
    b = dict(pmra_mas_yr=-801.551, pmdec_mas_yr=10362.394, parallax_mas=546.9759)
    sh = secular_hz_per_yr(b)
    checks.append((f"Barnard's Star secular term {sh:.2f} Hz/yr at 1420 MHz (expect ~21.5)", 21.0 < sh < 22.0))
    p = propagate(dict(ra=269.4520769586187, dec=4.693364966576667, rv_km_s=-110.11, **b), [57457.556])[0]
    off = sep_arcmin(p[0], p[1], 269.4483291667, 4.7391686111)     # the archive's 2016 pointing
    checks.append((f"Barnard's Star propagated to 2016-03-10 is {off:.2f}' from BL's pointing "
                   f"(J2000 position would be 2.8')", off < 0.2))
    checks.append(("date conversion", mjd_to_iso(57457.556) == "2016-03-10"))
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    if sys.argv[1:2] == ["selftest"]:
        return selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--unchecked", action="store_true", help="also stars targets.json never file-checked")
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--no-headers", action="store_true", help="skip the coverage header reads")
    ap.add_argument("--only", nargs="*", help="just these stars/tokens (no file is written)")
    a = ap.parse_args()
    rows, skipped, gj699, pick = run(a.unchecked, a.offline, not a.no_headers, a.only)
    if not a.only:
        write(rows, skipped, gj699, pick, a.unchecked)
    show(rows, gj699)
    if skipped:
        print(f"\nskipped {len(skipped)}: " + "; ".join(f"{x['star']} ({x['why']})" for x in skipped[:10]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
