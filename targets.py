#!/usr/bin/env python3
"""targets.py - the Breakthrough Listen star list, as one diffable file.

WHERE to look matters as much as HOW. This builds targets.json: every star of
BL's canonical nearby-star sample with its catalogue numbers, whether Earth
can be seen TRANSITING the Sun from it, the EXACT name the BL open-data
archive files it under, whether the archive holds a usable fine-resolution
ON/OFF cadence for it, and a simple priority score.

    python targets.py build             # fetch (cached in data/targets_raw/) + write targets.json
    python targets.py build --offline   # rebuild from the cache only
    python targets.py show GJ699        # one star (any name: GJ699, HIP87937, "Barnard")
    python targets.py top 20            # highest priority
    python targets.py etz               # Earth-Transit-Zone stars and their BL data
    python targets.py check             # cross-check a few well-known stars
    python targets.py selftest          # offline checks, no network

SOURCES
  * Star list: Isaacson et al. 2017, PASP 129, 054501, "The Breakthrough Listen
    Search for Intelligent Life: Target Selection of Nearby Stars and Galaxies".
    VizieR catalogue J/PASP/129/E4501, table1 (1,709 stars). The galaxies of the
    same paper are not in table1 and are not included here.
  * Astrometry, spectral type, radial velocity: SIMBAD (TAP), J2000/ICRS.
  * Archive names and files: BL open data (Lebofsky et al. 2019, PASP 131,
    124505), via bl_archive.py (exact names only - the archive matches by prefix).

EARTH TRANSIT ZONE (Heller & Pudritz 2016, Astrobiology 16, 259)
  From a star at ecliptic latitude |beta| < ~0.264 deg (the paper's value;
  arcsin(R_sun / 1 au) = 0.2665 deg; `build` reports any star in between), Earth is seen to cross the Sun once a year. Observers there could have found
  us the way we find exoplanets, which makes them a natural place to look for
  a deliberate signal. We flag  etz = |beta| < 0.264 deg  and the looser
  near_etz = |beta| < 2 deg, with beta in astropy's BarycentricTrueEcliptic
  frame (equinox J2000). Not implemented (only cited, because we could not
  verify it here): the moving ETZ over +-5,000 yr (Kaltenegger & Faherty 2021,
  Nature 594, 505) and the SETI-ellipsoid timing idea (Davenport et al. 2022,
  AJ 164, 117; Nilipour et al. 2023, AJ 166, 79), which needs a reference event
  and its distance error, not just a position.

PRIORITY SCORE (simple on purpose; every term is in targets.json)
    score = 40*etz + 10*near_etz_only + 25*bl_cadence + 20*closeness - 30*searched
      etz            1 if |beta| < 0.264 deg
      near_etz_only  1 if 0.264 <= |beta| < 2 deg
      bl_cadence     1 if the archive holds a valid fine-res (0000) cadence: either
                     bl_archive.cadences (OFFs named <name>_OFF) or abacad() below
                     (ON X ON Y ON Z, X/Y/Z other stars >= 0.2 deg away, same
                     session, band and product, gaps <= 20 min)
      closeness      max(0, 1 - distance_ly / 100)
      searched       1 if this repo already has an E1-or-better search (EVIDENCE.md);
                     VOID and E0 results do not count
  Stars whose archive cadence was NOT checked score bl_cadence = 0; the
  `bl.checked` / `bl.abacad_checked` fields say which.

WHAT WAS CHECKED (2026-11-06 build): archive names for all 1,709 stars (one
list-targets call); file lists for 375 stars (every near-ETZ star, the nearest
100, and every star with a <name>_OFF target); ABACAD time-window queries for
the near-ETZ stars and the nearest 100 (126 stars with fine-res ON scans).
Everything is cached in data/targets_raw/ (gitignored); `build --offline`
rebuilds targets.json from it byte for byte.
"""
import argparse
import csv
import io
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "targets.json"
RAW = HERE / "data" / "targets_raw"

CATALOGUE = "J/PASP/129/E4501/table1"
CATALOGUE_REF = "Isaacson et al. 2017, PASP 129, 054501"
VIZIER = ["https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
          "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"]
SIMBAD_TAP = ["https://simbad.cds.unistra.fr/simbad/sim-tap/sync",
              "https://simbad.harvard.edu/simbad/sim-tap/sync"]
BL_LIST = "http://seti.berkeley.edu/opendata/api/list-targets"

ETZ_DEG = 0.264          # Heller & Pudritz 2016 half-width; geometric asin(R_sun/1 au) = 0.2665
NEAR_ETZ_DEG = 2.0
LY_PER_PC = 3.261563777

# This repo's own searches, from EVIDENCE.md. Only E1 or better counts as searched.
OUR_SEARCHES = {
    "GJ699": "E2",        # GJ699_CADENCE_RESULT.md, full fine-res cadence 1.1-1.9 GHz
    "GJ144": "VOID",      # EPSILON_ERIDANI_RESULT.md (mid-res drift unresolved)
    "GJ411": "VOID", "GJ273": "VOID", "GJ15A": "VOID", "GJ71": "VOID",
    "GJ876": "VOID",      # SWEEP_RESULTS.md, single mid-res ON files
}
SEARCHED_LEVELS = ("E1", "E2", "E3")

W = dict(etz=40, near_etz_only=10, bl_cadence=25, closeness=20, searched=-30)

BANDS = [("L", 1000, 2000), ("S", 2000, 4000), ("C", 4000, 8000), ("X", 8000, 12000)]


def band_of(fc_mhz):
    for name, lo, hi in BANDS:
        if lo <= fc_mhz < hi:
            return name
    return "?"


# ----------------------------------------------------------------- fetching
def _get(url, data=None, timeout=120):
    req = urllib.request.Request(url, data=data, headers={"User-Agent": "seti-cookbook targets.py"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8")


_DEAD = set()


def _first_ok(urls, make, what):
    err = None
    for u in urls:
        if u in _DEAD and u != urls[-1]:
            continue
        try:
            return make(u)
        except Exception as e:          # mirror down: try the next one, and skip it from now on
            err = e
            _DEAD.add(u)
            print(f"  {what}: {u.split('/')[2]} failed ({type(e).__name__}); trying next mirror")
    raise RuntimeError(f"{what}: every mirror failed: {err}")


def fetch_isaacson(offline=False):
    path = RAW / "isaacson2017_table1.tsv"
    if path.exists() or offline:
        return path.read_text()
    q = urllib.parse.urlencode({"-source": CATALOGUE, "-out.max": "unlimited",
                                "-out.all": "", "-out.add": "_RA,_DE"})
    txt = _first_ok(VIZIER, lambda u: _get(u + "?" + q), "VizieR")
    RAW.mkdir(parents=True, exist_ok=True)
    path.write_text(txt)
    return txt


def parse_isaacson(txt):
    """VizieR tsv -> list of dicts (Star, ra/dec deg, Vmag, SpType, Dist pc, pm arcsec/yr)."""
    lines = [l for l in txt.splitlines() if l and not l.startswith("#")]
    head = lines[0].split("\t")
    out = []
    for l in lines[3:]:                 # header, units, dashes
        v = dict(zip(head, [x.strip() for x in l.split("\t")]))
        if not v.get("Star"):
            continue

        def f(k):
            try:
                return float(v.get(k, ""))
            except ValueError:
                return None
        out.append(dict(star=v["Star"], ra=f("_RAJ2000"), dec=f("_DEJ2000"), vmag=f("Vmag"),
                        sptype=v.get("SpType") or None, dist_pc=f("Dist"),
                        pmra=f("pmRA"), pmde=f("pmDE"), simbad_name=v.get("SimbadName") or v["Star"]))
    return out


def simbad_id(name):
    """Isaacson's compact names -> SIMBAD identifiers (HIP2 -> 'HIP 2', GJ15A -> 'GJ 15 A')."""
    m = re.fullmatch(r"HIP(\d+)", name)
    if m:
        return f"HIP {int(m.group(1))}"
    m = re.fullmatch(r"GJ(\d+(?:\.\d+)?)([A-Z]*)", name)
    if m:
        return f"GJ {m.group(1)}" + (f" {m.group(2)}" if m.group(2) else "")
    m = re.fullmatch(r"(LHS|LP)(\S+)", name)
    if m:
        return f"{m.group(1)} {m.group(2)}"
    return name


def _tap(query, fmt="csv"):
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL",
                                   "FORMAT": fmt, "QUERY": query}).encode()
    return _first_ok(SIMBAD_TAP, lambda u: _get(u, data=data, timeout=180), "SIMBAD")


def _q(s):
    return "'" + s.replace("'", "''") + "'"


def fetch_simbad(stars, offline=False, batch=300):
    """Two cached TAP queries per batch: the basic row for each id, then the
    HIP/GJ/NAME identifiers of every matched object."""
    RAW.mkdir(parents=True, exist_ok=True)
    ids = sorted({simbad_id(s["star"]) for s in stars})
    basic, idents = [], []
    for k in range(0, len(ids), batch):
        chunk = ids[k:k + batch]
        p1 = RAW / f"simbad_basic_{k // batch:02d}.csv"
        p2 = RAW / f"simbad_ident_{k // batch:02d}.csv"
        if not p1.exists():
            if offline:
                raise SystemExit(f"offline and {p1} is not cached")
            q = ("SELECT i.id AS qid, b.oid, b.main_id, b.ra, b.dec, b.pmra, b.pmdec, "
                 "b.plx_value, b.rvz_radvel, b.sp_type FROM ident AS i JOIN basic AS b "
                 "ON b.oid = i.oidref WHERE i.id IN (" + ",".join(_q(x) for x in chunk) + ")")
            p1.write_text(_tap(q))
            time.sleep(1.0)
        rows = list(csv.DictReader(io.StringIO(p1.read_text())))
        basic += rows
        if not p2.exists():
            if offline:
                raise SystemExit(f"offline and {p2} is not cached")
            oids = sorted({r["oid"] for r in rows})
            q = ("SELECT oidref, id FROM ident WHERE oidref IN (" + ",".join(oids) + ") AND "
                 "(id LIKE 'HIP %' OR id LIKE 'GJ %' OR id LIKE 'NAME %' OR id LIKE '* %' "
                 "OR id LIKE 'HD %')")
            p2.write_text(_tap(q))
            time.sleep(1.0)
        idents += list(csv.DictReader(io.StringIO(p2.read_text())))
    by_q = {r["qid"]: r for r in basic}
    ids_by_oid = {}
    for r in idents:
        ids_by_oid.setdefault(r["oidref"], []).append(re.sub(r"\s+", " ", r["id"]).strip())
    # Names SIMBAD does not know in Isaacson's compact form (SO0253, DEN1048, GJ725 ...):
    # one cached cone query at the catalogue position, nearest object within 30".
    miss = [x for x in stars if simbad_id(x["star"]) not in by_q]
    p3, p4 = RAW / "simbad_cone.csv", RAW / "simbad_cone_ident.csv"
    if miss and not p3.exists():
        if offline:
            raise SystemExit(f"offline and {p3} is not cached")
        circ = " OR ".join(f"CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', {x['ra']}, {x['dec']}, "
                           f"{CONE_DEG}))=1" for x in miss)
        q = ("SELECT b.oid, b.main_id, b.otype, b.ra, b.dec, b.pmra, b.pmdec, b.plx_value, "
             "b.rvz_radvel, b.sp_type FROM basic AS b WHERE (" + circ + ")")
        p3.write_text(_tap(q))
        time.sleep(1.0)
    cone = list(csv.DictReader(io.StringIO(p3.read_text()))) if miss else []
    for x in miss:
        best = None
        for r in cone:
            if r.get("otype", "").startswith("Pl") or r.get("otype") == "err":   # planets (and retracted ones) sit on their star
                continue
            try:
                d = _sep(x["ra"], x["dec"], float(r["ra"]), float(r["dec"]))
            except (TypeError, ValueError):
                continue
            if d <= CONE_DEG and (best is None or d < best[0]):
                best = (d, r)
        if best:
            by_q[simbad_id(x["star"]) + "|" + x["star"] + f"|{x['ra']:.4f}"] = dict(
                best[1], qid=None, by_position_arcsec=round(best[0] * 3600, 1))
    if cone and not p4.exists():
        if offline:
            raise SystemExit(f"offline and {p4} is not cached")
        oids = sorted({r["oid"] for r in cone})
        q = ("SELECT oidref, id FROM ident WHERE oidref IN (" + ",".join(oids) + ") AND "
             "(id LIKE 'HIP %' OR id LIKE 'GJ %' OR id LIKE 'NAME %' OR id LIKE '* %' "
             "OR id LIKE 'HD %')")
        p4.write_text(_tap(q))
        time.sleep(1.0)
    if cone:
        for r in csv.DictReader(io.StringIO(p4.read_text())):
            ids_by_oid.setdefault(r["oidref"], []).append(re.sub(r"\s+", " ", r["id"]).strip())
    return by_q, ids_by_oid


CONE_DEG = 30 / 3600


def _sep(ra1, de1, ra2, de2):
    r = math.radians
    c = (math.sin(r(de1)) * math.sin(r(de2))
         + math.cos(r(de1)) * math.cos(r(de2)) * math.cos(r(ra1 - ra2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def simbad_row(by_q, s):
    r = by_q.get(simbad_id(s["star"]))
    if r:
        return r, "SIMBAD (by identifier)"
    r = by_q.get(simbad_id(s["star"]) + "|" + s["star"] + f"|{s['ra']:.4f}")
    if r:
        return r, f"SIMBAD (nearest within 30 arcsec of the catalogue position: {r['by_position_arcsec']}\")"
    return None, "Isaacson2017 (SIMBAD not resolved)"


def fetch_bl_list(offline=False):
    path = RAW / "bl_list_targets.json"
    if not path.exists():
        if offline:
            return None
        path.write_text(_get(BL_LIST, timeout=90))
    return set(json.loads(path.read_text()))


def bl_rows(token, offline=False, limit=20000):
    """Prefix-matched archive reply for one token, cached. Returns (rows, truncated)."""
    import bl_archive
    d = RAW / "bl"
    d.mkdir(parents=True, exist_ok=True)
    path = d / (re.sub(r"[^A-Za-z0-9._+-]", "_", token) + ".json")
    if not path.exists():
        if offline:
            return None, False
        rows = bl_archive.query(token, limit=limit)
        path.write_text(json.dumps(rows))
        time.sleep(0.5)
    rows = json.loads(path.read_text())
    return rows, len(rows) >= limit


# ----------------------------------------------------------------- physics
def ecliptic_lat(ra_deg, dec_deg):
    """Barycentric true ecliptic latitude (deg, equinox J2000). Vectorised."""
    import numpy as np
    from astropy import units as u
    from astropy.coordinates import SkyCoord, BarycentricTrueEcliptic
    from astropy.time import Time
    c = SkyCoord(ra=np.atleast_1d(ra_deg) * u.deg, dec=np.atleast_1d(dec_deg) * u.deg, frame="icrs")
    e = c.transform_to(BarycentricTrueEcliptic(equinox=Time("J2000")))
    return e.lat.deg


def ecliptic_to_icrs(lon_deg, lat_deg):
    from astropy import units as u
    from astropy.coordinates import SkyCoord, BarycentricTrueEcliptic
    from astropy.time import Time
    e = SkyCoord(lon=lon_deg * u.deg, lat=lat_deg * u.deg,
                 frame=BarycentricTrueEcliptic(equinox=Time("J2000")))
    i = e.icrs
    return float(i.ra.deg), float(i.dec.deg)


def etz_window(beta, rate_per_century, horizon=5000):
    """[enter, exit] calendar years of the ETZ crossing for linear motion in beta
    (J2000 = year 2000), clipped to +-horizon years; None if never inside."""
    if abs(rate_per_century) < 1e-12:
        return [None, None] if abs(beta) < ETZ_DEG else None
    t1 = (-ETZ_DEG - beta) / rate_per_century * 100
    t2 = (ETZ_DEG - beta) / rate_per_century * 100
    lo, hi = sorted((t1, t2))
    if hi < -horizon or lo > horizon:
        return None
    return [round(2000 + lo) if lo > -horizon else None, round(2000 + hi) if hi < horizon else None]


def etz_flags(beta):
    return dict(etz=bool(abs(beta) < ETZ_DEG), near_etz=bool(abs(beta) < NEAR_ETZ_DEG))


# ----------------------------------------------------------------- archive token
def token_candidates(star, ids):
    """Names the archive might use: Isaacson's own first, then HIP/GJ forms."""
    c = [star]
    for x in ids:
        m = re.fullmatch(r"(HIP|GJ) (\d+(?:\.\d+)?)(?: ?([A-Z]+))?", x)
        if m:
            c.append(m.group(1) + m.group(2) + (m.group(3) or ""))
    seen, out = set(), []
    for x in c:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def resolve_token(cands, names):
    """First candidate that is EXACTLY an archive target name. A name that only
    starts some archive names (HIP2 -> HIP28614) is never accepted."""
    for t in cands:
        if t in names:
            return t
    return None


def bl_check(token, names, offline=False):
    """Archive facts for one exact token: ON/OFF scan counts per product, and
    every valid fine-res cadence (bl_archive rules)."""
    import bl_archive
    rows, trunc = bl_rows(token, offline=offline)
    if rows is None:
        return dict(checked=False)
    on, off, dropped = bl_archive.exact_rows(token, rows)
    good, rejected = bl_archive.cadences(token, rows, "0000")
    tel = {x["url"]: x.get("telescope", "?") for x in rows}
    cads = []
    for c in good:
        cads.append(dict(project=c["project"], mjd_day=c["day"], fc_mhz=c["fc"],
                         band=band_of(c["fc"]), telescope=tel.get(c["scans"][0]["url"], "?"),
                         gb=round(sum(s["size"] for s in c["scans"]) / 1e9, 1),
                         max_gap_min=round(c["max_gap_min"], 1),
                         off_sep_deg=round(c["off_sep_deg"], 2),
                         files=[s["name"] for s in c["scans"]]))
    n0 = len(bl_archive.scans(token, rows, "0000"))
    return dict(checked=True, rows=len(rows), truncated=trunc, exact_on=len(on),
                exact_off=len(off), prefix_dropped=len(dropped),
                off_name_listed=(token + "_OFF") in names,
                fine_scans=n0, cadences=cads, n_rejected_groups=len(rejected))


# --------------------------------------------- ABACAD cadences, OFFs = other stars
# Most BL Green Bank cadences of 2016-17 point ON the star three times and, in
# between, at three OTHER catalogue stars (Enriquez et al. 2017: "ABACAD").
# Those OFF scans are filed under the other stars' names, so bl_archive.cadences
# (which needs <name>_OFF) cannot see them (KNOWN_ISSUES #6). The archive does
# answer a TIME-window query with target="%" (SQL wildcard: every name), so we
# ask for exactly the minutes of each session and rebuild the sequence.
BL_FILES = "http://seti.berkeley.edu/opendata/api/query-files"


def window_rows(mjd0, mjd1, offline=False, cache_key=None):
    d = RAW / "bl_window"
    d.mkdir(parents=True, exist_ok=True)
    path = d / ((cache_key or f"{mjd0:.5f}_{mjd1:.5f}") + ".json")
    if not path.exists():
        if offline:
            return None
        q = urllib.parse.urlencode({"target": "%", "time-start": f"{mjd0:.6f}",
                                    "time-end": f"{mjd1:.6f}", "file-types": "HDF5", "limit": 2000})
        rows = json.loads(_get(BL_FILES + "?" + q, timeout=90)).get("data", [])
        path.write_text(json.dumps(rows))
        time.sleep(0.5)
    return json.loads(path.read_text())


def abacad(token, star_rows, window, prod="0000"):
    """ON X ON Y ON Z cadences of one session/band/product, X/Y/Z any other
    target names pointed >= MIN_OFF_DEG away, gaps <= MAX_GAP_MIN. `window`
    is every archive row (any name) in the session's time span."""
    import bl_archive
    ons = [x for x in bl_archive.scans(token, star_rows, prod) if x["kind"] == "ON"]
    found = []
    for (proj, day, fc) in sorted({(x["project"], x["day"], x["fc"]) for x in ons}):
        on_g = [x for x in ons if (x["project"], x["day"], x["fc"]) == (proj, day, fc)]
        ra0 = sorted(x["ra"] for x in on_g)[len(on_g) // 2]
        de0 = sorted(x["dec"] for x in on_g)[len(on_g) // 2]
        seq = {}
        for r in window:
            name = r["url"].split("/")[-1]
            if bl_archive.product(name) != prod:
                continue
            pj, dy = bl_archive.session(r["url"])
            if (pj, dy, round(float(r["center_freq"]), 1)) != (proj, day, fc):
                continue
            key = (r["target"], bl_archive.scan_no(name))
            seq.setdefault(key, dict(kind="ON" if r["target"] == token else "OFF", target=r["target"],
                                     name=name, mjd=float(r["mjd"]), ra=float(r["ra"]),
                                     dec=float(r["decl"]), size=r.get("size", 0),
                                     telescope=r.get("telescope", "?")))
        g = sorted(seq.values(), key=lambda x: x["mjd"])
        for i in range(len(g) - 5):
            w = g[i:i + 6]
            if [x["kind"] for x in w] != ["ON", "OFF"] * 3:
                continue
            gaps = [(b["mjd"] - a["mjd"]) * 1440 for a, b in zip(w, w[1:])]
            seps = [bl_archive.sep_deg(ra0, de0, x["ra"], x["dec"]) for x in w if x["kind"] == "OFF"]
            if max(gaps) > bl_archive.MAX_GAP_MIN or min(seps) < bl_archive.MIN_OFF_DEG:
                continue
            found.append(dict(kind="ABACAD-other-OFFs", project=proj, mjd_day=day, fc_mhz=fc,
                              band=band_of(fc), telescope=w[0]["telescope"],
                              gb=round(sum(x["size"] for x in w) / 1e9, 1),
                              max_gap_min=round(max(gaps), 1), off_sep_deg=round(min(seps), 2),
                              off_targets=[x["target"] for x in w if x["kind"] == "OFF"],
                              files=[x["name"] for x in w]))
            break
    return found


def abacad_check(token, offline=False):
    import bl_archive
    rows, _ = bl_rows(token, offline=offline)
    if rows is None:
        return None
    ons = [x for x in bl_archive.scans(token, rows, "0000") if x["kind"] == "ON"]
    out = []
    for key in sorted({(x["project"], x["day"], x["fc"]) for x in ons}):
        g = [x for x in ons if (x["project"], x["day"], x["fc"]) == key]
        m0, m1 = min(x["mjd"] for x in g) - 2 / 1440, max(x["mjd"] for x in g) + 25 / 1440
        ck = re.sub(r"[^A-Za-z0-9._+-]", "_", f"{token}_{key[0]}_{key[1]}_{key[2]}")
        win = window_rows(m0, m1, offline, ck)
        if win is None:
            return None
        out += abacad(token, rows, win)
    return out


# ----------------------------------------------------------------- score
def score(s):
    beta = s.get("ecl_lat_deg")
    etz = bool(s.get("etz"))
    near_only = bool(s.get("near_etz")) and not etz
    cad = bool(s.get("bl", {}).get("cadences")) or bool(s.get("bl", {}).get("abacad"))
    d = s.get("distance_ly")
    close = max(0.0, 1.0 - d / 100.0) if d else 0.0
    searched = OUR_SEARCHES.get(s.get("bl_token") or s["star"]) in SEARCHED_LEVELS \
        or OUR_SEARCHES.get(s["star"]) in SEARCHED_LEVELS
    parts = dict(etz=W["etz"] * etz, near_etz_only=W["near_etz_only"] * near_only,
                 bl_cadence=W["bl_cadence"] * cad, closeness=round(W["closeness"] * close, 2),
                 searched=W["searched"] * searched)
    return round(sum(parts.values()), 2), parts


# ----------------------------------------------------------------- build
def build(offline=False, check_all=False):
    import numpy as np
    stars = parse_isaacson(fetch_isaacson(offline))
    print(f"{CATALOGUE}: {len(stars)} stars")
    by_q, ids_by_oid = fetch_simbad(stars, offline)
    out = []
    for s in stars:
        r, src = simbad_row(by_q, s)
        e = dict(star=s["star"], isaacson=dict(ra=s["ra"], dec=s["dec"], vmag=s["vmag"],
                                               sptype=s["sptype"], dist_pc=s["dist_pc"],
                                               pmra_arcsec_yr=s["pmra"], pmde_arcsec_yr=s["pmde"]))

        def f(k):
            try:
                return float(r[k]) if r and r[k] not in ("", None) else None
            except ValueError:
                return None
        ids = ids_by_oid.get(r["oid"], []) if r else []
        e.update(simbad_main_id=r["main_id"] if r else None,
                 hip=next((x for x in ids if x.startswith("HIP ")), None),
                 gj=sorted({x for x in ids if x.startswith("GJ ")}) or None,
                 simbad_oid=int(r["oid"]) if r else None,
                 common=min((x[5:] for x in ids if x.startswith("NAME ")), key=len, default=None),
                 ra=f("ra") if r else s["ra"], dec=f("dec") if r else s["dec"],
                 pmra_mas_yr=f("pmra"), pmdec_mas_yr=f("pmdec"), parallax_mas=f("plx_value"),
                 rv_km_s=f("rvz_radvel"), sptype=(r["sp_type"] or None) if r else s["sptype"],
                 position_source=src)
        if e["parallax_mas"] and e["parallax_mas"] > 0:
            e["distance_ly"] = round(1000.0 / e["parallax_mas"] * LY_PER_PC, 3)
            e["distance_source"] = "SIMBAD parallax"
        elif s["dist_pc"]:
            e["distance_ly"] = round(s["dist_pc"] * LY_PER_PC, 3)
            e["distance_source"] = "Isaacson2017 Dist"
        else:
            e["distance_ly"], e["distance_source"] = None, None
        e["_ids"] = ids
        out.append(e)
    by_oid = {}
    for e in out:
        if e["simbad_oid"]:
            by_oid.setdefault(e["simbad_oid"], []).append(e["star"])
    for e in out:
        twins = [x for x in by_oid.get(e["simbad_oid"], []) if x != e["star"]]
        e["same_simbad_object_as"] = twins or None
    beta = ecliptic_lat(np.array([e["ra"] for e in out]), np.array([e["dec"] for e in out]))
    # proper motion moves a star across the ETZ edge: beta 100 years on (linear, small angle)
    cosd = np.cos(np.radians([e["dec"] for e in out]))
    pmra = np.array([(e["pmra_mas_yr"] or 0.0) for e in out]) / 3.6e6
    pmde = np.array([(e["pmdec_mas_yr"] or 0.0) for e in out]) / 3.6e6
    beta100 = ecliptic_lat(np.array([e["ra"] for e in out]) + 100 * pmra / np.maximum(cosd, 1e-6),
                           np.array([e["dec"] for e in out]) + 100 * pmde)
    for e, b, b100 in zip(out, beta, beta100):
        e["ecl_lat_deg"] = round(float(b), 5)
        e["ecl_lat_deg_per_century"] = round(float(b100 - b), 5)
        e.update(etz_flags(float(b)))
        e["etz_years"] = etz_window(float(b), float(b100 - b))

    names = fetch_bl_list(offline)
    nearest = sorted((e for e in out if e["distance_ly"]), key=lambda e: e["distance_ly"])[:100]
    near_set = {id(e) for e in nearest}
    for e in out:
        cands = token_candidates(e["star"], e.pop("_ids"))
        e["bl_token"] = resolve_token(cands, names) if names is not None else None
        e["bl_token_tried"] = cands
        want = check_all or e["near_etz"] or id(e) in near_set or \
            (e["bl_token"] and (e["bl_token"] + "_OFF") in names)
        if e["bl_token"] and want:
            e["bl"] = bl_check(e["bl_token"], names, offline)
            if (e["near_etz"] or id(e) in near_set or check_all) and e["bl"].get("fine_scans"):
                ab = abacad_check(e["bl_token"], offline)
                e["bl"]["abacad_checked"] = ab is not None
                e["bl"]["abacad"] = ab or []
        else:
            e["bl"] = dict(checked=False,
                           why=("no exact archive name" if not e["bl_token"] and names is not None
                                else "not in the checked subset"))
        e["score"], e["score_parts"] = score(e)
        e["our_search"] = OUR_SEARCHES.get(e["bl_token"] or e["star"])

    n_etz = sum(e["etz"] for e in out)
    meta = dict(
        catalogue=CATALOGUE, catalogue_ref=CATALOGUE_REF, n_stars=len(out),
        astrometry="SIMBAD TAP basic (ICRS, J2000 epoch); distance from SIMBAD parallax, "
                   "else Isaacson Dist",
        etz_deg=round(ETZ_DEG, 4), near_etz_deg=NEAR_ETZ_DEG,
        ecliptic_frame="astropy BarycentricTrueEcliptic, equinox J2000",
        score_formula="40*etz + 10*near_etz_only + 25*bl_cadence + 20*max(0,1-d_ly/100) - 30*searched(E1+)",
        bl_checked_subset="near-ETZ stars + nearest 100 + every star whose <token>_OFF is an archive "
                          "target name; fine-res product 0000; bl_archive.cadences rules",
        counts=dict(etz=n_etz,
                    etz_borderline_0264_02665=sum(ETZ_DEG <= abs(e["ecl_lat_deg"]) < 0.2665 for e in out), near_etz=sum(e["near_etz"] for e in out),
                    unique_simbad_objects=len({e["simbad_oid"] for e in out if e["simbad_oid"]}),
                    repeated_star_names=sum(c > 1 for c in __import__("collections").Counter(
                        e["star"] for e in out).values()),
                    rows_sharing_a_simbad_object=sum(bool(e["same_simbad_object_as"]) for e in out),
                    simbad_resolved=sum(e["simbad_main_id"] is not None for e in out),
                    bl_exact_name=sum(e["bl_token"] is not None for e in out),
                    bl_checked=sum(e["bl"].get("checked", False) for e in out),
                    bl_with_cadence=sum(bool(e["bl"].get("cadences")) for e in out),
                    etz_with_cadence=sum(bool(e["bl"].get("cadences")) and e["etz"] for e in out),
                    abacad_checked=sum(bool(e["bl"].get("abacad_checked")) for e in out),
                    with_abacad=sum(bool(e["bl"].get("abacad")) for e in out),
                    etz_with_abacad=sum(bool(e["bl"].get("abacad")) and e["etz"] for e in out),
                    near_etz_with_abacad=sum(bool(e["bl"].get("abacad")) and e["near_etz"] for e in out),
                    usable_any=sum(bool(e["bl"].get("abacad") or e["bl"].get("cadences")) for e in out),
                    near_etz_with_cadence=sum(bool(e["bl"].get("cadences")) and e["near_etz"]
                                              for e in out)),
        checks=crosscheck(out),
    )
    out.sort(key=lambda e: (-e["score"], e["distance_ly"] or 1e9))
    # one star per line: small, and a diff shows exactly which stars changed
    body = ",\n".join(json.dumps(e, separators=(",", ":")) for e in out)
    OUT.write_text('{"meta": ' + json.dumps(meta, indent=1) + ',\n"stars": [\n' + body + "\n]}\n")
    print(json.dumps(meta["counts"], indent=1))
    print(f"wrote {OUT.name}")
    return 0


# ----------------------------------------------------------------- cross-check
REFERENCE = [   # literature values, for a sanity check of the pipeline
    dict(star="GJ699", what="Barnard's Star", d_ly=5.96, rv=-110.5),
    dict(star="HIP8102", what="tau Ceti", d_ly=11.9, rv=-16.7),
    dict(star="HIP16537", what="eps Eridani", d_ly=10.5, rv=16.4),
]


def _find(stars, name):
    n = name.replace(" ", "").lower()
    for e in stars:
        keys = [e["star"], e.get("simbad_main_id") or "", e.get("hip") or "", e.get("common") or "",
                e.get("bl_token") or ""] + (e.get("gj") or [])
        if any(k.replace(" ", "").lower() == n for k in keys):
            return e
    for e in stars:
        if n in (e.get("common") or "").replace(" ", "").lower() or \
                n in (e.get("simbad_main_id") or "").replace(" ", "").lower():
            return e
    return None


def crosscheck(stars):
    res = []
    t = _find(stars, "SO0253")      # Kaltenegger & Faherty 2021: Teegarden's Star enters the ETZ in 2044
    if t:
        y = (t.get("etz_years") or [None])[0]
        res.append(dict(star="SO0253", what="Teegarden's Star ETZ entry year (K&F 2021: 2044)",
                        found=True, got=y, ok=bool(y and abs(y - 2044) <= 3)))
    for ref in REFERENCE:
        e = _find(stars, ref["star"])
        if not e:
            res.append(dict(ref, found=False))
            continue
        res.append(dict(ref, found=True, isaacson_name=e["star"], got_d_ly=e["distance_ly"],
                        got_rv=e["rv_km_s"],
                        ok=bool(e["distance_ly"] and abs(e["distance_ly"] - ref["d_ly"]) < 0.1
                                and e["rv_km_s"] is not None and abs(e["rv_km_s"] - ref["rv"]) < 1.0)))
    return res


# ----------------------------------------------------------------- CLI views
def load():
    if not OUT.exists():
        raise SystemExit("targets.json missing: run  python targets.py build")
    return json.loads(OUT.read_text())


def line(e):
    cad = (e["bl"].get("cadences") or []) + (e["bl"].get("abacad") or [])
    bands = ",".join(sorted({c["band"] for c in cad})) or ("-" if e["bl"].get("checked") else "?")
    label = e.get("common") or (e.get("simbad_main_id") or "").replace("NAME ", "")
    return (f"{e['score']:6.1f}  {e['star']:<10} {label:<18.18} "
            f"{(e['distance_ly'] or 0):7.2f} ly  beta {e['ecl_lat_deg']:+8.3f}"
            f"{' ETZ' if e['etz'] else (' near' if e['near_etz'] else '     ')}  "
            f"BL {str(e['bl_token'] or '-'):<10} cadence {bands}")


def show(name):
    d = load()
    e = _find(d["stars"], name)
    if not e:
        print(f"{name}: not in {CATALOGUE}")
        return 1
    print(json.dumps(e, indent=1))
    return 0


def top(n):
    for e in load()["stars"][:n]:
        print(line(e))
    return 0


def etz():
    d = load()
    m = d["meta"]
    rows = [e for e in d["stars"] if e["near_etz"]]
    rows.sort(key=lambda e: abs(e["ecl_lat_deg"]))
    print(f"{m['counts']['etz']} of {m['n_stars']} Isaacson stars in the ETZ (|beta| < {m['etz_deg']} deg); "
          f"{m['counts']['near_etz']} within {m['near_etz_deg']} deg")
    for e in rows:
        print(line(e))
        for c in (e["bl"].get("cadences") or []) + (e["bl"].get("abacad") or []):
            offs = f", OFFs {'/'.join(c['off_targets'])}" if c.get("off_targets") else ", OFFs <name>_OFF"
            print(f"        {c['telescope']} {c['project']} MJD {c['mjd_day']} {c['band']}-band "
                  f"fc {c['fc_mhz']} MHz, {c['gb']} GB{offs}")
    return 0


def check():
    for c in load()["meta"]["checks"]:
        print(c)
    return 0


# ----------------------------------------------------------------- selftest
FIXTURE_TSV = """#fixture: 4 rows in the VizieR tsv layout of J/PASP/129/E4501/table1
_RAJ2000\t_DEJ2000\trecno\tStar\tRAJ2000\tDEJ2000\tEp\tVmag\tSpType\tDist\tpmRA\tpmDE\tSimbadName
deg\tdeg\t \t \t"h:m:s"\t"d:m:s"\tyr\tmag\t \tpc\tarcsec/yr\tarcsec/yr\t
----------\t----------\t--------\t---------\t----------\t---------\t----\t-----\t-------\t-----\t-----\t-----\t------------------------
269.452083\t +4.693611\t     1\tGJ699    \t17 57 48.5\t+04 41 37\t2000\t 9.51\tM4V    \t 1.83\t-0.80\t10.33\tGJ699
000.004167\t-19.498611\t     2\tHIP2     \t00 00 01.0\t-19 29 55\t2000\t 9.27\tK3V    \t45.60\t 0.18\t 0.00\tHIP2
026.017083\t-15.937500\t     3\tHIP8102  \t01 44 04.1\t-15 56 15\t2000\t 3.50\tG8V    \t 3.65\t-1.72\t 0.85\tHIP8102
053.232917\t -9.458333\t     4\tGJ15A    \t03 32 55.9\t-09 27 30\t2000\t 3.73\tK2V    \t 3.22\t-0.98\t 0.02\tGJ15A
"""


def selftest():
    checks = []
    rows = parse_isaacson(FIXTURE_TSV)
    checks.append(("fixture parses 4 rows with distance and pm",
                   len(rows) == 4 and rows[0]["star"] == "GJ699" and abs(rows[0]["dist_pc"] - 1.83) < 1e-9
                   and abs(rows[0]["pmde"] - 10.33) < 1e-9))
    checks.append(("SIMBAD id forms", [simbad_id(x) for x in ("HIP2", "GJ15A", "GJ54.1", "GJ699", "LP944-020")]
                   == ["HIP 2", "GJ 15 A", "GJ 54.1", "GJ 699", "LP 944-020"]))
    # ETZ geometry: points built IN the ecliptic frame, transformed to ICRS and back
    on = ecliptic_to_icrs(123.0, 0.0)
    pole = ecliptic_to_icrs(0.0, 90.0)
    edge_in = ecliptic_to_icrs(200.0, 0.25)
    edge_out = ecliptic_to_icrs(200.0, 0.28)
    near = ecliptic_to_icrs(10.0, -1.5)
    b = ecliptic_lat([p[0] for p in (on, pole, edge_in, edge_out, near)],
                     [p[1] for p in (on, pole, edge_in, edge_out, near)])
    checks.append(("star on the ecliptic -> ETZ", etz_flags(b[0])["etz"] and abs(b[0]) < 1e-6))
    checks.append(("ecliptic pole -> not ETZ, not near",
                   not etz_flags(b[1])["etz"] and not etz_flags(b[1])["near_etz"] and abs(b[1] - 90) < 1e-6))
    checks.append(("ecliptic north pole is RA 18h, Dec +66.56", abs(pole[0] - 270) < 0.01 and abs(pole[1] - 66.56) < 0.01))
    checks.append(("beta 0.25 in, 0.28 out (ETZ = 0.264 deg)",
                   etz_flags(b[2])["etz"] and not etz_flags(b[3])["etz"] and ETZ_DEG == 0.264))
    checks.append(("ETZ window: beta 0.32 moving -0.13/century enters ~2044",
                   etz_window(0.321, -0.129)[0] == 2044 and etz_window(5.0, 0.001) is None))
    checks.append(("beta -1.5 near-ETZ only", etz_flags(b[4])["near_etz"] and not etz_flags(b[4])["etz"]))
    # Barnard's Star: beta ~ +28 deg (well off the ecliptic) - control that must NOT fire
    bb = ecliptic_lat([rows[0]["ra"]], [rows[0]["dec"]])[0]
    checks.append(("Barnard's Star beta +27..+29 deg, not ETZ", 27 < bb < 29 and not etz_flags(bb)["near_etz"]))
    # SIMBAD lookup: by identifier first, then the position fallback, else Isaacson
    s0, s1 = rows[0], dict(rows[1], star="SO0253")
    byq = {"GJ 699": dict(main_id="NAME Barnard's star"),
           "SO0253|SO0253|" + f"{s1['ra']:.4f}": dict(main_id="Teegarden", by_position_arcsec=1.0)}
    r0, src0 = simbad_row(byq, s0)
    r1, src1 = simbad_row(byq, s1)
    r2, src2 = simbad_row(byq, rows[2])
    checks.append(("SIMBAD row: identifier, then position, else catalogue",
                   r0["main_id"].endswith("star") and "identifier" in src0 and r1["main_id"] == "Teegarden"
                   and "nearest" in src1 and r2 is None and src2.startswith("Isaacson")))
    # Prefix trap: archive has HIP28614/HIP2_x but NOT HIP2 -> no token
    names = {"HIP28614", "HIP24284", "HIP2_S", "GJ699", "GJ699_OFF", "HIP8102"}
    checks.append(("prefix trap: HIP2 gets no token", resolve_token(["HIP2"], names) is None))
    checks.append(("exact token found", resolve_token(["GJ699", "HIP87937"], names) == "GJ699"))
    checks.append(("alternate HIP form used when Isaacson name absent",
                   resolve_token(token_candidates("GJ71", ["HIP 8102", "GJ 71"]), names) == "HIP8102"))
    # bl_check on a canned prefix reply (offline, via the cache)
    import bl_archive
    m = 5 / 1440

    def row(t, i, dec):
        return dict(target=t, mjd=57000 + i * m, ra=10.0, decl=dec, center_freq=1501.5, telescope="GBT",
                    size=2e9, url=f"https://x/pipeline/P1/h/spliced_guppi_57000_1_{t}_{i:04d}.gpuspec.0000.h5")
    good = [row("S", i, 5.0) if i % 2 == 0 else row("S_OFF", i, 7.0) for i in range(6)]
    trap = [row("S2", i + 10, 5.0) if i % 2 == 0 else row("S2_OFF", i + 10, 7.0) for i in range(6)]
    g, _ = bl_archive.cadences("S", trap)
    checks.append(("prefix-matched S2 cadence is not S's", len(g) == 0))
    g, _ = bl_archive.cadences("S", good + trap)
    checks.append(("S's own cadence still found beside S2's", len(g) == 1 and band_of(g[0]["fc"]) == "L"))
    # ABACAD with other stars as OFFs (time-window rows, any name)
    def wrow(t, i, ra, dec, mjd=None):
        return dict(target=t, mjd=mjd if mjd is not None else 57000 + i * m, ra=ra, decl=dec,
                    center_freq=1501.5, telescope="GBT", size=2e9,
                    url=f"https://x/pipeline/P1/h/spliced_guppi_57000_1_{t}_{i:04d}.gpuspec.0000.h5")
    sr = [wrow("S", i, 10.0, 5.0) for i in (0, 2, 4)]
    win = sr + [wrow("B1", 1, 12.0, 5.0), wrow("B2", 3, 9.0, 3.0), wrow("B3", 5, 11.0, 8.0)]
    g = abacad("S", sr, win)
    checks.append(("ABACAD with three other-star OFFs found",
                   len(g) == 1 and g[0]["off_targets"] == ["B1", "B2", "B3"]))
    win_off = sr + [wrow("S_OFF", i, 10.0, 7.0) for i in (1, 3, 5)]
    checks.append(("ABACAD finder also accepts <name>_OFF OFFs", len(abacad("S", sr, win_off)) == 1))
    win_close = sr + [wrow("B1", 1, 10.05, 5.0), wrow("B2", 3, 9.0, 3.0), wrow("B3", 5, 11.0, 8.0)]
    checks.append(("OFF 0.05 deg from the star rejected", len(abacad("S", sr, win_close)) == 0))
    win_pre = [wrow("S2", i, 10.0, 5.0) for i in (0, 2, 4)] + win[3:]
    checks.append(("prefix-twin S2 ON scans never count as S", len(abacad("S", win_pre, win_pre)) == 0))
    win_gap = sr[:2] + [wrow("S", 4, 10.0, 5.0, mjd=57000 + 4 * m + 0.05)] + win[3:5] + \
        [wrow("B3", 5, 11.0, 8.0, mjd=57000 + 5 * m + 0.05)]
    checks.append(("72-min gap rejected", len(abacad("S", win_gap[:3], win_gap)) == 0))
    # Score ordering
    base = dict(star="X", bl_token=None)
    a = dict(base, etz=True, near_etz=True, distance_ly=50, bl=dict(cadences=[1]))
    bq = dict(base, etz=True, near_etz=True, distance_ly=50, bl=dict(cadences=[]))
    c = dict(base, etz=False, near_etz=True, distance_ly=50, bl=dict(cadences=[1]))
    d = dict(base, etz=False, near_etz=False, distance_ly=5, bl=dict(cadences=[1]))
    e = dict(base, star="GJ699", bl_token="GJ699", etz=False, near_etz=False, distance_ly=5.96,
             bl=dict(cadences=[1]))
    f = dict(base, star="GJ144", bl_token="GJ144", etz=False, near_etz=False, distance_ly=5.96,
             bl=dict(cadences=[1]))
    sc = [score(x)[0] for x in (a, bq, c, d, e, f)]
    checks.append(("score: ETZ+cadence > ETZ > near+cadence > close non-ETZ", sc[0] > sc[1] > sc[2] > sc[3]))
    checks.append(("score: E2-searched star demoted; VOID search not counted", sc[4] < sc[3] and sc[5] > sc[4]))
    checks.append(("score of far star never negative closeness", score(dict(base, etz=False, near_etz=False,
                   distance_ly=500, bl={}))[0] == 0))
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(x[1] for x in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sp = p.add_subparsers(dest="cmd", required=True)
    b = sp.add_parser("build")
    b.add_argument("--offline", action="store_true", help="use data/targets_raw/ only")
    b.add_argument("--check-all", action="store_true", help="query the archive for every exact name")
    s = sp.add_parser("show")
    s.add_argument("name")
    t = sp.add_parser("top")
    t.add_argument("n", type=int, nargs="?", default=20)
    sp.add_parser("etz")
    sp.add_parser("check")
    sp.add_parser("selftest")
    a = p.parse_args()
    if a.cmd == "build":
        return build(a.offline, a.check_all)
    if a.cmd == "show":
        return show(a.name)
    if a.cmd == "top":
        return top(a.n)
    if a.cmd == "etz":
        return etz()
    if a.cmd == "check":
        return check()
    return selftest()


if __name__ == "__main__":
    sys.exit(main())
