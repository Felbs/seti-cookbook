#!/usr/bin/env python3
"""bl_archive.py - query the Breakthrough Listen open-data archive EXACTLY.

The archive's ?target= parameter matches by PREFIX: asking for "HIP2" returns
HIP28614, HIP24284, ... and asking for "GJ699" also returns "GJ699_OFF". A
pipeline that trusts the reply can download a different star, or build an
"ON/OFF cadence" out of two sessions years apart on different receivers.
(Found by the 2026-11-06 peer review; seti-cookbook law L12.)

This module is the one place that turns an archive reply into scans:

  ON  scan  = row["target"] == the name asked for, exactly
  OFF scan  = row["target"] == name + "_OFF" exactly, pointed >= 0.2 deg away
  cadence   = 6 consecutive scans of ONE session, ONE band (centre frequency),
              ONE data product, alternating ON OFF ON OFF ON OFF, no gap > 20 min

    python bl_archive.py GJ699            # list every cadence the archive has
    python bl_archive.py HIP2             # shows the prefix trap: 0 exact rows
    python bl_archive.py selftest         # offline checks of the rules above

Credit: Breakthrough Listen open data (Lebofsky et al. 2019 PASP 131:124505).
"""
import json
import math
import re
import sys
import urllib.parse
import urllib.request

API = "http://seti.berkeley.edu/opendata/api/query-files"
MAX_GAP_MIN = 20.0       # longest allowed pause between scans of one cadence
MIN_OFF_DEG = 0.2        # an OFF must point away from the target by at least this


def query(target, limit=2000, telescope=None):
    """Raw archive reply (prefix-matched - do not use directly)."""
    p = {"target": target, "file-types": "HDF5", "limit": limit}
    if telescope:
        p["telescopes"] = telescope
    with urllib.request.urlopen(API + "?" + urllib.parse.urlencode(p), timeout=90) as r:
        return json.loads(r.read()).get("data", [])


def product(name):
    """BL data product: 0000 fine frequency (~2.8 Hz), 0001 fine time, 0002 mid."""
    m = re.search(r"\.(000[012])\.h5$", name)
    return m.group(1) if m else "?"


def session(url):
    """(project, MJD day) from the URL, e.g. ('AGBT16A_999_81', 57457)."""
    parts = url.split("/")
    proj = parts[4] if len(parts) > 4 else "?"
    m = re.search(r"guppi_(\d{5})_", url)
    return proj, int(m.group(1)) if m else -1


def scan_no(name):
    m = re.search(r"_(\d{4})\.gpuspec", name)
    return int(m.group(1)) if m else -1


def sep_deg(ra1, de1, ra2, de2):
    r = math.radians
    c = (math.sin(r(de1)) * math.sin(r(de2))
         + math.cos(r(de1)) * math.cos(r(de2)) * math.cos(r(ra1 - ra2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def exact_rows(target, rows):
    """Split a prefix-matched reply into (on, off, dropped) by EXACT name."""
    on = [x for x in rows if x.get("target") == target]
    off = [x for x in rows if x.get("target") == target + "_OFF"]
    dropped = [x for x in rows if x.get("target") not in (target, target + "_OFF")]
    return on, off, dropped


def scans(target, rows, prod="0000"):
    """Exact-name scans of one product, tagged ON/OFF, one row per pointing."""
    on, off, _ = exact_rows(target, rows)
    out = {}
    for kind, group in (("ON", on), ("OFF", off)):
        for x in group:
            name = x["url"].split("/")[-1]
            if product(name) != prod:
                continue
            proj, day = session(x["url"])
            key = (proj, day, round(float(x["center_freq"]), 1), scan_no(name))
            out.setdefault(key, dict(kind=kind, url=x["url"], name=name,
                                     size=x.get("size", 0), mjd=float(x["mjd"]),
                                     ra=float(x["ra"]), dec=float(x["decl"]),
                                     project=proj, day=day,
                                     fc=round(float(x["center_freq"]), 1),
                                     md5=x.get("md5sum")))
    return sorted(out.values(), key=lambda s: s["mjd"])


def cadences(target, rows, prod="0000"):
    """Every valid ABACAD cadence in the reply, plus a list of rejected groups
    with the reason, so nothing is dropped silently."""
    groups = {}
    for s in scans(target, rows, prod):
        groups.setdefault((s["project"], s["day"], s["fc"]), []).append(s)
    good, rejected = [], []
    for key, g in sorted(groups.items()):
        ons = [s for s in g if s["kind"] == "ON"]
        if not ons:
            rejected.append((key, "no ON scans"))
            continue
        ra0 = sorted(s["ra"] for s in ons)[len(ons) // 2]
        de0 = sorted(s["dec"] for s in ons)[len(ons) // 2]
        found = False
        for i in range(len(g) - 5):
            w = g[i:i + 6]
            pattern = [s["kind"] for s in w]
            if pattern != ["ON", "OFF"] * 3:
                continue
            gaps = [(b["mjd"] - a["mjd"]) * 1440 for a, b in zip(w, w[1:])]
            if max(gaps) > MAX_GAP_MIN:
                continue
            offs = [sep_deg(ra0, de0, s["ra"], s["dec"]) for s in w if s["kind"] == "OFF"]
            if min(offs) < MIN_OFF_DEG:
                continue
            good.append(dict(project=key[0], day=key[1], fc=key[2], scans=w,
                             max_gap_min=max(gaps), off_sep_deg=min(offs)))
            found = True
            break
        if not found:
            pat = "".join("A" if s["kind"] == "ON" else "-" for s in g)
            rejected.append((key, f"no ON/OFF x3 window within {MAX_GAP_MIN:.0f} min "
                                  f"(scan pattern {pat})"))
    return good, rejected


def report(target, prod="0000"):
    rows = query(target)
    on, off, dropped = exact_rows(target, rows)
    print(f"{target}: archive returned {len(rows)} rows; exact ON {len(on)}, "
          f"exact OFF {len(off)}, dropped {len(dropped)} other targets (prefix match)")
    if dropped:
        names = sorted({x['target'] for x in dropped})
        print(f"  dropped names: {', '.join(names[:8])}{' ...' if len(names) > 8 else ''}")
    good, rejected = cadences(target, rows, prod)
    for c in good:
        gb = sum(s["size"] for s in c["scans"]) / 1e9
        print(f"  CADENCE {c['project']} MJD {c['day']} band {c['fc']} MHz product {prod}: "
              f"{gb:.1f} GB, max gap {c['max_gap_min']:.1f} min, OFF {c['off_sep_deg']:.2f} deg away")
        for s in c["scans"]:
            print(f"     {s['kind']:<3} {s['mjd']:.5f}  {s['name']}")
    for key, why in rejected:
        print(f"  rejected {key[0]} MJD {key[1]} band {key[2]}: {why}")
    return good


def selftest():
    def row(t, mjd, ra=10.0, dec=5.0, fc=1500.0, proj="AGBT16A_1", day=57000, n=1, p="0000"):
        url = (f"https://x/pipeline/{proj}/h/spliced_guppi_{day}_1_{t}_{n:04d}"
               f".gpuspec.{p}.h5")
        return dict(target=t, mjd=mjd, ra=ra, decl=dec, center_freq=fc, url=url, size=1)
    m = 5 / 1440
    base = [row("S", 57000 + i * m, dec=5.0 if i % 2 == 0 else 7.0, n=i,
                ) if i % 2 == 0 else row("S_OFF", 57000 + i * m, dec=7.0, n=i)
            for i in range(6)]
    checks = []
    g, _ = cadences("S", base)
    checks.append(("textbook ABACAD accepted", len(g) == 1))
    prefix = base + [row("S2", 57000.1, n=9), row("S2_OFF", 57000.11, n=10)]
    g, _ = cadences("S", prefix)
    checks.append(("prefix-matched other star ignored", len(g) == 1 and
                   all(s["name"].find("_S2") < 0 for s in g[0]["scans"])))
    split = base[:3] + [dict(r, url=r["url"].replace("AGBT16A_1", "AGBT17B_9")) for r in base[3:]]
    g, _ = cadences("S", split)
    checks.append(("two sessions never merged", len(g) == 0))
    band = base[:3] + [dict(r, center_freq=2250.0) for r in base[3:]]
    g, _ = cadences("S", band)
    checks.append(("two bands never merged", len(g) == 0))
    near = [dict(r, decl=5.05) if r["target"] == "S_OFF" else r for r in base]
    g, _ = cadences("S", near)
    checks.append(("OFF too close to target rejected", len(g) == 0))
    slow = [dict(r, mjd=r["mjd"] + (0.1 if i >= 3 else 0)) for i, r in enumerate(base)]
    g, _ = cadences("S", slow)
    checks.append(("2-hour gap rejected", len(g) == 0))
    res = base + [dict(r, url=r["url"].replace(".0000.h5", ".0002.h5")) for r in base]
    g, _ = cadences("S", res)
    checks.append(("one product only, no duplicates", len(g) == 1 and
                   all(s["name"].endswith(".0000.h5") for s in g[0]["scans"])))
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    a = sys.argv[1:] or ["GJ699"]
    if a[0] == "selftest":
        sys.exit(selftest())
    for t in a:
        report(t)
