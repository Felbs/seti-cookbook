#!/usr/bin/env python3
"""search_ledger.py — what have we actually searched, and how much is left?

A null result is only worth something if you can say what you covered. Without a
ledger you cannot quote a limit, you cannot avoid re-searching the same ground,
and "we found nothing" means nothing. This records every run so the search has a
denominator.

    python search_ledger.py record --target GJ699 --path data/x.h5 --recipe all
    python search_ledger.py coverage          # the headline percentages
    python search_ledger.py list              # what has been searched

WHAT "PERCENT SEARCHED" HONESTLY MEANS
--------------------------------------
There is no single true denominator, so this reports several and refuses to
collapse them into one number:

  * BYTES of the public archive we have actually opened. Breakthrough Listen's
    open-data release is the big one — of order 1 PB and growing (Lebofsky et al.
    2019 PASP 131:124505; Price et al. 2020 AJ 159:86). This is the number that
    is honest about how little any hobbyist has touched. Set the denominator
    yourself with SETI_ARCHIVE_BYTES if you have a better figure; the default is
    an ORDER-OF-MAGNITUDE estimate and is labelled as such everywhere.
  * TARGETS: distinct stars we have looked at, against the ~1,700 in the BL
    primary nearby-star sample.
  * The COVERAGE VOLUME that actually matters scientifically, which is not bytes
    at all: (targets) x (bandwidth) x (time) x (drift range) x (sensitivity). Two
    searches over the same file with different detectors are NOT duplicates —
    they cover different signal space. Bytes-searched is an upper bound on effort
    and a lower bound on insight.

So: treat the byte percentage as "how much of the haystack have we lifted", and
the per-recipe records as "which needles could we have seen". Both are in here.
"""
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEDGER = HERE / "search_ledger.json"

# Order-of-magnitude ONLY. BL's open data release passed ~1 PB around 2019-2020
# and has grown since; there is no single authoritative live figure, so this is
# deliberately a round number you can override rather than a false precision.
DEFAULT_ARCHIVE_BYTES = float(os.environ.get("SETI_ARCHIVE_BYTES", 1e15))
BL_PRIMARY_TARGETS = int(os.environ.get("SETI_ARCHIVE_TARGETS", 1702))


def _load():
    if not LEDGER.exists():
        return {"runs": [], "created": datetime.now(timezone.utc).isoformat()}
    try:
        return json.loads(LEDGER.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"runs": [], "created": datetime.now(timezone.utc).isoformat()}


def _save(led):
    LEDGER.write_text(json.dumps(led, indent=2), encoding="utf-8")


def record(path, recipe, target=None, n_candidates=0, n_unexplained=0,
           f_lo_mhz=None, f_hi_mhz=None, duration_s=None, res_hz=None,
           drift_max_hz_s=None, bytes_read=None, note=""):
    """Append one search. Idempotent-ish: the same (path, recipe) is recorded
    again because a re-run with different params IS a different search — the
    coverage maths de-duplicates by file for BYTES and keeps every run for
    signal-space accounting."""
    led = _load()
    p = Path(path)
    led["runs"].append(dict(
        when=datetime.now(timezone.utc).isoformat(),
        path=str(p), file=p.name, target=target, recipe=recipe,
        n_candidates=int(n_candidates), n_unexplained=int(n_unexplained),
        f_lo_mhz=f_lo_mhz, f_hi_mhz=f_hi_mhz,
        bandwidth_mhz=(None if f_lo_mhz is None or f_hi_mhz is None
                       else round(f_hi_mhz - f_lo_mhz, 6)),
        duration_s=duration_s, res_hz=res_hz, drift_max_hz_s=drift_max_hz_s,
        bytes_read=int(bytes_read) if bytes_read
        else (p.stat().st_size if p.exists() else None),
        note=note))
    _save(led)
    return led["runs"][-1]


def coverage():
    """The headline numbers, with every denominator stated rather than implied."""
    led = _load()
    runs = led.get("runs", [])
    by_file = {}
    for r in runs:
        if r.get("bytes_read"):
            by_file[r["file"]] = max(by_file.get(r["file"], 0), r["bytes_read"])
    bytes_searched = sum(by_file.values())
    targets = sorted({r["target"] for r in runs if r.get("target")})
    recipes = sorted({r["recipe"] for r in runs if r.get("recipe")})
    band_mhz = sum(r["bandwidth_mhz"] for r in runs if r.get("bandwidth_mhz"))
    hours = sum((r.get("duration_s") or 0) for r in runs) / 3600.0

    pct_bytes = 100.0 * bytes_searched / DEFAULT_ARCHIVE_BYTES
    pct_targets = 100.0 * len(targets) / BL_PRIMARY_TARGETS
    return {
        "runs": len(runs),
        "distinct_files": len(by_file),
        "bytes_searched": bytes_searched,
        "gb_searched": round(bytes_searched / 1e9, 3),
        "archive_bytes_assumed": DEFAULT_ARCHIVE_BYTES,
        "archive_note": "ORDER-OF-MAGNITUDE estimate of the BL open release "
                        "(~1 PB, Lebofsky+2019 / Price+2020). Override with "
                        "SETI_ARCHIVE_BYTES.",
        "percent_of_archive_bytes": pct_bytes,
        "percent_of_archive_bytes_str": _tiny(pct_bytes),
        "targets_searched": targets,
        "n_targets": len(targets),
        "bl_primary_targets": BL_PRIMARY_TARGETS,
        "percent_of_primary_targets": round(pct_targets, 4),
        "recipes_used": recipes,
        "integrated_bandwidth_mhz": round(band_mhz, 3),
        "integrated_hours": round(hours, 4),
        "signal_space_note": "bytes are effort, not insight: two recipes on the "
                             "same file cover DIFFERENT signal space. Use "
                             "recipes_used x bandwidth x drift range for that.",
    }


def _tiny(pct):
    """Render honestly small percentages without pretending they are zero."""
    if pct <= 0:
        return "0"
    if pct >= 0.01:
        return f"{pct:.4f}%"
    return f"{pct:.3e}%  (about 1 part in {int(round(100.0 / pct)):,})"


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    r = sub.add_parser("record")
    r.add_argument("--path", required=True)
    r.add_argument("--recipe", required=True)
    r.add_argument("--target")
    r.add_argument("--candidates", type=int, default=0)
    r.add_argument("--unexplained", type=int, default=0)
    r.add_argument("--note", default="")
    sub.add_parser("coverage")
    sub.add_parser("list")
    a = ap.parse_args()

    if a.cmd == "record":
        row = record(a.path, a.recipe, a.target, a.candidates, a.unexplained,
                     note=a.note)
        print(json.dumps(row, indent=2))
    elif a.cmd == "list":
        for x in _load().get("runs", []):
            print(f"{x['when'][:19]}  {x.get('target') or '-':<12} "
                  f"{x['recipe']:<20} {x['file'][:44]:<46} "
                  f"cand={x['n_candidates']}")
    else:
        c = coverage()
        print("=" * 66)
        print("  seti-cookbook search coverage — what we have actually looked at")
        print("=" * 66)
        print(f"  runs                  {c['runs']}")
        print(f"  distinct files        {c['distinct_files']}")
        print(f"  data searched         {c['gb_searched']} GB")
        print(f"  of the public archive {c['percent_of_archive_bytes_str']}")
        print(f"      (denominator: {c['archive_bytes_assumed']:.0e} bytes — "
              f"order-of-magnitude)")
        print(f"  targets               {c['n_targets']} of "
              f"{c['bl_primary_targets']} "
              f"({c['percent_of_primary_targets']}%)  {c['targets_searched']}")
        print(f"  recipes used          {', '.join(c['recipes_used']) or '-'}")
        print(f"  integrated bandwidth  {c['integrated_bandwidth_mhz']} MHz")
        print(f"  integrated time       {c['integrated_hours']} h")
        print()
        print("  " + c["signal_space_note"])
        print("=" * 66)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
