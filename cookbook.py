#!/usr/bin/env python3
"""cookbook.py - run, score and compare seti-cookbook RECIPES.

A recipe is one small python file in recipes/ that says "here is a new way a
civilization might be detectable, and here is how to look for it" (the contract
is recipe_api.py, the tutorial is RECIPES.md). This is the CLI that runs them.

  python cookbook.py list                       # every discovered recipe
  python cookbook.py selftest                   # run every recipe's own selftest
  python cookbook.py run all synth:drift        # all recipes on a synthetic scene
  python cookbook.py run all data/star_GJ699.h5 --f-start 1420.2 --f-stop 1420.7
  python cookbook.py run narrowband_drift data/Voyager1.single_coarse.fine_res.h5 \
        --f-start 8420.18 --f-stop 8420.25
  python cookbook.py bench                      # score everyone, write LEADERBOARD.md
  python cookbook.py cadence narrowband_drift ON=a.h5 OFF=b.h5 ON=c.h5 OFF=d.h5

Every candidate printed gets a VERDICT from recipe_api.explain() - known-RFI
band, natural spectral line, zero drift, band edge, or "unexplained - worth a
human". Nothing is deleted silently; that is the repo's honesty rail.
"""
import argparse
import json
import sys
from pathlib import Path

import recipe_api as R
import seti_io

HERE = Path(__file__).resolve().parent


def _open(path, args):
    kw = {}
    if args.fs:
        kw["fs"] = args.fs
    for a, k in (("f_start", "f_start"), ("f_stop", "f_stop"),
                 ("t_start", "t_start"), ("t_stop", "t_stop")):
        v = getattr(args, a, None)
        if v is not None and not str(path).startswith("synth:") \
                and Path(str(path)).suffix.lower() in (".h5", ".hdf5", ".fil"):
            kw[k] = v
    return seti_io.open_any(path, **kw)


def cmd_list(args):
    recs, broken = R.discover()
    print(f"{len(recs)} recipe(s) in {R.RECIPE_DIR}:\n")
    for r in recs:
        i = r.info()
        print(f"  {i['name']:20s} v{i['version']:5s} [{','.join(i['tags'])}]")
        print(f"      {i['description']}")
        print(f"      by {i['author']}  ({i['file']}"
              f"{', selftest' if i['has_selftest'] else ', NO selftest'})")
    for name, why in broken:
        print(f"  !! {name}: {why}")
    return 0


def cmd_selftest(args):
    recs, broken = R.discover()
    bad = list(broken)
    fails = []
    for r in recs:
        if not hasattr(r.module, "selftest"):
            print(f"  {r.name}: no selftest (allowed, but please add one)")
            continue
        print("-" * 70)
        try:
            rc = r.module.selftest()
        except Exception as e:
            rc = 1
            print(f"  {r.name} raised {type(e).__name__}: {e}")
        if rc:
            fails.append(r.name)
    print("=" * 70)
    print(f"{len(recs) - len(fails)}/{len(recs)} recipes green"
          + (f"; FAILING: {fails}" if fails else "")
          + (f"; BROKEN FILES: {bad}" if bad else ""))
    return 1 if fails or bad else 0


def cmd_run(args):
    spec = _open(args.data, args)
    print(f"data: {spec}")
    recs, broken = R.discover()
    for n, why in broken:
        print(f"  !! broken recipe {n}: {why}")
    if args.recipe != "all":
        recs = [r for r in recs if r.name == args.recipe or r.path.stem == args.recipe]
        if not recs:
            print(f"no recipe named {args.recipe!r}")
            return 1
    allc = []
    for r in recs:
        try:
            cands, dt = r.run(spec, json.loads(args.params) if args.params else None)
        except Exception as e:
            print(f"\n{r.name}: ERROR {type(e).__name__}: {e}")
            continue
        R.explain(cands, spec)
        allc += cands
        print(f"\n{r.name}  ({dt:.2f}s, {len(cands)} candidate(s))")
        for c in sorted(cands, key=lambda c: -c.score)[:args.top]:
            print(f"   {c.freq_mhz:14.6f} MHz  score {c.score:9.1f}  "
                  f"drift {c.drift_hz_s:+7.3f} Hz/s  [{c.kind}]")
            print(f"       {c.label}")
            print(f"       verdict: {c.verdict}")
    if args.json:
        Path(args.json).write_text(json.dumps(
            {"data": spec.summary(), "candidates": [c.to_dict() for c in allc]},
            indent=1), encoding="utf-8")
        print(f"\nwrote {args.json} ({len(allc)} candidates)")
    unexplained = [c for c in allc if c.verdict.startswith("unexplained")]
    print(f"\ntotal {len(allc)} candidate(s); {len(unexplained)} unexplained "
          f"after RFI / natural-line / zero-drift vetoes")
    return 0


def cmd_bench(args):
    print("=" * 70)
    print("seti-cookbook recipe benchmark - same scenes, same truth, for everyone")
    print("=" * 70)
    rows, scenes = R.bench()
    p = R.write_leaderboard(rows, scenes)
    print("=" * 70)
    print(f"wrote {p} and {p.with_suffix('.json')}")
    winner = rows[0] if rows else None
    if winner:
        print(f"top clean recipe: {winner['recipe']} "
              f"({winner['n_caught']} scenes, {winner['n_fa']} false alarms)")
    uncaught = set(s for s in scenes if not s.startswith("NULL")) - set(
        s for r in rows for s in r["caught"])
    if uncaught:
        print(f"NOBODY catches: {sorted(uncaught)}  <- open bounty, "
              "write that recipe (see RECIPES.md)")
    return 0


def cmd_cadence(args):
    """ON/OFF cadence verification - the standard SETI check, and the thing that
    killed BLC-1 (Sheikh et al. 2021)."""
    r = R.get(args.recipe)
    seq = []
    for item in args.pointings:
        role, _, path = item.partition("=")
        spec = _open(path, args)
        cands, _ = r.run(spec)
        R.explain(cands, spec)
        print(f"  {role.upper():4s} {Path(path).name}: {len(cands)} candidate(s)")
        seq.append((role, cands))
    survivors = R.cadence_verify(seq, tol_hz=args.tol_hz)
    print(f"\n{len(survivors)} candidate(s) present in the ON scans and absent "
          f"from every OFF scan:")
    for c in survivors:
        print(f"   {c.freq_mhz:14.6f} MHz score {c.score:8.1f}  {c.verdict}")
    if not survivors:
        print("   (none - which is what every honest SETI observation has "
              "reported so far)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")

    def common(p):
        p.add_argument("--fs", type=float, help="sample rate for raw IQ input")
        p.add_argument("--f-start", type=float, help="MHz")
        p.add_argument("--f-stop", type=float, help="MHz")
        p.add_argument("--t-start", type=float, help="seconds")
        p.add_argument("--t-stop", type=float, help="seconds")

    sub.add_parser("list")
    sub.add_parser("selftest")
    p = sub.add_parser("run")
    p.add_argument("recipe")
    p.add_argument("data")
    p.add_argument("--params", help="JSON dict of recipe parameters")
    p.add_argument("--json", help="write candidates to this JSON file")
    p.add_argument("--top", type=int, default=6)
    common(p)
    sub.add_parser("bench")
    p = sub.add_parser("cadence")
    p.add_argument("recipe")
    p.add_argument("pointings", nargs="+", help="ON=file OFF=file ...")
    p.add_argument("--tol-hz", type=float, default=200.0)
    common(p)

    args = ap.parse_args()
    if not args.cmd:
        ap.print_help()
        return 0
    return {"list": cmd_list, "selftest": cmd_selftest, "run": cmd_run,
            "bench": cmd_bench, "cadence": cmd_cadence}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
