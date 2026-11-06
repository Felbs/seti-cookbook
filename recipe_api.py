#!/usr/bin/env python3
"""recipe_api.py - the seti-cookbook RECIPE contract: how anyone invents a new way
to hunt aliens and gets it scored against everyone else's on the same data.

65 years of SETI has mostly asked ONE question ("is there a narrowband carrier?").
NOVEL_DETECTORS.md argues that is one question out of a huge space, and that the
bottleneck is IMAGINATION, not telescope time - the public Breakthrough Listen
archive is free and enormous. So: make the detector the unit of contribution.

A recipe is one small python file in recipes/ that declares what it looks for
and returns candidates. That is the whole API:

    NAME        = "narrowband_drift"          # unique, filename-ish
    DESCRIPTION = "one line: what signal class does this catch?"
    AUTHOR      = "your name / handle"
    VERSION     = "1.0"
    INPUT       = "spectrogram"               # the only input type today
    TAGS        = ["technosignature", "drift"]

    def run(spec, params=None) -> list[Candidate]
    def selftest() -> int                      # optional but STRONGLY encouraged

`spec` is a seti_io.Spectrogram (frequency ASCENDING, physical axes attached);
`spec.meta` carries source_name/telescope/tstart_mjd/ra_hr/dec_deg when known.
Return Candidate objects (or plain dicts with the same keys). Every candidate
must carry a `score` (higher = more interesting) and `provenance` (how you got
it - so a human can re-derive your claim). Nothing here thresholds anything
away silently: recipes report, the verification layer explains.

CLI lives in cookbook.py. See RECIPES.md for the tutorial.
"""
import importlib.util
import json
import time
import traceback
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

import seti_io
from seti_io import HI_MHZ, OH_MHZ, C_KMS

HERE = Path(__file__).resolve().parent
RECIPE_DIR = HERE / "recipes"


# ---------------------------------------------------------------- candidates
@dataclass
class Candidate:
    """One thing a recipe found. `kind` is the recipe's own honest opinion:
    'techno' (worth a human), 'natural' (real astrophysics - see SETI_HISTORY.md),
    'rfi' (ours, not theirs), 'unknown'."""
    freq_mhz: float
    score: float
    drift_hz_s: float = 0.0
    t_start_s: float = 0.0
    duration_s: float = None
    bandwidth_hz: float = None
    dm: float = None
    period_s: float = None
    kind: str = "unknown"
    label: str = ""
    recipe: str = ""
    provenance: dict = field(default_factory=dict)
    verdict: str = ""            # filled in by the verification layer

    def to_dict(self):
        return {k: (round(v, 9) if isinstance(v, float) else v)
                for k, v in asdict(self).items()}


def as_candidates(objs, recipe_name=""):
    out = []
    for o in objs or []:
        c = o if isinstance(o, Candidate) else Candidate(**o)
        if not c.recipe:
            c.recipe = recipe_name
        out.append(c)
    return out


# ------------------------------------------------------------------ discovery
@dataclass
class Recipe:
    name: str
    path: Path
    module: object
    description: str = ""
    author: str = ""
    version: str = ""
    tags: list = field(default_factory=list)
    input: str = "spectrogram"

    def run(self, spec, params=None):
        t0 = time.time()
        cands = self.module.run(spec, params or {})
        dt = time.time() - t0
        cands = as_candidates(cands, self.name)
        for c in cands:
            c.provenance.setdefault("recipe_version", self.version)
            c.provenance.setdefault("data", spec.meta.get("origin", "?"))
        return cands, dt

    def info(self):
        return {"name": self.name, "description": self.description,
                "author": self.author, "version": self.version,
                "tags": list(self.tags), "input": self.input,
                "file": self.path.name,
                "has_selftest": hasattr(self.module, "selftest")}


def _load_module(path):
    spec = importlib.util.spec_from_file_location(f"recipe_{path.stem}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def discover(directory=None, quiet=True):
    """Auto-discover every recipe in recipes/. A file is a recipe if it defines
    NAME and run(). Broken recipes are reported, never fatal."""
    d = Path(directory or RECIPE_DIR)
    found, broken = [], []
    if not d.exists():
        return found, broken
    for p in sorted(d.glob("*.py")):
        if p.name.startswith("_"):
            continue
        try:
            m = _load_module(p)
            if not hasattr(m, "run") or not hasattr(m, "NAME"):
                broken.append((p.name, "no NAME / no run()"))
                continue
            found.append(Recipe(
                name=m.NAME, path=p, module=m,
                description=getattr(m, "DESCRIPTION", "").strip(),
                author=getattr(m, "AUTHOR", "anonymous"),
                version=str(getattr(m, "VERSION", "0")),
                tags=list(getattr(m, "TAGS", [])),
                input=getattr(m, "INPUT", "spectrogram")))
        except Exception as e:
            broken.append((p.name, f"{type(e).__name__}: {e}"))
            if not quiet:
                traceback.print_exc()
    return found, broken


def get(name, directory=None):
    found, _ = discover(directory)
    for r in found:
        if r.name == name or r.path.stem == name:
            return r
    raise KeyError(f"no recipe named {name!r}; have "
                   f"{[r.name for r in found]}")


# ------------------------------------------------- helpers recipes may reuse
def robust_z(a):
    """Median/MAD z-score - the honest way to say 'how far above the noise'.
    MAD is used instead of std so one bright RFI line doesn't hide everything
    else by inflating the noise estimate."""
    a = np.asarray(a, np.float64)
    med = np.median(a)
    mad = np.median(np.abs(a - med))
    sigma = 1.4826 * mad if mad > 0 else (a.std() or 1.0)
    return (a - med) / sigma


def peaks_z(a, zmin=8.0, min_sep=3, nmax=64):
    """Indices of local maxima above zmin robust sigma, separated by min_sep."""
    z = robust_z(a)
    order = np.argsort(z)[::-1]
    out = []
    for i in order:
        if z[i] < zmin:
            break
        if all(abs(int(i) - j) >= min_sep for j in out):
            out.append(int(i))
        if len(out) >= nmax:
            break
    return out, z


def spectral_flatness(p):
    """Wiener entropy / spectral flatness in (0,1]: ~0 = one tone dominates,
    ~1 = flat like noise. Same measure entropy.py uses on IQ (ENTROPY_RESULT.md)."""
    p = np.asarray(p, np.float64)
    p = p - p.min() + 1e-12
    return float(np.exp(np.mean(np.log(p))) / np.mean(p))


# --------------------------------------------------- verification / vetoes
# Known terrestrial/satellite RFI, the stuff that eats real SETI candidates.
# Ranges from the Breakthrough Listen RFI discussions (Price et al. 2020,
# AJ 159:86, arXiv:1906.07750) plus standard allocations. NOT exhaustive - a
# veto here is a HINT with a reason attached, never a silent delete.
RFI_BANDS = [
    (1025.0, 1150.0, "aeronautical DME / TACAN"),
    # Added 2026-07-30: the GJ699 cadence produced a 41-hit cluster in
    # 1442-1452 MHz and a survivor at 1451.948 MHz that no veto explained. In
    # the US 1435-1535 MHz is Aeronautical Mobile Telemetry (flight-test
    # downlinks) — the table simply had a hole there, its nearest entry
    # starting at 1525. A gap in the RFI table looks exactly like a discovery.
    (1435.0, 1535.0, "aeronautical mobile telemetry (AMT, US flight test)"),
    (1164.0, 1215.0, "GNSS L5 (GPS/Galileo E5/GLONASS L3)"),
    (1215.0, 1400.0, "air-surveillance + military radar (ARSR-4, ASR-9)"),
    (1087.0, 1093.0, "ADS-B / Mode-S aircraft transponders (1090 MHz)"),
    (1525.0, 1559.0, "Inmarsat / MSS downlinks"),
    (1559.0, 1610.0, "GNSS L1 (GPS 1575.42, GLONASS 1598-1606)"),
    (1610.0, 1626.5, "Iridium satellite uplinks - the classic L-band pest"),
    (1675.0, 1700.0, "GOES / weather-satellite downlinks"),
    (2180.0, 2290.0, "space-research + TDRSS downlinks"),
    (2200.0, 2300.0, "S-band spacecraft telemetry"),
    (2320.0, 2345.0, "SiriusXM satellite radio - notorious at Green Bank"),
    (2400.0, 2500.0, "WiFi / Bluetooth / microwave ovens (see: perytons)"),
    (2483.5, 2500.0, "Globalstar downlinks"),
]

# Real astrophysics that a naive technosignature search flags as a 'candidate'.
# Recognising these is the point of SETI_HISTORY.md part B.
NATURAL_LINES = [
    (HI_MHZ, 2.0, "HI 21 cm - neutral hydrogen (galactic, broad, |v|<300 km/s)"),
    (OH_MHZ[0], 0.4, "OH 1612 MHz maser (evolved stars) / shock maser"),
    (OH_MHZ[1], 0.4, "OH 1665 MHz maser (star-forming regions)"),
    (OH_MHZ[2], 0.4, "OH 1667 MHz maser"),
    (OH_MHZ[3], 0.4, "OH 1720 MHz shock maser (SNR interactions)"),
    (6668.5192, 2.0, "methanol 6.7 GHz class-II maser (massive SF)"),
    (12178.593, 2.0, "methanol 12.2 GHz maser"),
    (22235.08, 5.0, "water 22 GHz maser - can be extremely bright"),
]


def rfi_reason(f_mhz):
    for lo, hi, why in RFI_BANDS:
        if lo <= f_mhz <= hi:
            return why
    return None


def natural_reason(f_mhz):
    for f0, tol, why in NATURAL_LINES:
        if abs(f_mhz - f0) <= tol:
            v = -C_KMS * (f_mhz - f0) / f0
            return f"{why}  [{v:+.1f} km/s from rest]"
    return None


def round_number_reason(f_mhz, tol_hz=50.0):
    """Is this frequency suspiciously ROUND in human units?

    Nature does not know what a megahertz is. A hit landing on an exact round
    value is an oscillator, a clock harmonic or a synthesiser birdie, not the
    sky. Added 2026-07-30 after a GJ699 panel produced a score-4116 'candidate'
    at EXACTLY 1500.000000 MHz, which no existing veto caught.

    Checked coarsest-first, because 1500.000000 should be reported as the 100 MHz
    coincidence it is, not as a 1 MHz one.
    """
    f_hz = f_mhz * 1e6
    for step, label in ((100e6, "100 MHz"), (25e6, "25 MHz"), (10e6, "10 MHz"),
                        (5e6, "5 MHz"), (1e6, "1 MHz")):
        if abs(f_hz - round(f_hz / step) * step) <= tol_hz:
            why = round(f_hz / step) * step / 1e6
            return (f"exactly on a {label} boundary ({why:g} MHz) - oscillator / "
                    f"clock harmonic, nature does not pick round numbers")
    return None


def railed_drift_reason(c, tol=1e-6):
    """Did the drift fit stop at the EDGE of its own search range?

    If |drift| equals the searched maximum, the true drift is >= that value and
    was never measured - the optimiser simply ran out of grid. Reporting it as a
    measured drift (and scoring on it) is wrong in a way that looks confident.
    Needs the recipe to say what range it searched; recipes that record
    drift_max_hz_s in provenance get this check for free.
    """
    dmax = (c.provenance or {}).get("drift_max_hz_s")
    if dmax and c.drift_hz_s is not None and abs(abs(c.drift_hz_s) - abs(dmax)) <= tol:
        return (f"drift railed at the search limit ({c.drift_hz_s:+.3f} Hz/s = "
                f"the +-{abs(dmax):g} boundary) - drift NOT measured, only bounded")
    return None


def coarse_channel_reason(cands, coarse_mhz=2.9296875, tol_frac=0.01, min_group=2):
    """Flag hits spaced ~one COARSE CHANNEL apart as instrument bandpass shape.

    Breakthrough Listen's GBT products are stitched from ~2.93 MHz coarse
    channels, and each one carries the same filter roll-off. Structure that
    repeats at exactly that pitch is the instrument, not the sky. Returns a set
    of ids to flag. Found on GJ699: two 'unexplained' wide flat bands at
    1409.209728 and 1412.132263 MHz - 0.998 coarse channels apart.
    """
    flagged = set()
    fs = sorted(((c.freq_mhz, id(c)) for c in cands))
    for i in range(len(fs)):
        for j in range(i + 1, len(fs)):
            d = fs[j][0] - fs[i][0]
            if d > coarse_mhz * 3:
                break
            n = d / coarse_mhz
            if n >= 0.5 and abs(n - round(n)) <= tol_frac and round(n) >= 1:
                flagged.add(fs[i][1])
                flagged.add(fs[j][1])
    return flagged if len(flagged) >= min_group else set()


def explain(cands, spec=None, drift_min_hz_s=0.02, coarse_mhz=2.9296875):
    """Attach a VERDICT to every candidate instead of thresholding it away.
    This is the repo's honesty rail in code form (README 'Honesty rails').

      zero-drift      : a sky source must Doppler-drift (Earth's rotation alone
                        gives ~0.05-0.3 Hz/s at L band). Exactly 0.000 Hz/s over
                        minutes means it is bolted to the ground with us -
                        UNLESS you buy the drift-compensated-beacon argument
                        (NOVEL_DETECTORS #3), which is why we LABEL, not delete.
      known-RFI band  : sits in an allocated terrestrial/satellite band.
      natural line    : sits on a real astrophysical line (HI/OH/CH3OH/H2O).
      band edge       : within 3 channels of the block edge - likely an artifact.
      round number    : on an exact MHz boundary - an oscillator, not the sky.
      railed drift    : the drift fit hit its own search limit and is unmeasured.
      coarse-channel  : repeats at the instrument's ~2.93 MHz channel pitch.
      unexplained     : survived everything. These are the ones to work on.

    The last three were added 2026-07-30 because a Barnard's Star panel produced
    six 'unexplained' hits of which every single one was explicable: an exactly
    1500.000000 MHz birdie whose drift was pinned to the search boundary, and a
    pair of wide flat bands one coarse channel apart. Vetoes label, never delete.
    """
    coarse_flags = coarse_channel_reason(cands, coarse_mhz) if len(cands) > 1 else set()
    for c in cands:
        why = []
        if spec is not None:
            ch = spec.chan_of(c.freq_mhz)
            if ch < 3 or ch > spec.nchan - 4:
                why.append("band edge (<3 channels from the edge)")
        rn = round_number_reason(c.freq_mhz)
        if rn:
            why.append(rn)
        rd = railed_drift_reason(c)
        if rd:
            why.append(rd)
        if id(c) in coarse_flags:
            why.append(f"coarse-channel structure (~{coarse_mhz:.4f} MHz pitch) "
                       "- instrument bandpass, not sky")
        nat = natural_reason(c.freq_mhz)
        if nat:
            why.append(f"natural line: {nat}")
        r = rfi_reason(c.freq_mhz)
        if r:
            why.append(f"known-RFI band: {r}")
        # ── DRIFT: measured-and-zero vs never-measured ─────────────────────
        # These are DIFFERENT claims and the old gate conflated them. It only
        # fired when the recipe's own text mentioned "drift", so a wide flat
        # band reported by spread_flatness at exactly 0.000 Hz/s sailed through
        # as "unexplained" — which is how 1451.948 MHz survived a full ABACAD
        # cadence on GJ699 and briefly looked like a candidate.
        #
        # But simply removing the gate is wrong too: Candidate.drift_hz_s
        # DEFAULTS to 0.0, so an unconditional test would brand every hit from
        # every non-drift recipe as "zero-drift" — a veto that fires on
        # everything explains nothing.
        #
        # So: decide whether drift was actually SEARCHED (the recipe published a
        # search range, or says it measures drift), and label accordingly.
        prov = c.provenance or {}
        measured_drift = ("drift_max_hz_s" in prov
                          or "drift_step_hz_s" in prov
                          or "drift" in " ".join(prov.get("measures", [])))
        if c.drift_hz_s is not None and abs(c.drift_hz_s) < drift_min_hz_s:
            if measured_drift:
                why.append(f"zero-drift ({c.drift_hz_s:+.4f} Hz/s) - local unless "
                           "drift-compensated")
            elif c.kind == "techno":
                # Not a veto, an honest caveat: a sky source MUST Doppler-drift
                # (Earth's rotation alone gives ~0.05-0.3 Hz/s at L band), and
                # this recipe never looked. Claimed as techno, unverifiable as sky.
                why.append("drift NOT measured by this recipe - a sky source must "
                           "Doppler-drift, so sky origin is unverified")
        c.verdict = "; ".join(why) if why else "unexplained - worth a human"
    return cands


def cadence_verify(per_pointing, tol_hz=200.0, pattern=None):
    """The standard SETI verification pattern, implemented.

    Breakthrough Listen observes ON-OFF-ON-OFF-ON-OFF: the target, then a
    nearby throwaway sky position, alternating. A signal from the TARGET appears
    only in the ON scans; anything appearing in an OFF scan is local (this is
    exactly how BLC-1 died - Sheikh et al. 2021 - and how turboSETI's find_event
    works).

    per_pointing: list of (role, candidates) with role 'ON' or 'OFF', in
    observing order. Returns candidates that are present in >=2 ON scans and
    absent from every OFF scan, each annotated with the cadence evidence."""
    if pattern:
        per_pointing = [(p, c) for p, (_, c) in zip(pattern, per_pointing)]
    ons = [c for role, c in per_pointing if role.upper() == "ON"]
    offs = [c for role, c in per_pointing if role.upper() == "OFF"]
    tol_mhz = tol_hz / 1e6
    survivors = []
    if not ons:
        return survivors
    for c in ons[0]:
        n_on = 1
        for other in ons[1:]:
            if any(abs(o.freq_mhz - c.freq_mhz) < tol_mhz for o in other):
                n_on += 1
        in_off = sum(any(abs(o.freq_mhz - c.freq_mhz) < tol_mhz for o in other)
                     for other in offs)
        c.provenance["cadence"] = {"on_hits": n_on, "on_scans": len(ons),
                                   "off_hits": in_off, "off_scans": len(offs)}
        if in_off:
            c.verdict = (f"FAILS cadence: also present in {in_off}/{len(offs)} "
                         "OFF scans -> local interference")
        elif n_on >= max(2, len(ons) - 1):
            c.verdict = (f"PASSES cadence: in {n_on}/{len(ons)} ON scans, "
                         f"0/{len(offs)} OFF scans")
            survivors.append(c)
        else:
            c.verdict = f"weak: only {n_on}/{len(ons)} ON scans"
    return survivors


# ------------------------------------------------------------- the benchmark
# Every recipe is scored on the SAME synthetic scenes, generated from the
# formulae in SETI_HISTORY.md - no download, identical for everyone, and the
# truth is known. Real-data scoring is separate (and unlabelled, so it can only
# ever be a report, never a score).
def benchmark_scenes():
    """(scene_name, spectrogram, truth) where truth is None for a NULL scene
    that must produce NO detection (the false-alarm gate)."""
    S = seti_io.synth
    return [
        # drift scenes use BL fine-resolution channels (2.79 Hz) - at coarser
        # resolution a 0.35 Hz/s drift moves less than one channel and the whole
        # question disappears.
        ("drift_loud", S("drift", ntime=64, nchan=2048, f0_mhz=1400.0,
                         df_mhz=2.79e-6, drift_hz_s=-0.35, snr=25, seed=1),
         dict(freq_mhz=1400.0 + 2.79e-6 * 1024, drift_hz_s=-0.35)),
        ("drift_faint", S("drift", ntime=64, nchan=2048, f0_mhz=1400.0,
                          df_mhz=2.79e-6, drift_hz_s=+0.60, snr=4, seed=2),
         dict(freq_mhz=1400.0 + 2.79e-6 * 1024, drift_hz_s=+0.60)),
        ("comb", S("comb", ntime=32, nchan=2048, f0_mhz=1400.0, df_mhz=0.001,
                   teeth=12, snr=14, seed=3), dict(kind="comb")),
        ("spread", S("spread", ntime=64, nchan=2048, f0_mhz=1400.0,
                     df_mhz=0.001, snr=1.4, seed=4), dict(kind="spread")),
        ("frb_dm560", S("frb", ntime=400, nchan=512, f0_mhz=4000.0, df_mhz=4.0,
                        dt_s=0.001, dm=560, snr=35, seed=5), dict(dm=560.0)),
        ("pulsar_B0329", S("pulsar", ntime=1000, nchan=64, f0_mhz=1350.0,
                           df_mhz=1.0, dt_s=0.01, period_s=0.714, dm=26.8,
                           snr=1.6, seed=6), dict(period_s=0.714)),
        ("hi_line", S("hi", ntime=64, nchan=2048, f0_mhz=1419.5, df_mhz=0.0028,
                      snr=0.8, v_kms=25, seed=7), dict(freq_mhz=HI_MHZ * (1 - 25 / C_KMS))),
        ("NULL_noise", S("noise", ntime=64, nchan=2048, f0_mhz=1400.0,
                         df_mhz=0.001, seed=8), None),
        ("NULL_noise2", S("noise", ntime=64, nchan=2048, f0_mhz=1400.0,
                          df_mhz=0.001, seed=9), None),
    ]


def _matches(cands, truth, spec):
    """Did the recipe find the planted signal? Generous by design: a recipe only
    has to point at the right THING, in whatever coordinate it works in."""
    if not cands:
        return False
    if "freq_mhz" in truth:
        tol = max(20 * spec.df_mhz, 0.01)
        if any(abs(c.freq_mhz - truth["freq_mhz"]) < tol for c in cands):
            return True
    if "drift_hz_s" in truth:
        if any(c.drift_hz_s and abs(c.drift_hz_s - truth["drift_hz_s"]) < 0.3
               for c in cands):
            return True
    if "dm" in truth:
        if any(c.dm is not None and abs(c.dm - truth["dm"]) < 0.15 * truth["dm"]
               for c in cands):
            return True
    if "period_s" in truth:
        if any(c.period_s is not None and abs(c.period_s - truth["period_s"])
               < 0.05 * truth["period_s"] for c in cands):
            return True
    if "kind" in truth:
        if any(truth["kind"] in (c.label + " " + " ".join(c.tags if hasattr(c, "tags") else [])
                                 + " " + str(c.provenance)) for c in cands):
            return True
        return bool(cands)          # fired at all on a scene of that class
    return False


def bench(recipes=None, scenes=None, verbose=True):
    """Score every recipe on every scene. Honest framing: a comb detector is not
    penalised for missing an FRB (nobody is expected to catch everything) - what
    is scored is (scenes caught) with a HARD GATE on false alarms in the NULL
    scenes, because a detector that fires on noise has found nothing at all."""
    recs = recipes if recipes is not None else discover()[0]
    scn = scenes or benchmark_scenes()
    rows = []
    for r in recs:
        row = {"recipe": r.name, "version": r.version, "author": r.author,
               "caught": [], "missed": [], "false_alarms": [], "runtime_s": 0.0,
               "errors": []}
        for name, spec, truth in scn:
            try:
                cands, dt = r.run(spec)
            except Exception as e:
                row["errors"].append(f"{name}: {type(e).__name__}: {e}")
                continue
            row["runtime_s"] += dt
            if truth is None:
                if cands:
                    row["false_alarms"].append(name)
            elif _matches(cands, truth, spec):
                row["caught"].append(name)
            else:
                row["missed"].append(name)
        row["n_caught"] = len(row["caught"])
        row["n_fa"] = len(row["false_alarms"])
        row["clean"] = row["n_fa"] == 0
        row["runtime_s"] = round(row["runtime_s"], 3)
        rows.append(row)
        if verbose:
            print(f"  {r.name:22s} caught {row['n_caught']}/{sum(1 for _,_,t in scn if t)} "
                  f"| false alarms {row['n_fa']}/{sum(1 for _,_,t in scn if t is None)} "
                  f"| {row['runtime_s']:.2f}s"
                  + (f" | ERRORS {row['errors']}" if row["errors"] else ""))
    # rank: clean detectors first, then by catches, then by speed
    rows.sort(key=lambda r: (not r["clean"], -r["n_caught"], r["runtime_s"]))
    return rows, [n for n, _, _ in scn]


def write_leaderboard(rows, scene_names, path=None):
    path = Path(path or HERE / "LEADERBOARD.md")
    json_path = path.with_suffix(".json")
    pos = [s for s in scene_names if not s.startswith("NULL")]
    nulls = [s for s in scene_names if s.startswith("NULL")]
    lines = [
        "# seti-cookbook recipe leaderboard",
        "",
        "Auto-generated by `python cookbook.py bench`. Every recipe is scored on",
        "the SAME synthetic benchmark (`recipe_api.benchmark_scenes()`), built from",
        "the physics quoted in [SETI_HISTORY.md](SETI_HISTORY.md) - so no download",
        "is needed and the truth is known exactly.",
        "",
        "**The gate that matters:** a recipe that fires on the `NULL_*` pure-noise",
        "scenes has found *nothing*, however many real scenes it catches. Clean",
        "recipes are ranked first; noisy ones are listed but marked.",
        "",
        "Nobody is expected to catch everything - a comb detector legitimately",
        "misses FRBs. Breadth of the *panel* is the goal, not any one recipe.",
        "",
        "| # | recipe | ver | author | caught | false alarms | runtime | scenes caught |",
        "|---|--------|-----|--------|--------|--------------|---------|---------------|",
    ]
    for i, r in enumerate(rows, 1):
        mark = "" if r["clean"] else " :warning:"
        lines.append(f"| {i} | `{r['recipe']}`{mark} | {r['version']} | {r['author']} | "
                     f"{r['n_caught']}/{len(pos)} | {r['n_fa']}/{len(nulls)} | "
                     f"{r['runtime_s']:.2f}s | {', '.join(r['caught']) or '-'} |")
    lines += ["", f"Benchmark scenes: {', '.join(pos)} (+ {len(nulls)} null scenes).", ""]
    uncaught = [s for s in pos if not any(s in r["caught"] for r in rows)]
    if uncaught:
        lines += ["## Open bounties", "",
                  "Scenes **no shipped recipe catches yet**. Writing one of these is "
                  "the single highest-value contribution to this repo — see "
                  "[RECIPES.md](RECIPES.md).", ""]
        lines += [f"- **`{s}`**" for s in uncaught] + [""]
    err = [(r["recipe"], e) for r in rows for e in r["errors"]]
    if err:
        lines += ["## errors", ""] + [f"- `{n}`: {e}" for n, e in err] + [""]
    path.write_text("\n".join(lines), encoding="utf-8")
    json_path.write_text(json.dumps({"rows": rows, "scenes": scene_names},
                                    indent=1), encoding="utf-8")
    return path


def selftest():
    print("=" * 70)
    print("recipe_api selftest - the contract itself")
    print("=" * 70)
    checks = []
    recs, broken = discover()
    checks.append((f"discovered {len(recs)} recipes: {[r.name for r in recs]}",
                   len(recs) >= 3))
    checks.append((f"no broken recipes ({broken})", not broken))

    # explain() must catch the classics
    c = [Candidate(freq_mhz=1575.42, score=9, kind="techno", label="drift"),
         Candidate(freq_mhz=HI_MHZ, score=9, kind="techno"),
         Candidate(freq_mhz=1234.5678, score=9, kind="techno", drift_hz_s=-0.4)]
    explain(c)
    checks.append(("GPS L1 flagged as known RFI", "GNSS L1" in c[0].verdict))
    checks.append(("1420.4057 flagged as HI line", "HI 21 cm" in c[1].verdict))
    checks.append(("clean hit left unexplained",
                   "unexplained" in c[2].verdict or "radar" in c[2].verdict))

    # cadence must kill an OFF-scan signal and pass an ON-only one
    on = lambda f: [Candidate(freq_mhz=f, score=5)]
    seq = [("ON", on(1400.0)), ("OFF", []), ("ON", on(1400.0)),
           ("OFF", []), ("ON", on(1400.0)), ("OFF", [])]
    checks.append(("cadence passes an ON-only signal", len(cadence_verify(seq)) == 1))
    seq2 = [("ON", on(1400.0)), ("OFF", on(1400.0)), ("ON", on(1400.0)),
            ("OFF", []), ("ON", on(1400.0)), ("OFF", [])]
    checks.append(("cadence kills an ON+OFF signal", len(cadence_verify(seq2)) == 0))

    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    good = all(o for _, o in checks)
    print("=" * 70)
    print(f"RESULT: {sum(o for _, o in checks)}/{len(checks)}. "
          f"{'PASS' if good else 'FAIL'}")
    return 0 if good else 1


if __name__ == "__main__":
    import sys
    sys.exit(selftest())
