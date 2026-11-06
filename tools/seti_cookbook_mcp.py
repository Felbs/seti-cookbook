"""seti_cookbook_mcp.py - seti-cookbook as an MCP server. OPTIONAL COMPANION, OFF BY DEFAULT.

Exposes seti-cookbook's real verbs to any MCP client (Claude Code, Claude Desktop, a
local LLM) as typed tools, so a language model can drive a SETI search: list the
data you have, run detectors and recipes on it, read the candidate list with its
verdicts, render a waterfall, and score recipes on the shared benchmark.

WHY: the bottleneck in this project is IMAGINATION, not telescope time
(NOVEL_DETECTORS.md). An LLM that can (a) read the recipe contract, (b) write a
new recipe file, and (c) immediately run and score it against the same public
data as everyone else, closes the hypothesis->experiment->conclusion loop that the
rest of the repo does by hand. That is what this is for.

PRIVACY / SCOPE - read this before you enable it:
  * Nothing here is required. seti-cookbook is a complete command-line tool without
    it; this file is a separate optional companion and no other file imports it.
  * It is a LOCAL process speaking MCP over stdio to a client on the same
    machine. It opens no network port and phones nothing home.
  * It sends nothing anywhere ITSELF. The only outbound traffic any tool here can
    cause is the same public Breakthrough Listen download that
    `python star_sweep.py` already does, and only when you call the fetch tool.
  * Your MCP CLIENT is a different question and not ours: if the client is a
    cloud LLM, then whatever these tools RETURN (file names, candidate lists,
    numbers) goes to that provider, exactly as if you had pasted it into a chat.
    Point it at a local model if that matters to you.
  * Filesystem access is confined to the repo by default: paths are resolved
    under the seti-cookbook directory and anything outside is refused unless you set
    SETI_COOKBOOK_MCP_ALLOW_ANY_PATH=1.
  * No telescope, no SDR, no radio hardware is touched by any tool here. seti-cookbook
    reads archive files.

Run standalone:      python tools/seti_cookbook_mcp.py
Claude Code config:  see .mcp.json.example at the repo root
Requires:            pip install fastmcp   (plus the repo's own requirements)
"""
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO))

from fastmcp import FastMCP                                       # noqa: E402

import recipe_api as R                                            # noqa: E402
import seti_io                                                    # noqa: E402

mcp = FastMCP("seti-cookbook")

ALLOW_ANY_PATH = os.environ.get("SETI_COOKBOOK_MCP_ALLOW_ANY_PATH", os.environ.get("SETITUNA_MCP_ALLOW_ANY_PATH", "0")) in ("1", "true")


def _safe(path):
    """Resolve a path inside the repo. Refuse anything outside unless the user
    explicitly opted out - an LLM should not be able to read /etc/passwd through
    a SETI tool."""
    s = str(path)
    if s.startswith("synth:"):
        return s
    p = Path(s)
    p = p if p.is_absolute() else (REPO / p)
    p = p.resolve()
    if not ALLOW_ANY_PATH:
        try:
            p.relative_to(REPO)
        except ValueError:
            raise ValueError(
                f"refusing to read outside the seti-cookbook repo ({REPO}). "
                "Set SETI_COOKBOOK_MCP_ALLOW_ANY_PATH=1 if you really mean to.")
    return str(p)


def _open(path, f_start=None, f_stop=None, t_start=None, t_stop=None, fs=None):
    kw = {}
    for k, v in (("f_start", f_start), ("f_stop", f_stop),
                 ("t_start", t_start), ("t_stop", t_stop), ("fs", fs)):
        if v is not None:
            kw[k] = v
    p = _safe(path)
    if str(p).startswith("synth:") or Path(str(p)).suffix.lower() in (".cs16", ".iq"):
        kw.pop("f_start", None)
        kw.pop("f_stop", None)
        kw.pop("t_start", None)
        kw.pop("t_stop", None)
    return seti_io.open_any(p, **kw)


# --------------------------------------------------------------- what's here
@mcp.tool
def list_data() -> dict:
    """Every data file seti-cookbook can open in this checkout (data/ is gitignored -
    these are files YOU downloaded from the public Breakthrough Listen archive).
    Also lists the synthetic scenes available without any download."""
    files = seti_io.list_data()
    return {
        "files": files,
        "note": "data/ is gitignored; no third-party data ships in this repo",
        "synthetic": [f"synth:{k}" for k in
                      ("noise", "drift", "zerodrift", "frb", "pulsar", "comb",
                       "spread", "maser", "hi")],
        "synthetic_usage": "synth:drift,drift_hz_s=-0.35,snr=25 - parameters are "
                           "comma-separated key=value pairs",
        "archive": "http://seti.berkeley.edu/opendata",
    }


@mcp.tool
def data_info(path: str, f_start: float = None, f_stop: float = None) -> dict:
    """Header and summary of one data product: source name, telescope, band,
    channel resolution, duration, MJD, and the pointing (RA/Dec) if present."""
    try:
        s = _open(path, f_start, f_stop)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    out = s.summary()
    out["repr"] = repr(s)
    out["ra_hr"] = s.meta.get("ra_hr")
    out["dec_deg"] = s.meta.get("dec_deg")
    return out


@mcp.tool
def list_recipes() -> dict:
    """Every detector RECIPE discovered in recipes/ - the pluggable
    'invent a new way to find aliens' format (see RECIPES.md). Includes broken
    files, with the reason, so nothing fails silently."""
    recs, broken = R.discover()
    return {"recipes": [r.info() for r in recs],
            "broken": [{"file": f, "why": w} for f, w in broken],
            "contract": "a recipe declares NAME/DESCRIPTION/AUTHOR/VERSION/TAGS "
                        "and defines run(spec, params) -> [Candidate]; see "
                        "RECIPES.md and recipe_api.py"}


@mcp.tool
def targets_available(target: str) -> dict:
    """Ask the public Breakthrough Listen open-data archive what files exist for
    a target name (e.g. GJ699, HIP54035, Voyager1). Read-only query; downloads
    nothing. This is the one tool that touches the network."""
    import urllib.request
    import urllib.parse
    url = ("http://seti.berkeley.edu/opendata/api/query-files?"
           + urllib.parse.urlencode({"target": target, "telescopes": "GBT",
                                     "file-types": "HDF5", "limit": 12}))
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            rows = json.loads(r.read()).get("data", [])
    except Exception as e:
        return {"error": f"BL archive not reachable ({e})", "url": url}
    return {"target": target, "n": len(rows),
            "files": [{"url": x.get("url"), "size_mb": round(x.get("size", 0) / 1e6, 1)}
                      for x in rows],
            "note": "fetch with star_sweep.py; files land in data/ (gitignored)"}


# ------------------------------------------------------------------ detecting
@mcp.tool
def run_recipe(recipe: str, path: str, f_start: float = None, f_stop: float = None,
               t_start: float = None, t_stop: float = None, fs: float = None,
               params: dict = None, top: int = 12) -> dict:
    """Run ONE recipe on one data product and return its candidates, each with a
    verdict from the verification layer (known-RFI band / natural spectral line /
    zero-drift / band edge / unexplained).

    recipe: a name from list_recipes(), e.g. narrowband_drift, comb_uniformity,
            spread_flatness, dispersion_sweep, hi_line_natural, pi_ratio
    path:   a file from list_data(), or a synth: scene
    f_start/f_stop in MHz narrow the search (needed for big BL files);
    fs is required for raw .cs16 IQ."""
    try:
        r = R.get(recipe)
    except KeyError as e:
        return {"error": str(e)}
    try:
        s = _open(path, f_start, f_stop, t_start, t_stop, fs)
        cands, dt = r.run(s, params or {})
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    R.explain(cands, s)
    cands.sort(key=lambda c: -c.score)
    return {"recipe": r.name, "version": r.version, "data": s.summary(),
            "runtime_s": round(dt, 3), "n_candidates": len(cands),
            "n_unexplained": sum(c.verdict.startswith("unexplained") for c in cands),
            "candidates": [c.to_dict() for c in cands[:top]]}


@mcp.tool
def run_panel(path: str, f_start: float = None, f_stop: float = None,
              t_start: float = None, t_stop: float = None, fs: float = None,
              top_per_recipe: int = 4) -> dict:
    """Run EVERY recipe on one data product - 'try every way of being
    artificial' (agent.py's ensemble idea, applied to the recipe cookbook).
    Returns per-recipe candidate lists plus a combined verdict count."""
    try:
        s = _open(path, f_start, f_stop, t_start, t_stop, fs)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    recs, broken = R.discover()
    out = {}
    total = unexplained = 0
    for r in recs:
        try:
            cands, dt = r.run(s)
        except Exception as e:
            out[r.name] = {"error": f"{type(e).__name__}: {e}"}
            continue
        R.explain(cands, s)
        cands.sort(key=lambda c: -c.score)
        total += len(cands)
        unexplained += sum(c.verdict.startswith("unexplained") for c in cands)
        out[r.name] = {"runtime_s": round(dt, 3), "n": len(cands),
                       "candidates": [c.to_dict() for c in cands[:top_per_recipe]]}
    return {"data": s.summary(), "per_recipe": out,
            "n_candidates": total, "n_unexplained": unexplained,
            "broken_recipes": [f for f, _ in broken]}


@mcp.tool
def novel_detector_panel(iq_path: str, fs: float) -> dict:
    """Run the ORIGINAL novel-detector suite (cyclo + comb + entropy) on a raw
    interleaved-int16 IQ capture. These work on voltages, not spectrograms, so
    they can see cyclostationary structure a waterfall cannot
    (CYCLO_RESULT.md / COMB_RESULT.md / ENTROPY_RESULT.md)."""
    import numpy as np
    import cyclo
    import comb
    import entropy
    try:
        p = _safe(iq_path)
        raw = np.fromfile(p, np.int16, count=2 * 4_000_000).astype(np.float32) / 32768.0
        x = (raw[0::2] + 1j * raw[1::2]).astype(np.complex64)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    if not len(x):
        return {"error": "no samples read"}
    a, sig, tau = cyclo.detect(x, fs)
    sp, csig, depth = comb.detect(x, fs)
    sfm, label = entropy.detect(x, fs)
    return {
        "file": Path(iq_path).name, "fs_hz": fs, "seconds": len(x) / fs,
        "cyclo": {"alpha_hz": round(float(a), 2), "significance": round(float(sig), 2),
                  "tau": int(tau), "fired": bool(sig >= 8.0),
                  "means": "hidden symbol/chip-rate periodicity = a digital signal, "
                           "even one spread below the noise floor"},
        "comb": {"spacing_hz": round(float(sp), 2), "significance": round(float(csig), 2),
                 "fired": bool(csig >= 10.0),
                 "means": "uniformly Hz-spaced tones = an engineered frequency reference"},
        "entropy": {"spectral_flatness": round(float(sfm), 4), "class": label,
                    "fired": bool(entropy.LO <= sfm <= entropy.HI),
                    "means": "structure between a dead carrier and pure noise"},
    }


@mcp.tool
def cadence_check(recipe: str, pointings: list, tol_hz: float = 200.0,
                  f_start: float = None, f_stop: float = None) -> dict:
    """The standard SETI verification test, as a tool. Give an ordered list like
    [["ON","a.h5"],["OFF","b.h5"],["ON","c.h5"],["OFF","d.h5"]] and it keeps only
    candidates present in the ON scans and absent from every OFF scan. A signal
    seen in an OFF scan is local interference - this is exactly how BLC-1 was
    resolved (Sheikh et al. 2021, Nature Astronomy)."""
    try:
        r = R.get(recipe)
    except KeyError as e:
        return {"error": str(e)}
    seq = []
    per = []
    for item in pointings:
        role, path = (item[0], item[1]) if isinstance(item, (list, tuple)) \
            else str(item).split("=", 1)
        try:
            s = _open(path, f_start, f_stop)
            cands, _ = r.run(s)
        except Exception as e:
            return {"error": f"{path}: {type(e).__name__}: {e}"}
        R.explain(cands, s)
        seq.append((role, cands))
        per.append({"role": role.upper(), "file": Path(path).name,
                    "n_candidates": len(cands)})
    survivors = R.cadence_verify(seq, tol_hz=tol_hz)
    return {"recipe": r.name, "pointings": per,
            "n_survivors": len(survivors),
            "survivors": [c.to_dict() for c in survivors],
            "note": "0 survivors is the normal, honest result"}


# ------------------------------------------------------------------ scoring
@mcp.tool
def benchmark(write_leaderboard: bool = True) -> dict:
    """Score every recipe on seti-cookbook's shared synthetic benchmark and return the
    ranking. Clean recipes (zero false alarms on the pure-noise NULL scenes) rank
    first; a detector that fires on noise has found nothing. Optionally rewrites
    LEADERBOARD.md."""
    rows, scenes = R.bench(verbose=False)
    path = str(R.write_leaderboard(rows, scenes)) if write_leaderboard else None
    uncaught = sorted(set(s for s in scenes if not s.startswith("NULL"))
                      - set(s for r in rows for s in r["caught"]))
    return {"scenes": scenes, "ranking": rows, "leaderboard": path,
            "uncaught_scenes": uncaught,
            "open_bounty": ("no recipe catches these yet - writing one is the "
                            "highest-value contribution" if uncaught else None)}


@mcp.tool
def leaderboard() -> dict:
    """The current recipe leaderboard as last written to LEADERBOARD.json."""
    p = REPO / "LEADERBOARD.json"
    if not p.exists():
        return {"error": "no leaderboard yet - call benchmark() first"}
    return json.loads(p.read_text())


@mcp.tool
def recipe_source(recipe: str) -> dict:
    """The full source of one recipe. Useful for an LLM that is about to write a
    new one: read a working example, then copy the shape."""
    try:
        r = R.get(recipe)
    except KeyError as e:
        return {"error": str(e)}
    return {"name": r.name, "file": str(r.path.relative_to(REPO)),
            "source": r.path.read_text(encoding="utf-8")}


@mcp.tool
def recipe_contract() -> dict:
    """The recipe API contract and the tutorial, verbatim - everything needed to
    write a new detector that this server can immediately run and score."""
    out = {"api_module": "recipe_api.py",
           "api_docstring": R.__doc__,
           "candidate_fields": list(R.Candidate.__dataclass_fields__.keys()),
           "recipe_dir": str(R.RECIPE_DIR.relative_to(REPO)),
           "how_to_add": "drop a .py file in recipes/ defining NAME and "
                         "run(spec, params); it is auto-discovered"}
    doc = REPO / "RECIPES.md"
    if doc.exists():
        out["RECIPES_md"] = doc.read_text(encoding="utf-8")
    return out


# --------------------------------------------------------------- looking at it
@mcp.tool
def render_waterfall(path: str, out: str = None, f_start: float = None,
                     f_stop: float = None, t_start: float = None,
                     t_stop: float = None, fs: float = None,
                     drift_hz_s: float = None, dm: float = None,
                     mark_mhz: list = None, log_freq: bool = False,
                     normalize: bool = False, view: str = "plot",
                     period_s: float = None) -> dict:
    """Render a spectrogram PNG and return its path.

    view='plot'   waterfall + integrated spectrum, with optional overlays:
                  drift_hz_s draws the technosignature drift line, dm draws the
                  FRB/pulsar dispersion sweep, mark_mhz drops candidate markers,
                  normalize divides out the instrument bandpass.
    view='hough'  drift-rate vs frequency plane: a sky tone sits off the
                  zero-drift row, local interference sits on it.
    view='fold'   phase-fold at period_s - how a pulsar becomes visible.

    Writes into figures/ by default. The file is local; nothing is uploaded."""
    import waterfall as W
    try:
        s = _open(path, f_start, f_stop, t_start, t_stop, fs)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    stem = (Path(str(path).replace(":", "_")).stem or "scene")
    target = Path(_safe(out)) if out else (REPO / "figures" / f"{stem}_{view}.png")
    try:
        if view == "hough":
            W.hough(s, target)
        elif view == "fold":
            if not period_s:
                return {"error": "view='fold' needs period_s"}
            W.fold(s, target, period_s, dm=dm or 0.0)
        else:
            W.plot(s, target, drift=[drift_hz_s] if drift_hz_s else None,
                   dm=[dm] if dm else None, marks=mark_mhz,
                   log_freq=log_freq, norm=normalize)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    return {"figure": str(target), "size_kb": round(target.stat().st_size / 1e3, 1),
            "view": view, "data": s.summary()}


@mcp.tool
def natural_signals(path: str) -> dict:
    """What ASTROPHYSICS is in this file? Measures the galactic 21 cm hydrogen
    line (with barycentric/LSR-corrected velocity when astropy is installed),
    identifies the spectrometer's own coarse-channel artifacts, and reports the
    maser lines the band covers. This is the 'what SETI finds instead of aliens'
    tool - see SETI_HISTORY.md."""
    import natural_signals as N
    try:
        return N.report(_safe(path), verbose=False)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool
def sky_or_instrument(paths: list, tol_khz: float = 20.0) -> dict:
    """Given several pointings, separate SKY features from INSTRUMENT features:
    a line at the same topocentric frequency in independent pointings on
    different days cannot be the sky. The oldest test in radio astronomy, and the
    same logic as the ON/OFF cadence."""
    import natural_signals as N
    try:
        verdicts, per = N.discriminate([_safe(p) for p in paths],
                                       tol_khz=tol_khz, verbose=False)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    return {"verdicts": verdicts,
            "pointings": [{"file": n, "source": m.get("source_name"),
                           "l_deg": f.get("l_deg"), "b_deg": f.get("b_deg"),
                           "n_lines": len(c)} for n, m, f, c in per]}


@mcp.tool
def selftest() -> dict:
    """Run the whole repo's regression gates (seti_io, recipe_api, every recipe's
    own selftest, and the original novel-detector suite). Proof the tools work
    before you trust an answer."""
    res = {}
    for cmd, label in (([sys.executable, "seti_io.py", "selftest"], "seti_io"),
                       ([sys.executable, "recipe_api.py"], "recipe_api"),
                       ([sys.executable, "cookbook.py", "selftest"], "recipes"),
                       ([sys.executable, "agent.py", "selftest"], "novel_detectors")):
        try:
            p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                               timeout=600)
            res[label] = {"pass": p.returncode == 0,
                          "tail": (p.stdout or p.stderr)[-700:]}
        except Exception as e:
            res[label] = {"pass": False, "tail": f"{type(e).__name__}: {e}"}
    res["all_pass"] = all(v["pass"] for v in res.values() if isinstance(v, dict))
    return res


# ------------------------------------------------------- getting more haystack
@mcp.tool
def fetch_data(what: str = "list", target: str = None, timeout_s: int = 1800) -> dict:
    """DOWNLOAD public telescope data into data/ so there is something to search.

    This is the tool that actually moves the needle: the repo ships no data, and a
    detector with one file is a demo, not a search. Wraps fetch_public_data.py,
    which only ever pulls from the public archives it credits.

    what: 'list' (default, downloads nothing and shows the catalogue) |
          'voyager' (the standing calibration signal) |
          'bl' (a Breakthrough Listen nearby-star pointing; pass target=) |
          'callisto' (solar radio bursts, small) | 'frb' | 'pulsar'

    NETWORK: yes, this one downloads. Everything else in this server is local.
    Files land in data/, which is gitignored. Sizes are in the catalogue — check
    with what='list' before pulling a BL pointing.
    """
    cmd = [sys.executable, "fetch_public_data.py", what]
    if target:
        cmd += ["--target", target]
    try:
        pr = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                            timeout=timeout_s)
    except subprocess.TimeoutExpired:
        return {"error": f"timed out after {timeout_s}s — a BL pointing can be "
                         f"several GB; raise timeout_s or fetch from a shell"}
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    have = sorted(x.name for x in (REPO / "data").glob("*")
                  if x.suffix in (".h5", ".fil", ".dat", ".cs16"))
    return {"what": what, "target": target, "ok": pr.returncode == 0,
            "stdout_tail": (pr.stdout or "")[-1500:],
            "stderr_tail": (pr.stderr or "")[-600:],
            "data_dir_now": have,
            "next": "run_recipe / run_panel on the new file, then ledger_record "
                    "so it counts toward coverage"}


# ------------------------------------------------------- authoring a detector
@mcp.tool
def propose_recipe(name: str, source: str, overwrite: bool = False) -> dict:
    """WRITE a new detector recipe and immediately prove it on the controls.

    The gap this closes: an LLM could read the contract and run recipes, but not
    AUTHOR one, so every new detector had to be written outside the loop. Now the
    hypothesis -> experiment -> conclusion cycle closes inside MCP.

    THE HOUSE RULE IS ENFORCED MECHANICALLY, NOT SUGGESTED: the file is staged
    under a temporary name, its selftest is run, and it is installed into
    recipes/ ONLY if that selftest passes. A detector that fires on pure noise
    has found nothing (RECIPES.md), so this refuses to install it and says why.
    Nothing is ever silently kept.

    name:   recipe name, becomes recipes/<name>.py
    source: the complete Python file (see recipe_contract() for the format)
    """
    import re as _re
    if not _re.fullmatch(r"[a-z][a-z0-9_]{2,40}", name or ""):
        return {"installed": False,
                "why": "name must be lower_snake_case, 3-41 chars"}
    final = REPO / "recipes" / f"{name}.py"
    if final.exists() and not overwrite:
        return {"installed": False,
                "why": f"recipes/{name}.py exists; pass overwrite=True to replace"}
    staged = REPO / "recipes" / f"_staged_{name}.py"
    try:
        staged.write_text(source, encoding="utf-8")
    except Exception as e:
        return {"installed": False, "why": f"could not write: {e}"}
    try:
        pr = subprocess.run([sys.executable, staged.name],
                            cwd=str(REPO / "recipes"),
                            capture_output=True, text=True, timeout=900)
        out = ((pr.stdout or "") + (pr.stderr or ""))[-1800:]
        if pr.returncode != 0:
            return {"installed": False,
                    "why": "selftest FAILED — not installed. A recipe must find "
                           "its own signal AND stay silent on noise and on the "
                           "WRONG signal (RECIPES.md, 'the one rule').",
                    "selftest_output": out}
        staged.replace(final)
        recs, _broken = R.discover()
        return {"installed": True, "path": f"recipes/{name}.py",
                "selftest_output": out,
                "discovered": any(r.info().get("name") == name for r in recs),
                "next": "benchmark() to score it on the shared scenes, then "
                        "run_recipe on real data"}
    except subprocess.TimeoutExpired:
        return {"installed": False, "why": "selftest timed out (>900 s)"}
    finally:
        if staged.exists():
            try:
                staged.unlink()
            except OSError:
                pass


# ------------------------------------------------------------ what is covered
@mcp.tool
def ledger_record(path: str, recipe: str, target: str = None,
                  n_candidates: int = 0, n_unexplained: int = 0,
                  note: str = "") -> dict:
    """Record that a file was searched with a recipe, so the search has a
    denominator. A null result is only worth something if you can say what you
    covered — this is what turns 'we found nothing' into a limit."""
    sys.path.insert(0, str(REPO))
    import search_ledger
    f = _safe(path)
    meta = {}
    try:
        sp = seti_io.open_any(str(f))
        meta = dict(f_lo_mhz=float(sp.f_lo), f_hi_mhz=float(sp.f_hi),
                    duration_s=float(sp.duration_s),
                    res_hz=float(abs(sp.res_hz)),
                    target=target or (sp.meta or {}).get("source_name"))
    except Exception:
        meta = dict(target=target)
    tgt = meta.pop("target", target)
    row = search_ledger.record(str(f), recipe, target=tgt,
                               n_candidates=n_candidates,
                               n_unexplained=n_unexplained, note=note, **meta)
    return {"recorded": row, "coverage_now": search_ledger.coverage()}


@mcp.tool
def coverage() -> dict:
    """HOW MUCH OF THE HAYSTACK HAVE WE LIFTED? Percent of the public archive
    searched, targets covered, integrated bandwidth and time.

    Reports several denominators and deliberately refuses to collapse them into
    one number: bytes-of-archive is honest about scale but measures EFFORT, not
    insight — two recipes on the same file cover different signal space. The
    archive size is an order-of-magnitude estimate (~1 PB; Lebofsky et al. 2019,
    Price et al. 2020) and is overridable via SETI_ARCHIVE_BYTES."""
    sys.path.insert(0, str(REPO))
    import search_ledger
    return search_ledger.coverage()


@mcp.tool
def gpu_status() -> dict:
    """Is a GPU actually usable, and is it worth it? seti-cookbook's array module
    switches to cupy only when SETI_COOKBOOK_GPU=1 AND a device really exists."""
    import numpy as _np
    info = {"env_SETI_COOKBOOK_GPU": os.environ.get("SETI_COOKBOOK_GPU", os.environ.get("SETITUNA_GPU", "0")),
            "measured_speedup_note":
                "doppler_coherence on the 16x1048576 Voyager file: 92.1 s CPU "
                "vs 4.4 s GPU = 21x, identical result. NOTE the CPU and GPU "
                "paths use DIFFERENT algorithms on purpose — an index-gather "
                "beats np.roll on a GPU and loses to it badly on a CPU (204 s)."}
    try:
        import cupy
        cupy.zeros(1)
        prop = cupy.cuda.runtime.getDeviceProperties(0)
        info.update(available=True, cupy=cupy.__version__,
                    device=prop["name"].decode(),
                    active=seti_io.xp() is not _np)
    except Exception as e:
        info.update(available=False, why=f"{type(e).__name__}: {e}", active=False)
    return info


if __name__ == "__main__":
    mcp.run()
