#!/usr/bin/env python3
"""reach.py - how loud would aliens have to be for THIS instrument to hear them?

Turns "we found nothing" into a number: the minimum transmitter power (EIRP,
watts, aimed at us) each instrument could detect at each nearby star. Same
radiometer physics Breakthrough Listen uses (Enriquez et al. 2017, ApJ 849,
104, eqs. 2, 4, 5):

    SEFD   = 2 k Tsys / Aeff                       (Jy)
    S_min  = SNR * SEFD/dnu_t * sqrt(dnu / (npol * t))   (Jy; narrowband,
             signal width dnu_t <= channel width dnu)
    EIRP   = 4 pi d^2 * S_min * dnu_t              (W)

GATE (runs first, refuses to print anything if it fails): with BL's own GBT
numbers (SEFD 10 Jy, SNR 25, 3 Hz channels, 300 s, 2 pols, 1 Hz signal) we
must reproduce Enriquez 2017 Table 3 - S_min 17 Jy and EIRP 5.2e12 W at 50 pc
- within 20 %. (The paper's text quotes 28 GW at 10 ly, which implies a ~1.5 Hz
signal and does not match its own Table 3; we gate on the table.)

    python reach.py                 # the table
    python reach.py --hours 1       # integrate the hobby rows for an hour
    python reach.py --etap ETAP_SUMMARY.json --n-stars 1   # transmitter-fraction bound

TRANSMITTER FRACTION (Margot et al. 2023 / Choza et al. 2024 style): with n
stars searched, a search that recovers a fraction eta_P of plausible signals
(etap.py) over a fraction D of the drift/band space, finding nothing,
bounds the fraction of stars with a detectable transmitter at confidence CL:

    f_TX < [1 - (1 - CL)^(1/n)] / (eta_P * D)

Hobby rows are ESTIMATES with every assumption printed. They become
measurements only when a real capture (hydrogen line, Cas A) confirms the
system temperature - until then treat them as the order of magnitude.
"""
import argparse
import math
import sys

K_B = 1.380649e-23        # J/K
C = 299_792_458.0          # m/s
JY = 1e-26                 # W m^-2 Hz^-1
LY = 9.4607e15             # m
PC = 3.0857e16             # m

F_HZ = 1420.405751768e6    # quote everything at the hydrogen line
LAMBDA = C / F_HZ
L_AO = 1e13                # Arecibo planetary radar EIRP, as Enriquez 2017 uses it


def sefd_jy(tsys_k, aeff_m2):
    return 2 * K_B * tsys_k / aeff_m2 / JY


def s_min_jy(sefd, snr, dnu, t, npol=2, dnu_t=1.0):
    return snr * sefd / dnu_t * math.sqrt(dnu / (npol * t))


def eirp_min_w(s_jy, d_m, dnu_t=1.0):
    return 4 * math.pi * d_m ** 2 * s_jy * JY * dnu_t


def reach_m(s_jy, eirp_w, dnu_t=1.0):
    """Distance at which a transmitter of eirp_w is just detectable."""
    return math.sqrt(eirp_w / (4 * math.pi * s_jy * JY * dnu_t))


def f_tx_upper(n_stars, eta_p, D, cl=0.95):
    """Upper bound on the fraction of stars hosting a transmitter (in the
    searched band/drift space, above the injection strength), from a null
    over n_stars with recovery rate eta_p and space fraction D. A value >= 1
    means no constraint."""
    if n_stars < 1 or eta_p <= 0 or D <= 0:
        return math.inf
    return (1 - (1 - cl) ** (1.0 / n_stars)) / (eta_p * D)


def print_ftx(res, n_stars=1, cl=0.95):
    """res: an etap.py summary (eta_all, D_all, eta_outside_rfi, D_outside_rfi)."""
    print(f"\nTransmitter-fraction bound, CL {cl:.0%}, n = {n_stars} star(s) "
          f"(eta_P = UCLA match, |drift| <= 8.86 Hz/s, one scan, post-channelisation plants)")
    for lab, e, d, ci in (("all band", res["eta_all"], res["D_all"], res.get("eta_all_ci")),
                          ("outside RFI", res["eta_outside_rfi"], res["D_outside_rfi"],
                           res.get("eta_outside_rfi_ci"))):
        f = f_tx_upper(n_stars, e, d, cl)
        rng = ""
        if ci:
            f_hi, f_lo = f_tx_upper(n_stars, ci[0], d, cl), f_tx_upper(n_stars, ci[1], d, cl)
            rng = f"  (eta_P CI -> {f_lo:.3g} .. {f_hi:.3g})"
        verdict = f"f_TX < {f:.3g}" if f < 1 else f"bound {f:.3g} >= 1: NO constraint"
        print(f"  {lab:12s} eta_P {e:.4f}  D {d:.3f}  ->  {verdict}{rng}")
        if lab == "all band" and res.get("eta_all_note"):
            print(f"               note: {res['eta_all_note']}")
    for n in (100, 1000, 10000):
        if n > n_stars:
            print(f"  for scale: the same eta_P, D over {n:>5d} stars -> f_TX < "
                  f"{f_tx_upper(n, res['eta_all'], res['D_all'], cl):.3g} (all band), "
                  f"{f_tx_upper(n, res['eta_outside_rfi'], res['D_outside_rfi'], cl):.3g} (outside RFI)")


def gate():
    s = s_min_jy(10, 25, 3, 300, 2, 1.0)
    e = eirp_min_w(s, 50 * PC)
    ok = abs(s / 17 - 1) < 0.2 and abs(e / 5.2e12 - 1) < 0.2
    print(f"GATE vs Enriquez 2017 Table 3: S_min {s:.1f} Jy (paper 17), "
          f"EIRP@50pc {e:.2e} W (paper 5.2e12)  ->  {'PASS' if ok else 'FAIL'}")
    return ok


def aeff_from_gain_dbi(g_dbi):
    return 10 ** (g_dbi / 10) * LAMBDA ** 2 / (4 * math.pi)


# name, Aeff m^2, Tsys K, channel Hz, how the numbers were chosen
INSTRUMENTS = [
    ("GBT 100 m (BL survey)", None, None, 3.0,
     "SEFD 10 Jy as published (Tsys ~20 K, eta 0.72)"),
    ("RSPdx + roof discone", aeff_from_gain_dbi(0.0), 500.0, 1.0,
     "0 dBi; Tsys ~500 K = RSPdx NF ~4.5 dB at 1.4 GHz + warm ground, no LNA"),
    ("21 cm horn + LNA", aeff_from_gain_dbi(17.0), 60.0, 1.0,
     "17 dBi horn; Tsys ~60 K = ~40 K LNA (NF ~0.5 dB) + sky + spillover"),
    ("used 3 m TVRO dish + LNA", 0.5 * math.pi * 1.5 ** 2, 50.0, 1.0,
     "eta 0.5; Tsys ~50 K"),
]

STARS = [("Barnard's Star", 5.96), ("Wolf 359", 7.86), ("Lalande 21185", 8.31),
         ("Sirius", 8.6), ("Eps Eridani", 10.5), ("Tau Ceti", 11.9)]


def fmt_w(w):
    for unit, k in (("EW", 1e18), ("PW", 1e15), ("TW", 1e12), ("GW", 1e9), ("MW", 1e6)):
        if w >= k:
            return f"{w / k:6.1f} {unit}"
    return f"{w:8.0f} W"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--hours", type=float, default=None,
                    help="integration for the hobby rows (default: 5 min, like BL)")
    ap.add_argument("--snr", type=float, default=25.0)
    ap.add_argument("--etap", metavar="JSON",
                    help="etap.py summary JSON: print the f_TX upper bound (all band and outside RFI)")
    ap.add_argument("--n-stars", type=int, default=1, help="stars searched, for --etap")
    ap.add_argument("--cl", type=float, default=0.95, help="confidence level, for --etap")
    args = ap.parse_args()

    if not gate():
        sys.exit("refusing to print a reach table the gate does not trust")
    if args.etap:
        import json
        with open(args.etap) as fh:
            print_ftx(json.load(fh), n_stars=args.n_stars, cl=args.cl)
        return

    print(f"\nMinimum detectable EIRP at {F_HZ/1e6:.3f} MHz, SNR {args.snr:g}, 2 pols, 1 Hz-wide signal.")
    print(f"Reference: Arecibo planetary radar ~{fmt_w(L_AO).strip()} (the loudest beam Earth has sent).\n")
    head = f"{'instrument':28s} {'SEFD':>10s} {'t':>6s} " + " ".join(f"{n[:13]:>13s}" for n, _ in STARS) + "   hears Arecibo to"
    print(head)
    print("-" * len(head))
    for name, aeff, tsys, dnu, why in INSTRUMENTS:
        sefd = 10.0 if aeff is None else sefd_jy(tsys, aeff)
        t = 300.0 if (aeff is None or args.hours is None) else args.hours * 3600
        s = s_min_jy(sefd, args.snr, dnu, t)
        cells = " ".join(f"{fmt_w(eirp_min_w(s, d * LY)):>13s}" for _, d in STARS)
        r_ly = reach_m(s, L_AO) / LY
        rtxt = f"{r_ly:,.0f} ly" if r_ly >= 1 else f"{r_ly * 63241:,.0f} AU"
        tl = f"{t/60:.0f}m" if t < 3600 else f"{t/3600:.0f}h"
        print(f"{name:28s} {sefd:9.3g}J {tl:>6s} {cells}   {rtxt}")
    print("\nAssumptions:")
    for name, aeff, tsys, dnu, why in INSTRUMENTS:
        a = "" if aeff is None else f"Aeff {aeff:.3g} m^2, "
        print(f"  {name}: {a}{why}")
    print("\nRead it as: to be heard by that instrument from that star, a transmitter")
    print("aimed at Earth must be at least this powerful. Hobby rows are estimates until")
    print("a hydrogen-line or Cas A capture measures the real system temperature.")


if __name__ == "__main__":
    main()
