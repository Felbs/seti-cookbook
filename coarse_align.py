#!/usr/bin/env python3
"""coarse_align.py - cut data on TRUE coarse-channel boundaries.

Green Bank (BL) files are channelised in two stages: a polyphase filterbank
makes coarse channels of 2.9296875 MHz, then each is split into fine channels.
Every coarse channel has its own bandpass roll-off at both edges and a DC spike
at its centre. A search slice that crosses a coarse boundary therefore holds a
roll-off notch in the middle and a DC spike where the search does not expect
it (KNOWN_ISSUES #8).

How many fine channels make one coarse channel (the rule blimpy uses):
  * nchans a multiple of 2^20 and >= 2^20     -> 2^20 per coarse channel;
  * otherwise, Green Bank (telescope_id 6)    -> 2.9296875 MHz / |foff|
    (1,033,216 for the 2016 GJ699 fine-res files, 1,024 for their mid-res
    .0002 files);
  * from foff alone (an npz slice without nchans) the same 2.9296875 MHz / |foff|
    when it is an integer, else 2^20.

    python coarse_align.py selftest
"""
import sys

import numpy as np

GBT_COARSE_MHZ = 2.9296875


def nfine_per_coarse(foff_mhz, nchans=None, telescope_id=6):
    """Fine channels per coarse channel (blimpy's rule, see module doc)."""
    if nchans is not None and nchans >= 2 ** 20 and nchans % 2 ** 20 == 0:
        return 2 ** 20
    r = GBT_COARSE_MHZ / abs(float(foff_mhz))
    if (telescope_id in (6, None)) and abs(r - round(r)) < 1e-3:
        n = int(round(r))
        if nchans is None or nchans % n == 0:
            return n
    return 2 ** 20


def segment_bounds(ch0, nf, nfpc):
    """Largest piece of a slice [ch0, ch0+nf) (file channel numbers) that lies inside
    ONE coarse channel. Returns (lo, hi, coarse_index, dc_in_segment or None),
    lo/hi relative to the slice."""
    first_cc = ch0 // nfpc
    pieces = []
    c = first_cc
    while c * nfpc < ch0 + nf:
        a = max(ch0, c * nfpc)
        b = min(ch0 + nf, (c + 1) * nfpc)
        if b > a:
            pieces.append((b - a, a - ch0, b - ch0, c))
        c += 1
    _, lo, hi, cc = max(pieces)
    dc = cc * nfpc + nfpc // 2 - ch0 - lo
    return lo, hi, cc, (dc if 0 <= dc < hi - lo else None)


def coarse_segment(data, ch0, foff_mhz, nchans=None, repair_dc=True):
    """data (nt, nf) slice starting at file channel ch0 -> (seg, info).
    seg holds only the channels inside one coarse channel; its DC channel (if
    inside) is replaced by the mean of its neighbours, as turboSETI's blank_dc does.
    info: lo, hi (in the slice), coarse, dc (in seg), nfpc, dropped (channels cut),
    straddled (bool)."""
    nf = data.shape[1]
    nfpc = nfine_per_coarse(foff_mhz, nchans)
    lo, hi, cc, dc = segment_bounds(int(ch0), nf, nfpc)
    seg = np.ascontiguousarray(data[:, lo:hi])
    if repair_dc and dc is not None and 0 < dc < seg.shape[1] - 1:
        seg[:, dc] = 0.5 * (seg[:, dc - 1] + seg[:, dc + 1])
    info = dict(lo=lo, hi=hi, coarse=cc, dc=dc, nfpc=nfpc, dropped=nf - (hi - lo),
                straddled=(hi - lo) < nf)
    return seg, info


def describe(info, nf):
    if not info["straddled"]:
        return f"slice lies inside coarse channel {info['coarse']} ({nf} ch)"
    return (f"slice straddles coarse channels: using {info['hi'] - info['lo']} of {nf} channels "
            f"(slice ch {info['lo']}-{info['hi']}) inside coarse channel {info['coarse']} of "
            f"{info['nfpc']} ch; DC spike at segment ch {info['dc']}")


def selftest():
    checks = []
    foff = -2.835503418452676e-06
    checks.append(("GJ699 fine-res: 1,033,216 per coarse from foff",
                   nfine_per_coarse(foff) == 1033216))
    checks.append(("GJ699 fine-res with nchans 264,503,296", nfine_per_coarse(foff, 264503296) == 1033216))
    checks.append(("GJ699 mid-res: 1,024 per coarse", nfine_per_coarse(-0.00286102294921875, 262144) == 1024))
    checks.append(("2^20 multiple stays 2^20", nfine_per_coarse(-2.7939677238464355e-06, 2 ** 20 * 64) == 2 ** 20))
    # the real hi_starframe slice: ch0 = 160,432,128, 2^20 channels
    lo, hi, cc, dc = segment_bounds(160432128, 1048576, 1033216)
    checks.append(("HI slice -> coarse 155, 749,568 channels from slice ch 0",
                   (lo, hi, cc) == (0, 749568, 155)))
    checks.append(("HI slice DC spike at slice channel 232,960", dc == 232960))
    d = np.ones((4, 1048576), np.float32)
    d[:, 232960] = 50.0
    seg, info = coarse_segment(d, 160432128, foff)
    checks.append(("DC spike repaired, segment shape right", seg.shape == (4, 749568) and seg[0, 232960] == 1.0))
    seg2, info2 = coarse_segment(np.ones((2, 1000), np.float32), 1033216 * 3 + 100, foff)
    checks.append(("a slice inside one coarse channel is kept whole", seg2.shape[1] == 1000 and not info2["straddled"]))
    print(describe(info, 1048576))
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    ok = all(c[1] for c in checks)
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    if sys.argv[1:2] == ["selftest"]:
        sys.exit(selftest())
    print(__doc__)
