#!/usr/bin/env python3
"""patch_turboseti.py - make turboSETI 2.3.x run on numpy >= 2.x.

turboSETI formats a 1-element numpy array with "%i" in a log line
(find_doppler.py, search_coarse_channel). numpy 2 refuses that, so every
search dies with "TypeError: %i format: a real number is required". The fix
is one line: format the array's first element. Idempotent; run it once after
`pip install turbo_seti`.

    python tools/patch_turboseti.py
"""
import importlib.util
import sys
from pathlib import Path

OLD = '" is: %i" % max_val.total_n_hits)'
NEW = '" is: %i" % int(max_val.total_n_hits[0]))'

spec = importlib.util.find_spec("turbo_seti")
if spec is None:
    sys.exit("turbo_seti is not installed in this python")
p = Path(spec.submodule_search_locations[0]) / "find_doppler" / "find_doppler.py"
s = p.read_text()
if NEW in s:
    print(f"already patched: {p}")
elif OLD in s:
    p.write_text(s.replace(OLD, NEW))
    print(f"patched: {p}")
else:
    sys.exit(f"pattern not found in {p} - turboSETI changed; check by hand")
