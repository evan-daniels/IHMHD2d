"""x-mirror symmetry check for GEM runs (problem 1).

GEM is mirror-symmetric about x = 0. Under x -> -x (B axial, v polar) the IC
(Bx = tanh(y) - cos(kx) sin(ky) terms, By ~ sin(kx) cos(ky), Bz = 0, v = 0) gives:
    even:  rho, p, Bx (b1), vy (u2), vz (u3)
    odd :  By (b2), Bz (b3), vx (u1)
The cell centres are symmetric about x = 0, so f[::-1, :] is the mirror image.
Error per field: max|f - s f[::-1]| / max|f| (s = +1 even, -1 odd); for fields that
are zero at t = 0 (b3, u*) the normalisation falls back to max|B| or 1.

Usage: python mirror_check.py RUN_DIR [--fields p,rho,b1,b2,b3] [--loc N]
  --loc N  prints the (x, y) of the largest error for the first N outputs.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gem_load import Run  # noqa: E402

PARITY = {"p": 1, "rho": 1, "b1": 1, "b2": -1, "b3": -1, "u1": -1, "u2": 1, "u3": 1}


def mirror_err(f, s, scale):
    d = np.abs(f - s * f[::-1, :])
    i = np.unravel_index(np.argmax(d), d.shape)
    return d[i] / scale, i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--fields", default="p,rho,b1,b2,b3")
    ap.add_argument("--loc", type=int, default=4)
    a = ap.parse_args()
    run = Run(a.run)
    fields = a.fields.split(",")
    xm = np.abs(run.x + run.x[::-1]).max()
    print(f"{a.run}: {len(run.x)}x{len(run.y)}, {run.nranks} ranks, "
          f"mesh asymmetry max|x + x[::-1]| = {xm:.1e}")
    print("   t   " + "".join(f"{v:>10s}" for v in fields) + "   loc(p): x, y")
    for n in range(len(run.t)):
        bmax = max(np.abs(run.field("b1", n)).max(), np.abs(run.field("b2", n)).max())
        row, loc = [], ""
        for v in fields:
            f = run.field(v, n)
            scale = np.abs(f).max()
            if v in ("b3",):
                scale = bmax
            elif v.startswith("u"):
                scale = scale if scale > 1e-12 else 1.0
            e, i = mirror_err(f, PARITY[v], scale)
            row.append(e)
            if v == fields[0] and n < a.loc and e > 0:
                loc = f"   {run.x[i[0]]:+7.3f}, {run.y[i[1]]:+6.3f}"
        print(f"{run.t[n]:6.2f} " + "".join(f"{e:10.1e}" for e in row) + loc)


if __name__ == "__main__":
    main()
