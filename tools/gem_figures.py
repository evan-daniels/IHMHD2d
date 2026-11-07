"""gem_figures.py -- two GEM summary figures for IHMHD2d runs (analysis only).

Figure 1 (--out1): total-energy drift from EDRIFT lines (left) and reconnection rate against
reconnected flux on the rising branch (right). Figure 2 (--out2): reconnected flux against t
(left) and u_y mirror error under x -> -x (right).

Every curve is LABEL=RUN_DIR. A run that is missing, has no EDRIFT lines, or has no output
frames is reported on stdout and skipped (never silently).

Example:
  python3 gem_figures.py \
      --energy "1x=runs/gem_1x" --rate "1/2x=runs/gem_halfx" "1x=runs/gem_1x" \
      --flux "1x=runs/gem_1x" --mirror "1x=runs/gem_1x" \
      --out1 fig1_energy_and_rate.png --out2 fig2_dphi_symmetry.png
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gem_load import Run  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def pairs(items):
    out = []
    for it in items or []:
        if "=" not in it:
            sys.exit(f"expected LABEL=RUN_DIR, got {it!r}")
        lab, path = it.split("=", 1)
        out.append((lab, path))
    return out


def skip(fig, lab, path, why):
    print(f"[{fig}] skipped {lab} ({path}): {why}")


def edrift(path):
    """(t, relative Etot drift) from the EDRIFT lines of run.log, or None."""
    f = os.path.join(path, "run.log")
    if not os.path.exists(f):
        return None
    t, e = [], []
    for line in open(f, errors="ignore"):
        if line.startswith("EDRIFT"):
            w = line.split()
            d = dict(zip(w[1::2], w[2::2]))
            t.append(float(d["t"])); e.append(float(d["Etot"]))
    if not t:
        return None
    t, e = np.array(t), np.array(e)
    return t, (e - e[0]) / abs(e[0])


def flux_rate(path):
    r = Run(path)
    if len(r.t) < 3:
        raise ValueError("fewer than 3 output frames")
    d = r.phi()
    return r.t, d, np.gradient(d, r.t)


def mirror_u2(path):
    r = Run(path)
    e = np.full(len(r.t), np.nan)
    for n in range(len(r.t)):
        u = r.field("u2", n)
        m = np.abs(u).max()
        if m > 0:
            e[n] = np.abs(u - u[::-1]).max() / m
    return r.t, e


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--energy", nargs="*", help="LABEL=RUN for the energy-drift panel")
    ap.add_argument("--rate", nargs="*", help="LABEL=RUN for the rate-vs-flux panel")
    ap.add_argument("--flux", nargs="*", help="LABEL=RUN for the flux-vs-t panel")
    ap.add_argument("--mirror", nargs="*", help="LABEL=RUN for the mirror-error panel")
    ap.add_argument("--title", default="GEM", help="title prefix")
    ap.add_argument("--out1", default="fig1_energy_and_rate.png")
    ap.add_argument("--out2", default="fig2_dphi_symmetry.png")
    a = ap.parse_args()
    plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for i, (lab, p) in enumerate(pairs(a.energy)):
        if not os.path.isdir(p):
            skip("fig1", lab, p, "no such folder"); continue
        r = edrift(p)
        if r is None:
            skip("fig1", lab, p, "no EDRIFT lines (run without EDRIFT=1)"); continue
        ax[0].plot(r[0], 100 * r[1], lw=1.8, label=lab, zorder=3 - 0.1 * i)
        print(f"[fig1] {lab}: max |Etot drift| {np.abs(r[1]).max():.1e}")
    ax[0].axhline(0, color="k", lw=0.6, zorder=1)
    ax[0].set(xlabel="t", ylabel=r"$(E_{tot}-E_{tot,0})/E_{tot,0}$  [%]", title=f"{a.title}: total energy")
    ax[0].legend()
    for lab, p in pairs(a.rate):
        try:
            t, d, r = flux_rate(p)
        except Exception as e:  # noqa: BLE001
            skip("fig1", lab, p, str(e)); continue
        m = t <= t[np.argmax(d)]
        idx = np.where(m)[0]
        if len(idx) > 4:
            m[idx[-2:]] = False                       # drop the turnover points of the centred derivative
        ax[1].plot(d[m], r[m], lw=1.8, label=lab)
        i = np.argmax(r); ax[1].plot(d[i], r[i], "o", color=ax[1].lines[-1].get_color())
    ax[1].set(xlabel=r"$\Delta\Phi$", ylabel=r"$d\Delta\Phi/dt$", title=f"{a.title}: rate vs reconnected flux")
    ax[1].legend()
    fig.tight_layout(); fig.savefig(a.out1, dpi=150); print("wrote", a.out1)

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for lab, p in pairs(a.flux):
        try:
            t, d, _ = flux_rate(p)
        except Exception as e:  # noqa: BLE001
            skip("fig2", lab, p, str(e)); continue
        ax[0].plot(t, d, lw=1.8, label=lab)
    ax[0].set(xlabel="t", ylabel=r"$\Delta\Phi$", title=f"{a.title}: reconnected flux")
    ax[0].legend(fontsize=9)
    for lab, p in pairs(a.mirror):
        try:
            t, e = mirror_u2(p)
        except Exception as ex:  # noqa: BLE001
            skip("fig2", lab, p, str(ex)); continue
        m = np.isfinite(e) & (e > 0)          # t = 0 is exactly symmetric (error 0): not drawn on a log axis
        ax[1].semilogy(t[m], e[m], lw=1.8, label=lab)
    ax[1].axhline(1e-4, color="gray", lw=0.6)
    ax[1].set(xlabel="t", ylabel=r"$u_y$ mirror error", ylim=(1e-16, 1), title=r"Departure from $x\to-x$ symmetry")
    ax[1].legend(fontsize=9)
    fig.tight_layout(); fig.savefig(a.out2, dpi=150); print("wrote", a.out2)


if __name__ == "__main__":
    main()
