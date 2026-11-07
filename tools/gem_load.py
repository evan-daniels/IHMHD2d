"""Shared loaders and diagnostics for IHMHD2d (CentPack-format) runs.

CentPack layout: <run>/{b1,b2,b3,rho,u1,u2,u3,p}_files/<v>_pNNN_FFF, one file per
rank per frame, J rows (x) by K columns (y). The rank count is taken from
mesh_files/x_p*, so 5-rank and 10-rank runs load the same way.
"""
import glob
import os

import numpy as np

GAMMA = 5.0 / 3.0   # defaults for energies(); pass the run's own values if they differ
ETA = 0.005


class Run:
    def __init__(self, path):
        self.path = path
        xs = sorted(glob.glob(os.path.join(path, "mesh_files", "x_p*")))
        if not xs:
            raise FileNotFoundError(f"no mesh_files/x_p* in {path}")
        self.nranks = len(xs)
        self.x = np.concatenate([np.loadtxt(f).ravel() for f in xs])
        self.y = np.loadtxt(os.path.join(path, "mesh_files", "y")).ravel()
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        tf = sorted(glob.glob(os.path.join(path, "t_files", "t_*")))
        self.t = np.array([float(open(f).read().split()[0]) for f in tf])
        # the last write can repeat t_final to round-off; keep only increasing times
        n_ok = 1
        while n_ok < len(self.t) and self.t[n_ok] - self.t[n_ok - 1] > 1e-6:
            n_ok += 1
        self.t = self.t[:n_ok]
        # drop frames whose field files are incomplete (run still going or crashed mid-write)
        n = len(self.t)
        while n > 0 and not self._complete(n - 1):
            n -= 1
        self.t = self.t[:n]
        self._cache = {}

    def _file(self, var, r, n):
        return os.path.join(self.path, f"{var}_files", f"{var}_p{r:03d}_{n:03d}")

    def _complete(self, n):
        return all(os.path.getsize(self._file("b2", r, n)) > 0
                   if os.path.exists(self._file("b2", r, n)) else False
                   for r in range(self.nranks))

    def field(self, var, n):
        """Field on the full grid, shape (Nx, Ny)."""
        key = (var, n)
        if key not in self._cache:
            a = np.concatenate([np.loadtxt(self._file(var, r, n), ndmin=2)
                                for r in range(self.nranks)], axis=0)
            self._cache[key] = a
        return self._cache[key]

    def frame(self, t):
        return int(np.argmin(np.abs(self.t - t)))

    def has(self, t):
        return len(self.t) and abs(self.t[self.frame(t)] - t) < 1e-6

    def finite(self, n):
        return np.isfinite(self.field("b1", n)).all() and np.isfinite(self.field("b2", n)).all()

    # ---- derived quantities ----
    def jz(self, n):
        bx, by = self.field("b1", n), self.field("b2", n)
        # periodic in x, one-sided at the walls in y
        dbydx = (np.roll(by, -1, 0) - np.roll(by, 1, 0)) / (2 * self.dx)
        dbxdy = np.gradient(bx, self.dy, axis=1)
        return dbydx - dbxdy

    def psi(self, n):
        """Flux function: Bx = dpsi/dy, psi = 0 on the bottom wall."""
        bx = self.field("b1", n)
        return (np.cumsum(bx, axis=1) - 0.5 * bx) * self.dy

    def mid(self, a):
        """Average of the two rows straddling y = 0."""
        k0 = np.searchsorted(self.y, 0.0) - 1
        return 0.5 * (a[:, k0] + a[:, k0 + 1])

    def phi(self):
        """GEM reconnected flux: psi_max - psi_min along the midplane, minus its initial value."""
        p = np.array([phi_from_row(self.mid(self.field("b2", n)), self.x)
                      for n in range(len(self.t))])
        return p - p[0]

    def energy_log(self):
        """energy_log.txt written by the code (None for runs made before it existed)."""
        f = os.path.join(self.path, "energy_log.txt")
        if not os.path.exists(f):
            return None
        names = open(f).readline().lstrip("# ").split()
        a = np.loadtxt(f, ndmin=2)
        return {k: a[:, i] for i, k in enumerate(names)}

    def energies(self, n, gamma=GAMMA, eta=ETA):
        rho, p = self.field("rho", n), self.field("p", n)
        v2 = sum(self.field(c, n) ** 2 for c in ("u1", "u2", "u3"))
        b2 = sum(self.field(c, n) ** 2 for c in ("b1", "b2", "b3"))
        a = self.dx * self.dy
        return dict(kin=0.5 * np.sum(rho * v2) * a, th=np.sum(p) / (gamma - 1) * a,
                    mag=0.5 * np.sum(b2) * a, ohm=eta * np.sum(self.jz(n) ** 2) * a)


def phi_from_row(by, x):
    """psi_max - psi_min of psi = int By dx along a row (the original GEM definition).
    Unlike 0.5 int |By| dx it does not count the flux of a secondary island."""
    psi = np.concatenate([[0.0], np.cumsum(0.5 * (by[1:] + by[:-1]) * np.diff(x))])
    return psi.max() - psi.min()


def yoon_phi(path):
    """Same Phi for the lower sheet of a HallMHDcode double-GEM .h5 output (y_phys = y_file - 6.4)."""
    import h5py
    with h5py.File(path, "r") as f:
        t, x, y = f["time"][:], f["x"][:], f["y"][:] - 6.4
        by = f["fields/By"]                    # (time, y, x)
        k0 = np.searchsorted(y, 0.0) - 1
        p = np.array([phi_from_row(0.5 * (by[n, k0, :] + by[n, k0 + 1, :]), x)
                      for n in range(len(t))])
    return t, p - p[0]


def local_extrema(s, prom):
    """Indices of interior local maxima and minima of s that stand out by more than prom."""
    out = []
    for i in range(1, len(s) - 1):
        if (s[i] > s[i - 1] and s[i] >= s[i + 1]) or (s[i] < s[i - 1] and s[i] <= s[i + 1]):
            out.append(i)
    # drop pairs of adjacent extrema whose difference is below prom (grid noise)
    keep = list(out)
    changed = True
    while changed and len(keep) > 1:
        changed = False
        for a, b in zip(keep[:-1], keep[1:]):
            if abs(s[a] - s[b]) < prom:
                keep.remove(a); keep.remove(b); changed = True
                break
    return keep
