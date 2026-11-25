# IHMHD2d

IHMHD2d is a 2.5D Hall MHD solver with electron inertia, written for studying collisionless magnetic reconnection.
It is built on CentPack, the central-scheme package of Balbás and Tadmor. The magnetic field is advanced with
constrained transport, which keeps ∇·B = 0 to round-off. The code runs in parallel with MPI, with the domain split
along x.

## Physics

The electric field is given by the generalized Ohm's law

$$\mathbf{E} = -\mathbf{v}\times\mathbf{B}' + \frac{d_i}{n}\,\mathbf{J}\times\mathbf{B}' + \eta\,\mathbf{J} - \eta_H\nabla^2\mathbf{J}, \qquad \mathbf{B}' = \mathbf{B} - d_e^2\nabla^2\mathbf{B},$$

where d_i and d_e are the ion and electron inertial lengths, η is the resistivity and η_H the hyper-resistivity.
Electron inertia caps the whistler frequency at the electron cyclotron frequency ω_ce. The scheme conserves the
total energy exactly, including the electron kinetic energy ½d_e²|J|².

With d_e = 0 the code reduces to Hall MHD.

## Building

Requirements:

- an MPI C++ compiler (`mpicxx`) and `mpirun`,
- FFTW3, linked statically from `$(FFTW_DIR)/lib/libfftw3.a` (default `FFTW_DIR=/opt/homebrew/opt/fftw`).

Tested on macOS (Apple Silicon) with Apple clang 14 and MPICH.

```sh
make                                  # builds ./IHMHD2d
make FFTW_DIR=/path/to/fftw           # if FFTW is installed elsewhere
make test                             # unit tests and the d_e = 0 regression (NP=10 ranks by default)
```

## Running

The executable reads a file named `input` from the current directory and writes its output there. The output
folders must exist before the run starts:

```sh
mkdir -p runs/gem_1x && cd runs/gem_1x
cp ../../samples/input_gem_1x input
mkdir -p b1_files b2_files b3_files mesh_files p_files rho_files t_files u1_files u2_files u3_files
mpirun -np 10 ../../IHMHD2d > run.log 2>&1
```

Sample inputs in `samples/`:

| file | problem |
|---|---|
| `input_gem_1x` | GEM reconnection, d_e = 0.2, to t = 40 (10 ranks) |
| `input_gem_1x_de0` | GEM reconnection with d_e = 0, to t = 1 |
| `input_whistler` | linear whistler dispersion test, d_i = 0.5, d_e = 0.1 |
| `input_ot` | Orszag–Tang vortex, ideal MHD |

### Input file

The input file alternates between a label line and a value line. Only the order of the lines matters; the labels
are read and ignored.

| entry | meaning |
|---|---|
| `x_left`, `x_right`, `y_bottom`, `y_top` | domain |
| `J` | cells in x **per rank** (the total is J × number of ranks) |
| `K` | cells in y (total) |
| `L` | number of components (8) |
| `t_final`, `dt_out` | end time and output interval |
| `cfl`, `alpha` | CFL number and minmod3 limiter parameter |
| `par_count` | number of physics parameters that follow (4 for ideal runs, 10 otherwise) |
| `gamma` | ratio of specific heats |
| `problem` | 0 Orszag–Tang, 1 GEM reconnection, 2 linear whistler, 3 force-free resistive decay |
| `eta`, `di` | resistivity, ion inertial length |
| `s_B`, `s_Bz` | scale factors (default 1) for the upwind dissipation in the CT electric field: E_z (B_x, B_y update) and E_x, E_y (B_z update) |
| `eta_H` | hyper-resistivity |
| `c_whistler` | removed; read but ignored (keep 0) |
| `d_e` | electron inertial length |
| `c_whistler_inertia` | capped whistler upwinding fraction, used only when d_e > 0 (default 0.5) |

### Output

- `<var>_files/<var>_pNNN_FFF`: one file per rank (NNN) per output frame (FFF), for `rho`, `u1`–`u3`, `b1`–`b3` and
  `p`. Each file holds J rows (x) by K columns (y).
- `t_files/`: the output times.
- `mesh_files/`: the cell centres (`x_pNNN` per rank, `y`).
- `energy_log.txt`: energy budget at every output.

Setting the environment variable `EDRIFT=1` prints the total energy at every step to stdout. This is a diagnostic
only and doesn't change the solution.

`tools/gem_load.py` loads a run directory in Python. `tools/gem_figures.py` and `tools/mirror_check.py` make the
GEM figures and the x → −x symmetry check.

## Validation

The full numbers are in [docs/results/results.md](docs/results/results.md).

- **d_e = 0 regression** (`tests/regress_de0.sh`, `make regress`): GEM 1x to t = 1 on 10 ranks is bit-identical to
  the reference in `tests/reference/` (the CentPackL MHD2d Hall MHD code), including all field files,
  `energy_log.txt` and the time steps.
- **Unit tests** (`make test_helmholtz test_energy`): the Helmholtz solver with walls and the discrete energy
  identities hold to round-off on 1 and 4 ranks.
- **Resistive decay with electron inertia**: the decay-rate ratio converges to the theoretical 0.8 as N = 64, 128,
  256 (errors 1.7e-2, 1.9e-3, 1.7e-4). See [decay_ratio.png](docs/results/decay_ratio.png).
- **Whistler dispersion**: second-order convergence to the inertia-limited branch, with a maximum relative error of
  2.8e-3 at N = 400. See [whistler_dispersion.png](docs/results/whistler_dispersion.png).
- **GEM reconnection** (d_e = 0.2, η_H = 1e-4) at ½x, 1x and 2x resolution: total energy constant (drift 0), and
  similar reconnected flux and peak rates (0.260, 0.294, 0.276). See
  [gem_energy_rate.png](docs/results/gem_energy_rate.png) and [gem_flux_mirror.png](docs/results/gem_flux_mirror.png),
  and results.md for the x → −x symmetry discussion.

## Repository layout

```
Makefile
src/main.cc
src/core/        CentPack core: central scheme, reconstruction, time step, input, arrays
src/physics/     Hall MHD with electron inertia: fluxes, CT electric field, Helmholtz solver,
                 time evolution, initial and boundary conditions, output
include/         headers
samples/         example input files
tests/           unit tests, d_e = 0 regression script and its reference output
tools/           Python analysis scripts
docs/results/    verification results and figures
```

## Author

Evan Daniels, M.S. Physics, California State University, Northridge (advisor: Dr. Jorge Balbás).

## License

CentPack is copyright (C) Jorge Balbás and Eitan Tadmor and is distributed for research and instructional use only;
its notice is in [src/core/disclaimer.cc](src/core/disclaimer.cc) and is printed at the start of every run.
The license for the rest of this code is to be decided.

---

I'm not versed in coding standards and documentation conventions, so this README was written with help from Claude (Anthropic).
