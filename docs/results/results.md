# IHMHD2d verification results

All runs were made with the IHMHD2d binary at this commit's sources (built on macOS, MPICH, 10 ranks unless noted).
The run output is not in the repository; the numbers below are taken from it.

## Regression and unit tests (`make test`)

| test | result |
|---|---|
| `test_helmholtz_walls` (1 and 4 ranks) | worst relative residual 1.0e-13 (tolerance 1e-12) |
| `test_energy_identities` (1 and 4 ranks) | worst relative residual 1.7e-15 / 7.6e-16 (tolerance 1e-12) |
| d_e = 0 regression, GEM 1x to t = 1 | bit-identical to CentPackL MHD2d (e2626bb, c_whistler = 0): 92 field files, energy_log, time steps |
| d_e = 0 force-free decay, N064 | bit-identical to CentPackL MHD2d: 401 field files and energy_log |

## Resistive decay with electron inertia (`decay_ratio.png`)

Force-free field, k = 1, η = 0.01, d_i = 0. Amplitude decay rate fitted from the magnetic energy.
Theory: γ = ηk²/(1 + d_e²k²), so γ(d_e = 0.5)/γ(d_e = 0) = 0.8.

| N | γ(d_e = 0) | γ(d_e = 0.5) | ratio | error vs 0.8 |
|---|---|---|---|---|
| 64 | 0.010090 | 0.008212 | 0.81388 | +1.7e-2 |
| 128 | 0.010007 | 0.008020 | 0.80148 | +1.9e-3 |
| 256 | 0.010000 | 0.008002 | 0.80014 | +1.7e-4 |

Total energy (EDRIFT) is constant to the 13 printed digits in every d_e > 0 run.

## Whistler dispersion (`whistler_dispersion.png`)

d_i = 0.5, d_e = 0.1, c_whistler_inertia = 0, modes k = 1–12. Measured ω(k) (matrix pencil on FFT[B_y + iB_z])
against the R branch (1 + d_e²k²)ω² − d_i k²ω − k² = 0.

| N | max relative error (k = 12) | error at k = 6 |
|---|---|---|
| 200 | 1.2e-2 | 3.1e-3 |
| 400 | 2.8e-3 | 7.4e-4 |

Second-order convergence. At large k the branch approaches d_i B k²/(n(1 + d_e²k²)).

## GEM reconnection, d_e = 0.2, η_H = 1e-4, c_whistler_inertia = 0.5, to t = 40

Figures: `gem_energy_rate.png` (total energy; rate against ΔΦ), `gem_flux_mirror.png` (ΔΦ(t); u_y mirror error
under x → −x, with the uncapped runs for comparison).

| resolution | ΔΦ at t = 20 / 30 / 40 | peak rate | u_y mirror error at t = 30 / 40 | Etot drift | wall (10 ranks) |
|---|---|---|---|---|---|
| ½x | 1.832 / 3.022 / 2.730 | 0.260 at t = 17.5 | 4.0e-12 / 4.7e-10 | 0 | 410 s |
| 1x | 1.451 / 2.979 / 2.782 | 0.294 at t = 18.0 | 3.4e-11 / 4.6e-10 | 0 | 1772 s |
| 2x | pending | pending | pending | pending | pending |

Without the capped upwinding (c_whistler_inertia = 0) a grid-scale (λ_x ≈ 3.3Δx) Hall mode breaks the x → −x
symmetry from t ≈ 19; the u_y mirror error reaches 3.2e-4 (½x) and 1.6e-3 (1x) by t = 30.
With c = 0.5 the mirror error stays at round-off until the plasmoid forms (t ≈ 25) and then grows slowly
at physical wavelengths. c = 0.25 does not remove the mode (mirror defect 1.7e-4 by t = 22 at 1x).

The ½x and 1x runs are bit-identical to the development build they were ported from (all field files,
energy_log, EDRIFT).
