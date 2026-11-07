////////////////////////////////////////////////////////////////////////////////
//
// NEW (2026-10-31): INERTIA2d copy of src/2d/SD2/evolution_2d_SD2.cc with electron inertia.
//   d_e = parameters(8) (default 0). For d_e > 0 the EM state is the generalized field
//   B' = (1 - d_e^2 lap_h) B (faces Bp1_xf, Bp2_yf; cell Bpz), advanced by the same CT /
//   SSP-RK3 sub-cycle with E* = -v x B' + (di/rho) J x B' + eta J + E_hyper (J from B);
//   B is recovered after every stage by helmholtz_solve_walls. Uniform-n inertia, ion
//   vorticity dropped (tools/design/electron_inertia.md). For d_e = 0 every operation is
//   the original one, so INERTIA2d at d_e = 0 is bit-identical to MHD2d.
//   This file defines CENTPACK::evolution_2d_SD2 and replaces SD2/evolution_2d_SD2.o in
//   the INERTIA2d build (Makefile target INERTIA2d_SD2).
//
// evolution_2d_SD2.cc — interleaved fluid / EM time step
//
// Fluid (rho, rho*v, energy) advances with the SSP-RK2 central scheme.
// B advances with a constrained-transport EM solver that sub-cycles on the
// whistler CFL. The EM sub-cycles are interleaved with the fluid RK stages:
//
//   Strang / Verlet splitting:
//   EM sub-cycle  0    -> dt/2   (rho, v frozen at t^n interface values)
//   fluid SSP-RK2 full dt        (B frozen at t^{n+1/2})
//   EM sub-cycle  dt/2 -> dt     (rho, v frozen at t^{n+1} interface values)
//   EM sub-cycles use SSP-RK3 (whistler-stable without upwind damping)
//
// un layout (8 components):
//   0 rho, 1-3 rho*v, 4-6 B (cell centered), 7 total energy
//   Fluid update touches 0-3 and 7 only. B is owned by the EM solver:
//   Bx, By live on faces (B1_xf, B2_yf) and are recovered to cell centers;
//   Bz lives in un(:,:,6) and is updated by a finite-volume EM flux.
//
// Face conventions:
//   B1_xf(j,k) : Bx on the east face of cell j   (x_{j+1/2}, y_k)
//   B2_yf(j,k) : By on the north face of cell k  (x_j, y_{k+1/2})
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"
#include <cstdlib>
#include <cstdio>

// NEW (2026-11-06): IHMHD2d. For d_e > 0 the step uses Formulation B (energy-conserving electron inertia,
//   see README.md): the fluid evolves e_f = 1/2 rho v^2 + p/(gamma-1) with no field source, the EM
//   sub-steps hand the exact field-energy change to e_f cell by cell, and the momentum gets the
//   electron-inertia stress -div T_e, T_e = d_e^2 (J J - 1/2 J^2 I), so the force is J x B'. The total
//   energy sum e_f + 1/2 <B, B'> is conserved to round-off (EDRIFT=1 prints it every step).
//   The EM upwind speed includes c_whistler_inertia d_i|B| khat/(rho (1 + d_e^2 khat^2)), khat = pi/h
//   (parameters(9), default 0.5), and the EM time step includes the same speed.
//   d_e = 0 runs the original MHD2d step unchanged (bit-identical to CentPackL MHD2d with c_whistler = 0).
//   Opt-in diagnostics (environment): EDRIFT=1, EFB_LOG=1, EFB_CHECK=N (README.md).

void CENTPACK::evolution_2d_SD2(doublearray3d& un, doublearray3d& B1_xf, doublearray3d& B2_yf,
    const doublearray1d& lambda, const doublearray1d& mu,
    const doublearray1d& dx_cell, const doublearray1d& dx_interface,
    const doublearray1d& dy_cell, const doublearray1d& dy_interface,
    const double& alpha, const doublearray1d& parameters, const double& dt_fluid,
    doublearray1d& diss, const int& id, const int& p)
{
    using namespace std;

    long j, k, l, J, K, L;
    J = un.getIndex1Size() - 4;
    K = un.getIndex2Size() - 4;
    L = un.getIndex3Size();

    const double eta    = parameters(2);
    const double di     = parameters(3);
    // upwind-dissipation scale factors for the EM flux (default 1 = unscaled)
    const long   n_par  = parameters.getIndex1Size();
    const double s_B    = (n_par > 4) ? parameters(4) : 1.0;   // Ez jumps -> Bx, By
    const double s_Bz   = (n_par > 5) ? parameters(5) : 1.0;   // Ey, -Ex jumps -> Bz
    const double eta_H  = (n_par > 6) ? parameters(6) : 0.0;   // hyperresistivity, E += -eta_H lap J
    // ORIGINAL (pre-2026-11-06):
    // const double c_whistler = (n_par > 7) ? parameters(7) : 0.0;   // whistler fraction in the EM upwinding speed
    // NEW (2026-11-06): the uncapped c_whistler option is removed (README.md); slot 7 is read but ignored
    const double c_whistler_removed = (n_par > 7) ? parameters(7) : 0.0;
    const bool   walls  = ((int) parameters(1) == 1);   // reconnection: conducting walls in y
    const int    yw_J   = walls ? 3 : 0;                  // ghost-fill modes for boundary_conditions_ct
    const int    yw_Bx  = walls ? 1 : 0;
    const int    yw_By  = walls ? 2 : 0;
    // NEW (2026-10-31): electron skin depth (electron inertia), parameters(8), default 0 = MHD2d
    const double de     = (n_par > 8) ? parameters(8) : 0.0;
    // NEW (2026-11-06): c_whistler_inertia, parameters(9), d_e > 0 only (default 0.5 when par_count < 10)
    const double c_wi   = (de > 0.0) ? ((n_par > 9) ? parameters(9) : 0.5) : 0.0;
    {
        static bool told = false;
        if (!told && id == 0)
        {
            if (c_whistler_removed != 0.0)
                printf("IHMHD2d: WARNING parameters(7) = c_whistler = %g is ignored (the uncapped c_whistler option was removed)\n", c_whistler_removed);
            if (de > 0.0)
                printf("IHMHD2d: c_whistler_inertia = %g%s\n", c_wi, (n_par > 9) ? "" : " (default; par_count < 10)");
            // NEW (2026-11-07): the EM step no longer includes the capped upwind speed
            if (de > 0.0 && c_wi > 0.5)
                printf("IHMHD2d: WARNING c_whistler_inertia = %g: values above the tested 0.5 may need a smaller EM step.\n", c_wi);
        }
        told = true;
    }
    // NEW (2026-11-06): Formulation B with the conservative fluid energy is the d_e > 0 scheme
    const bool   efb    = (de > 0.0);
    doublearray3d* efb_edge = 0;   // edge E* export from C_flux_em, set during Formulation B sub-steps
    const int    hz_x   = walls ? HZ_NEU : HZ_PER;           // Helmholtz kind per component
    const int    hz_y   = walls ? HZ_DIR : HZ_PER;
    const int    hz_z   = walls ? HZ_NEU : HZ_PER;
    // const double cfl_em = 0.4;
    const double cfl_em = .2;

    // Fluid components updated by the central scheme (B = 4,5,6 excluded)
    const int n_fluid = 5;
    const int fluid_comp[n_fluid] = {0, 1, 2, 3, 7};

    // ---------- fluid arrays ----------
    doublearray3d u_N(J+4,K+4,L), u_S(J+4,K+4,L), u_E(J+4,K+4,L), u_W(J+4,K+4,L);
    doublearray3d C0(J+4,K+4,L), C1(J+4,K+4,L);

    // ---------- frozen fluid state at interfaces ----------
    doublearray2d rho_N(J+4,K+4), rho_S(J+4,K+4), rho_E(J+4,K+4), rho_W(J+4,K+4);
    doublearray3d v_N(J+4,K+4,3), v_S(J+4,K+4,3), v_E(J+4,K+4,3), v_W(J+4,K+4,3);

    // ---------- EM arrays ----------
    doublearray3d Jc(J+4,K+4,3);
    doublearray3d B_N(J+4,K+4,3), B_S(J+4,K+4,3), B_E(J+4,K+4,3), B_W(J+4,K+4,3);
    doublearray3d J_N(J+4,K+4,3), J_S(J+4,K+4,3), J_E(J+4,K+4,3), J_W(J+4,K+4,3);
    doublearray3d E_N(J+4,K+4,3), E_S(J+4,K+4,3), E_E(J+4,K+4,3), E_W(J+4,K+4,3);
    doublearray2d Cx_em(J+4,K+4), Cy_em(J+4,K+4), Cz_em(J+4,K+4);
    doublearray1d lambda_em(J+4), mu_em(K+4);
    doublearray1d u_cell(L);
    
    doublearray2d Bx0(J+4,K+4), By0(J+4,K+4), Bz0(J+4,K+4);   // RK3 base state
    doublearray3d E_hyper(J+4,K+4,3);                           // -eta_H lap J on the CT lattice

    double t_em = 0.0;
    double dt_em = 0.0;

    // NEW (2026-10-31): electron inertia state, persistent across calls (B is frozen during the
    // fluid step, so B' carries over). Built from B on the first call with d_e > 0.
    static doublearray3d *Bp1 = 0, *Bp2 = 0, *Bpz = 0, *Bz_tmp = 0, *un_p = 0;
    static doublearray3d *Bp_N = 0, *Bp_S = 0, *Bp_E = 0, *Bp_W = 0, *Jd_N = 0, *Jd_S = 0, *Jd_E = 0, *Jd_W = 0;
    static double hx = 0.0, hy = 0.0;
    if (de > 0.0 && Bp1 == 0)
    {
        // uniform mesh required by the FFT solve
        hx = dx_cell(2); hy = dy_cell(2);
        double dev = 0.0;
        for (j = 2; j < J+2; j++) dev = std::max(dev, std::fabs(dx_cell(j) - hx));
        for (k = 2; k < K+2; k++) dev = std::max(dev, std::fabs(dy_cell(k) - hy));
        if (dev > 1e-12*std::min(hx, hy))
        {
            if (id == 0) printf("INERTIA2d: d_e > 0 needs a uniform mesh\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        Bp1 = new doublearray3d(J+4,K+4,1);  Bp2 = new doublearray3d(J+4,K+4,1);
        Bpz = new doublearray3d(J+4,K+4,1);  Bz_tmp = new doublearray3d(J+4,K+4,1);
        un_p = new doublearray3d(J+4,K+4,L);
        Bp_N = new doublearray3d(J+4,K+4,3); Bp_S = new doublearray3d(J+4,K+4,3);
        Bp_E = new doublearray3d(J+4,K+4,3); Bp_W = new doublearray3d(J+4,K+4,3);
        Jd_N = new doublearray3d(J+4,K+4,3); Jd_S = new doublearray3d(J+4,K+4,3);
        Jd_E = new doublearray3d(J+4,K+4,3); Jd_W = new doublearray3d(J+4,K+4,3);
        for (j = 0; j < J+4; j++)
            for (k = 0; k < K+4; k++)
            {
                (*Bp1)(j,k,0) = (*Bp2)(j,k,0) = (*Bpz)(j,k,0) = 0.0;
                (*Bz_tmp)(j,k,0) = un(j,k,6);
            }
        boundary_conditions_ct(B1_xf, id, p, yw_Bx);
        boundary_conditions_ct(B2_yf, id, p, yw_By);
        boundary_conditions_ct(*Bz_tmp, id, p, yw_Bx);       // Bz: cell-centred, even at walls
        helmholtz_apply_walls(B1_xf, 0, *Bp1, 0, hz_x, de, hx, hy);
        helmholtz_apply_walls(B2_yf, 0, *Bp2, 0, hz_y, de, hx, hy);
        helmholtz_apply_walls(*Bz_tmp, 0, *Bpz, 0, hz_z, de, hx, hy);
        boundary_conditions_ct(*Bp1, id, p, yw_Bx);
        boundary_conditions_ct(*Bp2, id, p, yw_By);
        boundary_conditions_ct(*Bpz, id, p, yw_Bx);
        if (id == 0) printf("INERTIA2d: electron inertia on, d_e = %g (d_e/dx = %.2f)\n", de, de/hx);
    }

    // ------------------------------------------------------------------------
    // J = curl B on the interior; ghosts filled by y-boundary (periodic or wall) / MPI-x exchange
    // ------------------------------------------------------------------------
    auto compute_J = [&]()
    {
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            {
                double dx = dx_cell(j), dy = dy_cell(k);
                Jc(j,k,0) =  (un(j,k+1,6) - un(j,k-1,6))/(2.0*dy);
                Jc(j,k,1) = -(un(j+1,k,6) - un(j-1,k,6))/(2.0*dx);
                Jc(j,k,2) =  (un(j+1,k,5) - un(j-1,k,5))/(2.0*dx)
                           - (un(j,k+1,4) - un(j,k-1,4))/(2.0*dy);
            }
        boundary_conditions_ct(Jc, id, p, yw_J);
    };

    // ------------------------------------------------------------------------
    // Pull rho and v from the fluid reconstruction; frozen during a sub-cycle
    // ------------------------------------------------------------------------
    auto freeze_fluid = [&]()
    {
        for (k = 1; k < K+3; k++)
            for (j = 1; j < J+3; j++)
            {
                rho_N(j,k) = u_N(j,k,0);  rho_S(j,k) = u_S(j,k,0);
                rho_E(j,k) = u_E(j,k,0);  rho_W(j,k) = u_W(j,k,0);
                for (l = 0; l < 3; l++)
                {
                    v_N(j,k,l) = u_N(j,k,l+1)/u_N(j,k,0);
                    v_S(j,k,l) = u_S(j,k,l+1)/u_S(j,k,0);
                    v_E(j,k,l) = u_E(j,k,l+1)/u_E(j,k,0);
                    v_W(j,k,l) = u_W(j,k,l+1)/u_W(j,k,0);
                }
            }
    };

    // ------------------------------------------------------------------------
    // EM right-hand side: Cx_em, Cy_em, Cz_em from the current B and frozen rho, v
    // (lambda_em, mu_em must already hold dt_em/dx, dt_em/dy)
    // w = the stage's SSP-RK3 weight; diss accumulates w*dt_em times the stage's
    //   (0) Ohmic   sum eta |Jc|^2 dA
    //   (1) Hyper   sum eta_H J.(-lap J) dA  (= eta_H sum |grad_h J|^2 dA)
    //   (2) Lorentz sum v.(Jc x B) dA        (magnetic -> kinetic work)
    // over owned cells (rank-local; reduced at output). Diagnostics only.
    // ------------------------------------------------------------------------
    auto em_rhs = [&](double w)
    {
        compute_J();
        double ohm = 0.0, lor = 0.0, hyp = 0.0;
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            {
                double dA = dx_cell(j)*dy_cell(k);
                double jx = Jc(j,k,0), jy = Jc(j,k,1), jz = Jc(j,k,2);
                double bx = un(j,k,4), by = un(j,k,5), bz = un(j,k,6);
                double r  = un(j,k,0);
                ohm += eta*(jx*jx + jy*jy + jz*jz)*dA;
                lor += (un(j,k,1)*(jy*bz - jz*by) + un(j,k,2)*(jz*bx - jx*bz)
                      + un(j,k,3)*(jx*by - jy*bx))/r*dA;
            }
        if (eta_H > 0.0)
            hyp = hyper_em(B1_xf, B2_yf, un, dx_interface, dy_interface, eta_H, E_hyper);
        diss(0) += w*dt_em*ohm;
        diss(1) += w*dt_em*hyp;
        diss(2) += w*dt_em*lor;

        // NEW (2026-10-31): d_e = 0 keeps the original calls unchanged (bit-identical to MHD2d)
        if (de == 0.0)
        {
        reconstruction_em(un, Jc, B1_xf, B2_yf, B_N, B_S, B_E, B_W,
                          J_N, J_S, J_E, J_W,
                          dx_cell, dx_interface, dy_cell, dy_interface, alpha);
        electric_field(rho_N, rho_S, rho_E, rho_W, v_N, v_S, v_E, v_W,
                       B_N, B_S, B_E, B_W, J_N, J_S, J_E, J_W,
                       E_N, E_S, E_E, E_W, di, eta);
        C_flux_em(rho_N, rho_S, rho_E, rho_W, B_N, B_S, B_E, B_W,
                  E_N, E_S, E_E, E_W, lambda_em, mu_em,
                  dx_interface, dy_interface, di, Cx_em, Cy_em, Cz_em, walls ? 1 : 0, s_B, s_Bz,
                  // ORIGINAL (pre-2026-11-06):
                  // eta_H > 0.0 ? &E_hyper : 0, c_whistler);
                  // NEW (2026-11-06): c_whistler removed (0 = the value of every CentPackL reference and paper run)
                  eta_H > 0.0 ? &E_hyper : 0, 0.0);
        }
        else
        {
            // NEW (2026-10-31): E* = -v x B' + (di/rho) J x B' + eta J (+ E_hyper), J from B.
            // B' interface states from the same reconstruction applied to a copy of un whose
            // cell field (4..6) is B'; J states come from Jc (curl of B) unchanged. The LxF
            // upwinding in C_flux_em then acts on the jumps of the evolved variable B'.
            doublearray3d& up = *un_p;
            for (k = 0; k < K+4; k++)
                for (j = 0; j < J+4; j++)
                    for (l = 0; l < L; l++) up(j,k,l) = un(j,k,l);
            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                {
                    up(j,k,4) = 0.5*((*Bp1)(j,k,0) + (*Bp1)(j-1,k,0));
                    up(j,k,5) = 0.5*((*Bp2)(j,k,0) + (*Bp2)(j,k-1,0));
                    up(j,k,6) = (*Bpz)(j,k,0);
                }
            boundary_conditions(up, parameters, id, p);
            reconstruction_em(up, Jc, *Bp1, *Bp2, *Bp_N, *Bp_S, *Bp_E, *Bp_W,
                              *Jd_N, *Jd_S, *Jd_E, *Jd_W,
                              dx_cell, dx_interface, dy_cell, dy_interface, alpha);
            electric_field(rho_N, rho_S, rho_E, rho_W, v_N, v_S, v_E, v_W,
                           *Bp_N, *Bp_S, *Bp_E, *Bp_W, *Jd_N, *Jd_S, *Jd_E, *Jd_W,
                           E_N, E_S, E_E, E_W, di, eta);
            C_flux_em(rho_N, rho_S, rho_E, rho_W, *Bp_N, *Bp_S, *Bp_E, *Bp_W,
                      E_N, E_S, E_E, E_W, lambda_em, mu_em,
                      dx_interface, dy_interface, di, Cx_em, Cy_em, Cz_em, walls ? 1 : 0, s_B, s_Bz,
                      // ORIGINAL (pre-2026-11-06):
                      // eta_H > 0.0 ? &E_hyper : 0, c_whistler);
                      // NEW (2026-11-06): capped whistler upwind (c_whistler_inertia); edge E* exported for Formulation B
                      eta_H > 0.0 ? &E_hyper : 0, 0.0, c_wi, de, efb_edge);
        }
    };

    // ------------------------------------------------------------------------
    // SSP-RK stage:  B <- a0 * B0 + a1 * (B + C)     (each stage is a CT update,
    // so div B is preserved exactly), then ghosts and cell-centered Bx, By
    // ------------------------------------------------------------------------
    auto em_stage = [&](double a0, double a1)
    {
        // NEW (2026-10-31): d_e > 0 -- the CT stage acts on B', then B = (1 - d_e^2 lap_h)^-1 B'
        // per component (div B' = 0 by CT => div B = 0, design doc section 2.3)
        if (de > 0.0)
        {
            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                {
                    (*Bp1)(j,k,0) = a0*Bx0(j,k) + a1*((*Bp1)(j,k,0) + Cx_em(j,k));
                    (*Bp2)(j,k,0) = a0*By0(j,k) + a1*((*Bp2)(j,k,0) + Cy_em(j,k));
                    (*Bpz)(j,k,0) = a0*Bz0(j,k) + a1*((*Bpz)(j,k,0) + Cz_em(j,k));
                }
            boundary_conditions_ct(*Bp1, id, p, yw_Bx);
            boundary_conditions_ct(*Bp2, id, p, yw_By);
            boundary_conditions_ct(*Bpz, id, p, yw_Bx);

            helmholtz_solve_walls(*Bp1, 0, B1_xf, 0, hz_x, de, hx, hy, id, p);
            helmholtz_solve_walls(*Bp2, 0, B2_yf, 0, hz_y, de, hx, hy, id, p);
            helmholtz_solve_walls(*Bpz, 0, *Bz_tmp, 0, hz_z, de, hx, hy, id, p);
            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                    un(j,k,6) = (*Bz_tmp)(j,k,0);

            boundary_conditions_ct(B1_xf, id, p, yw_Bx);
            boundary_conditions_ct(B2_yf, id, p, yw_By);

            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                {
                    un(j,k,4) = 0.5*(B1_xf(j,k,0) + B1_xf(j-1,k,0));
                    un(j,k,5) = 0.5*(B2_yf(j,k,0) + B2_yf(j,k-1,0));
                }

            boundary_conditions(un, parameters, id, p);
            return;
        }

        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            {
                B1_xf(j,k,0) = a0*Bx0(j,k) + a1*(B1_xf(j,k,0) + Cx_em(j,k));
                B2_yf(j,k,0) = a0*By0(j,k) + a1*(B2_yf(j,k,0) + Cy_em(j,k));
                un(j,k,6)    = a0*Bz0(j,k) + a1*(un(j,k,6)    + Cz_em(j,k));
            }

        boundary_conditions_ct(B1_xf, id, p, yw_Bx);
        boundary_conditions_ct(B2_yf, id, p, yw_By);

        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            {
                un(j,k,4) = 0.5*(B1_xf(j,k,0) + B1_xf(j-1,k,0));
                un(j,k,5) = 0.5*(B2_yf(j,k,0) + B2_yf(j,k-1,0));
            }

        boundary_conditions(un, parameters, id, p);
    };

    // ------------------------------------------------------------------------
    // NEW (2026-11-06): Formulation B (d_e > 0).
    //   EM sub-step: total energy E = e_f + w, w per cell = 1/2 B.B' (face products split half-and-half
    //     to the cells); E^{n+1} = E^n - dt_em sum_i b_i A[div (S_i + 1/2 d_e^2 G_i)], S = E* x B from the exact
    //     edge E* of C_flux_em, G = (bdot grad b - b grad bdot) on each B lattice, bdot_i from the Shu-Osher
    //     stage algebra; so e_f += -dt_em sum_i b_i A[div F_i] - (w^{n+1} - w^n) per cell (exact by construction).
    //   Fluid step: e_f has no source; momentum += -div T_e, T_e = d_e^2 (J J - 1/2 J^2 I), compact J averaged
    //     to the fluid faces, so the force is J x B'.
    //   EFB_CHECK=N (environment): for the first N EM sub-steps print the per-cell exchange check
    //     Q - [dt sum b_i A(J_i.E_i) + R_RK], R_RK = -(dw - dt sum b_i wdot_i), wdot_i = A[1/2(bdot.B'_i + B_i.Bdot'_i)].
    //   EFB_LOG=1 (environment): per fluid step, an EFB line with sum Q, sum dt b J.E, sum R_RK (x dA).
    // ------------------------------------------------------------------------
    static const long efb_check = getenv("EFB_CHECK") ? atol(getenv("EFB_CHECK")) : 0;
    static const bool efb_log = getenv("EFB_LOG") != 0 || efb_check > 0;
    static long efb_nsub = 0;
    static double efb_sQ = 0.0, efb_sJE = 0.0, efb_sR = 0.0;      // rank-local, this fluid step
    static doublearray3d *eE = 0, *eB0 = 0, *eBs = 0, *eBps = 0, *eBd = 0, *eCd = 0, *eT = 0, *eZ = 0;
    static doublearray2d *eDiv = 0, *eW0 = 0, *eQ = 0, *eJE = 0, *eWd = 0;
    if (efb && eE == 0)
    {
        doublearray3d** a3[8] = {&eE, &eB0, &eBs, &eBps, &eBd, &eCd, &eT, &eZ};
        for (int q = 0; q < 8; q++) *a3[q] = new doublearray3d(J+4, K+4, q == 7 ? 1 : 3);
        doublearray2d** a2[5] = {&eDiv, &eW0, &eQ, &eJE, &eWd};
        for (int q = 0; q < 5; q++) { *a2[q] = new doublearray2d(J+4, K+4);
            for (k = 0; k < K+4; k++) for (j = 0; j < J+4; j++) (**a2[q])(j,k) = 0.0; }
        for (int q = 0; q < 7; q++) for (k = 0; k < K+4; k++) for (j = 0; j < J+4; j++) for (l = 0; l < 3; l++) (**a3[q])(j,k,l) = 0.0;
        if (id == 0) printf("IHMHD2d: Formulation B on (d_e = %g), EFB_CHECK %ld, EFB_LOG %d\n", de, efb_check, (int) efb_log);
    }
    const int efb_yw[3] = {yw_Bx, yw_By, yw_Bx};
    // CT ghosts for each component of a 3-component lattice array (Bx: x-faces, By: y-faces, Bz: cells)
    auto efb_ct3 = [&](doublearray3d& T) {
        for (int c = 0; c < 3; c++)
        {
            for (k = 0; k < K+4; k++) for (j = 0; j < J+4; j++) (*eZ)(j,k,0) = T(j,k,c);
            boundary_conditions_ct(*eZ, id, p, efb_yw[c]);
            for (k = 0; k < K+4; k++) for (j = 0; j < J+4; j++) T(j,k,c) = (*eZ)(j,k,0);
        }
    };
    // snapshot of B (prime = false) or B' (prime = true; = B at d_e = 0) on the CT lattices
    auto efb_snap = [&](doublearray3d& T, bool prime) {
        bool pr = prime && de > 0.0;
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            {
                T(j,k,0) = pr ? (*Bp1)(j,k,0) : B1_xf(j,k,0);
                T(j,k,1) = pr ? (*Bp2)(j,k,0) : B2_yf(j,k,0);
                T(j,k,2) = pr ? (*Bpz)(j,k,0) : un(j,k,6);
            }
        efb_ct3(T);
    };
    // w per cell, (3): faces split half-and-half
    auto efb_w = [&](const doublearray3d& B, const doublearray3d& P, long a, long b) {
        return 0.25*(B(a-1,b,0)*P(a-1,b,0) + B(a,b,0)*P(a,b,0))
             + 0.25*(B(a,b-1,1)*P(a,b-1,1) + B(a,b,1)*P(a,b,1)) + 0.5*B(a,b,2)*P(a,b,2); };
    // A[div S] per cell from B (3-comp, CT ghosts) and the edge E* (0 Ez nodes, 1 Ey x-faces, 2 Ex y-faces)
    auto efb_divS = [&](const doublearray3d& B, const doublearray3d& E, long a, long b) {
        const double ddx = hx > 0.0 ? hx : dx_cell(2), ddy = hy > 0.0 ? hy : dy_cell(2);
        auto F1y = [&](long jj, long kk) { return  E(jj,kk,0)*0.5*(B(jj,kk,0) + B(jj,kk+1,0)); };
        auto F2x = [&](long jj, long kk) { return -E(jj,kk,0)*0.5*(B(jj,kk,1) + B(jj+1,kk,1)); };
        auto F3x = [&](long jj, long kk) { return  E(jj,kk,1)*0.5*(B(jj,kk,2) + B(jj+1,kk,2)); };
        auto F3y = [&](long jj, long kk) { return -E(jj,kk,2)*0.5*(B(jj,kk,2) + B(jj,kk+1,2)); };
        auto dX  = [&](long jj, long kk) { return (F1y(jj,kk) - F1y(jj,kk-1))/ddy; };
        auto dY  = [&](long jj, long kk) { return (F2x(jj,kk) - F2x(jj-1,kk))/ddx; };
        return 0.5*(dX(a-1,b) + dX(a,b)) + 0.5*(dY(a,b-1) + dY(a,b))
             + (F3x(a,b) - F3x(a-1,b))/ddx + (F3y(a,b) - F3y(a,b-1))/ddy; };
    // A[J.E] per cell, compact J = C^T B on the E lattices
    auto efb_JE = [&](const doublearray3d& B, const doublearray3d& E, long a, long b) {
        const double ddx = hx > 0.0 ? hx : dx_cell(2), ddy = hy > 0.0 ? hy : dy_cell(2);
        auto qa = [&](long jj, long kk) { return -(B(jj,kk+1,0) - B(jj,kk,0))/ddy*E(jj,kk,0); };
        auto qb = [&](long jj, long kk) { return  (B(jj+1,kk,1) - B(jj,kk,1))/ddx*E(jj,kk,0); };
        auto qc = [&](long jj, long kk) { return -(B(jj+1,kk,2) - B(jj,kk,2))/ddx*E(jj,kk,1); };
        auto qd = [&](long jj, long kk) { return  (B(jj,kk+1,2) - B(jj,kk,2))/ddy*E(jj,kk,2); };
        auto jX = [&](long jj, long kk) { return 0.5*(qa(jj,kk) + qa(jj,kk-1)); };
        auto jY = [&](long jj, long kk) { return 0.5*(qb(jj,kk) + qb(jj-1,kk)); };
        return 0.5*(jX(a-1,b) + jX(a,b)) + 0.5*(jY(a,b-1) + jY(a,b))
             + 0.5*(qc(a,b) + qc(a-1,b)) + 0.5*(qd(a,b) + qd(a,b-1)); };
    // A[div G(f, g)] per cell, G = (f grad g - g grad f) on each lattice (f = bdot, g = b)
    auto efb_divG = [&](const doublearray3d& Fd, const doublearray3d& G, long a, long b) {
        const double ddx = hx > 0.0 ? hx : dx_cell(2), ddy = hy > 0.0 ? hy : dy_cell(2);
        auto gx = [&](int c, long jj, long kk) { return (Fd(jj,kk,c)*G(jj+1,kk,c) - G(jj,kk,c)*Fd(jj+1,kk,c))/ddx; };
        auto gy = [&](int c, long jj, long kk) { return (Fd(jj,kk,c)*G(jj,kk+1,c) - G(jj,kk,c)*Fd(jj,kk+1,c))/ddy; };
        auto dl = [&](int c, long jj, long kk) { return (gx(c,jj,kk) - gx(c,jj-1,kk))/ddx + (gy(c,jj,kk) - gy(c,jj,kk-1))/ddy; };
        return 0.5*(dl(0,a-1,b) + dl(0,a,b)) + 0.5*(dl(1,a,b-1) + dl(1,a,b)) + dl(2,a,b); };
    // direct quadratic-form rate A[1/2 (bdot.B' + B.Bdot')] per cell
    auto efb_wd = [&](const doublearray3d& Bd, const doublearray3d& P, const doublearray3d& B, const doublearray3d& Cd, long a, long b) {
        auto v = [&](int c, long jj, long kk) { return 0.5*(Bd(jj,kk,c)*P(jj,kk,c) + B(jj,kk,c)*Cd(jj,kk,c)); };
        return 0.5*(v(0,a-1,b) + v(0,a,b)) + 0.5*(v(1,a,b-1) + v(1,a,b)) + v(2,a,b); };

    // one SSP-RK3 EM sub-step with the Formulation B energy exchange (replaces the three stage lines)
    auto efb_substep = [&]() {
        const double A0[3] = {0.0, 0.75, 1.0/3.0}, A1[3] = {1.0, 0.25, 2.0/3.0}, BW[3] = {1.0/6.0, 1.0/6.0, 2.0/3.0};
        doublearray3d& Pn = *eBps;
        efb_snap(*eB0, false);
        efb_snap(Pn, true);
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            { (*eW0)(j,k) = efb_w(*eB0, Pn, j, k); (*eDiv)(j,k) = 0.0; (*eJE)(j,k) = 0.0; (*eWd)(j,k) = 0.0; }
        for (int s = 0; s < 3; s++)
        {
            efb_snap(*eBs, false);
            if (efb_log) efb_snap(Pn, true);
            efb_edge = eE;
            em_rhs(BW[s]);
            efb_edge = 0;
            if (efb_log)
            {
                for (k = 2; k < K+2; k++)
                    for (j = 2; j < J+2; j++)
                    { (*eCd)(j,k,0) = Cx_em(j,k)/dt_em; (*eCd)(j,k,1) = Cy_em(j,k)/dt_em; (*eCd)(j,k,2) = Cz_em(j,k)/dt_em; }
                efb_ct3(*eCd);
            }
            em_stage(A0[s], A1[s]);
            efb_snap(*eT, false);
            for (k = 0; k < K+4; k++)
                for (j = 0; j < J+4; j++)
                    for (l = 0; l < 3; l++)
                        (*eBd)(j,k,l) = ((*eT)(j,k,l) - A0[s]*(*eB0)(j,k,l) - A1[s]*(*eBs)(j,k,l))/(A1[s]*dt_em);
            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                {
                    double d = efb_divS(*eBs, *eE, j, k);
                    if (de > 0.0) d += 0.5*de*de*efb_divG(*eBd, *eBs, j, k);
                    (*eDiv)(j,k) += BW[s]*d;
                    if (efb_log)
                    {
                        (*eJE)(j,k) += BW[s]*efb_JE(*eBs, *eE, j, k);
                        (*eWd)(j,k) += BW[s]*efb_wd(*eBd, Pn, *eBs, *eCd, j, k);
                    }
                }
        }
        // end of sub-step: exchange to e_f, per cell
        efb_snap(*eT, false);
        efb_snap(Pn, true);
        double rmax = 0.0, qmax = 0.0, rrmax = 0.0, sE = 0.0, sEa = 0.0, sQ = 0.0, sJE = 0.0, sR = 0.0;
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            {
                double dA = dx_cell(j)*dy_cell(k);
                double dw = efb_w(*eT, Pn, j, k) - (*eW0)(j,k);
                double dE = -dt_em*(*eDiv)(j,k);
                double Q  = dE - dw;
                (*eQ)(j,k) += Q;
                sE += dE*dA; sEa += fabs(dE)*dA; sQ += Q*dA;
                if (efb_log)
                {
                    double R = -(dw - dt_em*(*eWd)(j,k));
                    double res = Q - (dt_em*(*eJE)(j,k) + R);
                    rmax = std::max(rmax, fabs(res)); qmax = std::max(qmax, fabs(Q)); rrmax = std::max(rrmax, fabs(R));
                    sJE += dt_em*(*eJE)(j,k)*dA; sR += R*dA;
                }
            }
        efb_sQ += sQ; efb_sJE += sJE; efb_sR += sR;
        if (efb_nsub < efb_check)
        {
            double a[3] = {rmax, qmax, rrmax}, ga[3], b[5] = {sE, sEa, sQ, sJE, sR}, gb[5];
            MPI_Allreduce(a, ga, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(b, gb, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (id == 0)
                printf("EFBCHK sub %ld dt_em %.3e  cell: max|Q - dt bJE - R_RK|/max|Q| %.2e  max|R_RK|/max|Q| %.2e  "
                       "global: sum dE/sum|dE| %.2e  sum Q %+.6e  sum dt bJE %+.6e  sum R_RK %+.6e\n",
                       efb_nsub, dt_em, ga[0]/ga[1], ga[2]/ga[1], gb[0]/gb[1], gb[2], gb[3], gb[4]);
        }
        efb_nsub++;
    };

    // ------------------------------------------------------------------------
    // EM sub-cycle: SSP-RK3 steps from the current t_em up to t_target
    // ------------------------------------------------------------------------
    auto em_subcycle = [&](double t_target)
    {
        int nsub = 0;
        while (t_em < t_target)
        {
            // --- whistler CFL, per cell, global min ---
            double dtp_em = 1.0e30;
            double rx, ry;
            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                {
                    for (l = 0; l < L; l++) u_cell(l) = un(j,k,l);
                    spectral_radii_em(u_cell, dx_cell(j), dy_cell(k), di, rx, ry);
                    double rate = rx/dx_cell(j) + ry/dy_cell(k);
                    // NEW (2026-10-31): whistler capped at omega_ce by electron inertia:
                    // omega(k) = di|B| k^2 / (rho (1 + de^2 k^2)) at k^2 = khat^2_max = 4/dx^2 + 4/dy^2
                    if (de > 0.0)
                        rate /= 1.0 + de*de*(4.0/(dx_cell(j)*dx_cell(j)) + 4.0/(dy_cell(k)*dy_cell(k)));
                    // NEW (2026-11-06): + the c_whistler_inertia upwind speed c d_i|B| khat/(rho (1 + d_e^2 khat^2)) / h
                    // per direction, khat = pi/h (d_e > 0 only), at the same cfl_em
                    // ORIGINAL (pre-2026-11-07):
                    // if (de > 0.0 && c_wi > 0.0)
                    // {
                    //     double Bm = sqrt(u_cell(4)*u_cell(4) + u_cell(5)*u_cell(5) + u_cell(6)*u_cell(6)), rr = std::max(u_cell(0), 0.1);
                    //     for (int q = 0; q < 2; q++)
                    //     {
                    //         double h = q ? dy_cell(k) : dx_cell(j);
                    //         double kh = M_PI/h;
                    //         rate += c_wi*di*Bm*kh/(rr*(1.0 + de*de*kh*kh))/h;
                    //     }
                    // }
                    // NEW (2026-11-07): term removed from the EM step (no effect on results at 1/2x, 1x, 2x
                    // for c = 0.5; doubled wall time). The upwind speed stays in C_flux_em.
                    if (rate > 0.0)
                        dtp_em = std::min(dtp_em, cfl_em/rate);
                    // parabolic limit of the explicit 4th-order term (SSP-RK3 real-axis
                    // stability ~2.5; max eigenvalue of lap^2 on the 5-point lattice 64/h^4)
                    if (eta_H > 0.0)
                    {
                        double h = std::min(dx_cell(j), dy_cell(k));
                        // ORIGINAL (pre-2026-11-06):
                        // dtp_em = std::min(dtp_em, 2.5*h*h*h*h/(64.0*eta_H));
                        // NEW (2026-11-06): the eta_H term acts on B', so B relaxes at eta_H lam^2/(1 + d_e^2 lam),
                        // lam_max = 8/h^2; the factor 2.0 keeps |R_RK3| <= 1 together with |omega dt| <= 1
                        if (de > 0.0)
                            dtp_em = std::min(dtp_em, 2.0*(1.0 + 8.0*de*de/(h*h))*h*h*h*h/(64.0*eta_H));
                        else
                            dtp_em = std::min(dtp_em, 2.5*h*h*h*h/(64.0*eta_H));
                    }
                }

            if (std::isnan(dtp_em) || std::isinf(dtp_em))
                MPI_Abort(MPI_COMM_WORLD, 1);

            MPI_Allreduce(&dtp_em, &dt_em, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            dt_em = std::min(dt_em, t_target - t_em);   // land exactly on target

            for (j = 2; j < J+2; j++) lambda_em(j) = dt_em/dx_cell(j);
            for (k = 2; k < K+2; k++) mu_em(k)     = dt_em/dy_cell(k);

            // --- SSP-RK3 (Shu-Osher): stable for the imaginary whistler eigenvalues
            //     without upwind damping, |omega*dt_em| <= sqrt(3) ---
            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                {
                    Bx0(j,k) = B1_xf(j,k,0);
                    By0(j,k) = B2_yf(j,k,0);
                    Bz0(j,k) = un(j,k,6);
                }
            // NEW (2026-10-31): with d_e > 0 the RK3 base state is B'
            if (de > 0.0)
                for (k = 2; k < K+2; k++)
                    for (j = 2; j < J+2; j++)
                    {
                        Bx0(j,k) = (*Bp1)(j,k,0);
                        By0(j,k) = (*Bp2)(j,k,0);
                        Bz0(j,k) = (*Bpz)(j,k,0);
                    }

            // (Butcher weights of this scheme, for the dissipation integrals: 1/6, 1/6, 2/3)
            // ORIGINAL (pre-2026-11-06):
            // em_rhs(1.0/6.0);  em_stage(0.0,     1.0);       // u1 = u0 + L(u0)
            // em_rhs(1.0/6.0);  em_stage(0.75,    0.25);      // u2 = 3/4 u0 + 1/4 (u1 + L(u1))
            // em_rhs(2.0/3.0);  em_stage(1.0/3.0, 2.0/3.0);   // u3 = 1/3 u0 + 2/3 (u2 + L(u2))
            // NEW (2026-11-06): d_e > 0 runs the same three stages through efb_substep (energy exchange)
            if (efb) efb_substep();
            else
            {
            em_rhs(1.0/6.0);  em_stage(0.0,     1.0);       // u1 = u0 + L(u0)
            em_rhs(1.0/6.0);  em_stage(0.75,    0.25);      // u2 = 3/4 u0 + 1/4 (u1 + L(u1))
            em_rhs(2.0/3.0);  em_stage(1.0/3.0, 2.0/3.0);   // u3 = 1/3 u0 + 2/3 (u2 + L(u2))
            }

            nsub++;
            t_em += dt_em;
        }
        
        if (id == 0) printf("  EM sub-steps: %d\n", nsub);

    };

    // ========================================================================
    // Strang / Verlet splitting (symmetric, neutrally stable for Alfven waves):
    //   EM   0    -> dt/2   with rho, v frozen at t^n
    //   fluid RK2 full dt   with B frozen at t^{n+1/2}
    //   EM   dt/2 -> dt     with rho, v frozen at t^{n+1}
    // (the old interleaving had |amplification| > 1 for every k: growth ~ k^2 dt / 8)
    // ========================================================================

    // ORIGINAL (pre-2026-11-06): the step sequence below, now with the Formulation B calls for d_e > 0
    // // ---- EM half step 1: v at t^n ----
    // reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    // freeze_fluid();
    // em_subcycle(0.5*dt_fluid);
    //
    // // ---- fluid RK2, B frozen at t^{n+1/2} ----
    // compute_J();
    // reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    // C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Jc, C0);
    // for (int n = 0; n < n_fluid; n++)
    // {
    //     l = fluid_comp[n];
    //     for (k = 2; k < K+2; k++)
    //         for (j = 2; j < J+2; j++)
    //             un(j,k,l) += C0(j,k,l);
    // }
    // boundary_conditions(un, parameters, id, p);
    //
    // reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    // C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Jc, C1);   // same Jc: B unchanged
    // for (int n = 0; n < n_fluid; n++)
    // {
    //     l = fluid_comp[n];
    //     for (k = 2; k < K+2; k++)
    //         for (j = 2; j < J+2; j++)
    //             un(j,k,l) += 0.5*(C1(j,k,l) - C0(j,k,l));
    // }
    // boundary_conditions(un, parameters, id, p);
    //
    // // ---- EM half step 2: v at t^{n+1} ----
    // reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    // freeze_fluid();
    // em_subcycle(dt_fluid);

    // NEW (2026-11-06): Formulation B helpers (d_e > 0). e_f = 1/2 rho v^2 + p/(gamma-1); slot 7 stays the total
    // energy E = e_f + 1/2|B|^2 (cell B) and is rebuilt from e_f after every fluid update and EM half step.
    const double gam_p = parameters(0);
    doublearray2d Pc(J+4,K+4), Pl(J+4,K+4), Ekeep(J+4,K+4), Et0(J+4,K+4);
    auto ef_c = [&](long a, long b) { return un(a,b,7) - 0.5*(un(a,b,4)*un(a,b,4) + un(a,b,5)*un(a,b,5) + un(a,b,6)*un(a,b,6)); };
    auto setEf = [&](long a, long b, double e) { un(a,b,7) = e + 0.5*(un(a,b,4)*un(a,b,4) + un(a,b,5)*un(a,b,5) + un(a,b,6)*un(a,b,6)); };
    // e_f rate: d_t e_f + div[(e_f + p) v] = 0 (the field-fluid exchange is applied in the EM step); flux on the
    // SD2 interface states, central part Rc and LxF part Rl with the fluid spectral radius (as Hx/Hy_flux_2d_SD2)
    auto efl_rhs = [&](doublearray2d& Rc, doublearray2d& Rl) {
        doublearray1d uv(L);
        auto st = [&](const doublearray3d& U, long a, long b, double& e, double& fx, double& fy, double& rx, double& ry) {
            for (long q = 0; q < L; q++) uv(q) = U(a,b,q);
            double r = uv(0), vx = uv(1)/r, vy = uv(2)/r, vz = uv(3)/r;
            double kin = 0.5*r*(vx*vx + vy*vy + vz*vz);
            double pp = (gam_p - 1.0)*(uv(7) - kin - 0.5*(uv(4)*uv(4) + uv(5)*uv(5) + uv(6)*uv(6)));
            e = kin + pp/(gam_p - 1.0); fx = (e + pp)*vx; fy = (e + pp)*vy;
            spectral_radii(uv, parameters, rx, ry); };
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            {
                double el, fxl, fyl, rxl, ryl, er, fxr, fyr, rxr, ryr, Fc[4], Fl[4];
                st(u_E, j, k, el, fxl, fyl, rxl, ryl); st(u_W, j+1, k, er, fxr, fyr, rxr, ryr);
                Fc[0] = 0.5*(fxl + fxr); Fl[0] = -0.5*std::max(rxl, rxr)*(er - el);
                st(u_E, j-1, k, el, fxl, fyl, rxl, ryl); st(u_W, j, k, er, fxr, fyr, rxr, ryr);
                Fc[1] = 0.5*(fxl + fxr); Fl[1] = -0.5*std::max(rxl, rxr)*(er - el);
                st(u_N, j, k, el, fxl, fyl, rxl, ryl); st(u_S, j, k+1, er, fxr, fyr, rxr, ryr);
                Fc[2] = 0.5*(fyl + fyr); Fl[2] = -0.5*std::max(ryl, ryr)*(er - el);
                st(u_N, j, k-1, el, fxl, fyl, rxl, ryl); st(u_S, j, k, er, fxr, fyr, rxr, ryr);
                Fc[3] = 0.5*(fyl + fyr); Fl[3] = -0.5*std::max(ryl, ryr)*(er - el);
                Rc(j,k) = -(Fc[0] - Fc[1])/dx_cell(j) - (Fc[2] - Fc[3])/dy_cell(k);
                Rl(j,k) = -(Fl[0] - Fl[1])/dx_cell(j) - (Fl[2] - Fl[3])/dy_cell(k);
            }
    };
    // e_f is held through an EM half step (B changes, e_f does not), then the exchange of that half is added
    auto keep_ef = [&]() { for (k = 2; k < K+2; k++) for (j = 2; j < J+2; j++) Ekeep(j,k) = ef_c(j,k); };
    auto restore_ef = [&]() { for (k = 2; k < K+2; k++) for (j = 2; j < J+2; j++) setEf(j, k, Ekeep(j,k));
                              boundary_conditions(un, parameters, id, p); };
    auto efb_apply = [&]() {
        for (k = 2; k < K+2; k++) for (j = 2; j < J+2; j++) { setEf(j, k, ef_c(j,k) + (*eQ)(j,k)); (*eQ)(j,k) = 0.0; }
        boundary_conditions(un, parameters, id, p); };
    // electron-inertia stress -div T_e, T_e = d_e^2 (J J - 1/2 |J|^2 I), compact J (C^T B) averaged to the fluid
    // faces; B is frozen in the fluid step, so one evaluation serves both Heun stages (x dt_fluid)
    static doublearray3d *eDm = 0;
    if (efb && eDm == 0) eDm = new doublearray3d(J+4, K+4, 3);
    auto efb_stress = [&]() {
        efb_snap(*eT, false);
        const doublearray3d& B = *eT;
        const double ddx = dx_cell(2), ddy = dy_cell(2), d2 = de*de;
        auto Jxy = [&](long a, long b) { return  (B(a,b+1,2) - B(a,b,2))/ddy; };                     // Jx on y-face (a,b+1/2)
        auto Jyx = [&](long a, long b) { return -(B(a+1,b,2) - B(a,b,2))/ddx; };                     // Jy on x-face (a+1/2,b)
        auto Jzn = [&](long a, long b) { return (B(a+1,b,1) - B(a,b,1))/ddx - (B(a,b+1,0) - B(a,b,0))/ddy; };   // node
        auto fx = [&](long a, long b, double* f) {   // x-face (a+1/2, b)
            double jx = 0.25*(Jxy(a,b) + Jxy(a,b-1) + Jxy(a+1,b) + Jxy(a+1,b-1));
            double jy = Jyx(a,b), jz = 0.5*(Jzn(a,b) + Jzn(a,b-1)), j2 = jx*jx + jy*jy + jz*jz;
            f[0] = d2*(jx*jx - 0.5*j2); f[1] = d2*jx*jy; f[2] = d2*jx*jz; };
        auto fy = [&](long a, long b, double* f) {   // y-face (a, b+1/2)
            double jy = 0.25*(Jyx(a,b) + Jyx(a-1,b) + Jyx(a,b+1) + Jyx(a-1,b+1));
            double jx = Jxy(a,b), jz = 0.5*(Jzn(a,b) + Jzn(a-1,b)), j2 = jx*jx + jy*jy + jz*jz;
            f[0] = d2*jy*jx; f[1] = d2*(jy*jy - 0.5*j2); f[2] = d2*jy*jz; };
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            {
                double fe[3], fw_[3], fn[3], fs[3];
                fx(j, k, fe); fx(j-1, k, fw_); fy(j, k, fn); fy(j, k-1, fs);
                for (l = 0; l < 3; l++)
                    (*eDm)(j,k,l) = -dt_fluid*((fe[l] - fw_[l])/ddx + (fn[l] - fs[l])/ddy);
            }
    };

    // NEW (2026-11-06): step sequence; every efb branch is d_e > 0 only, so d_e = 0 is the original sequence
    // ---- EM half step 1: v at t^n ----
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    freeze_fluid();
    if (efb) keep_ef();
    em_subcycle(0.5*dt_fluid);
    if (efb) { restore_ef(); efb_apply(); }   // B changed, e_f did not; then e_f += exchange of this EM half

    // ---- fluid RK2, B frozen at t^{n+1/2} ----
    compute_J();
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Jc, C0);
    if (efb)   // J x B' = J x B + d_e^2 J x curl J as -div T_e; e_f0 and its rate from the same state
    {
        efb_stress();
        for (k = 2; k < K+2; k++) for (j = 2; j < J+2; j++) for (l = 0; l < 3; l++) C0(j,k,l+1) += (*eDm)(j,k,l);
        efl_rhs(Pc, Pl);
        for (k = 2; k < K+2; k++) for (j = 2; j < J+2; j++)
        { Ekeep(j,k) = ef_c(j,k); Et0(j,k) = Pc(j,k) + Pl(j,k); }
    }
    for (int n = 0; n < n_fluid; n++)
    {
        l = fluid_comp[n];
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
                un(j,k,l) += C0(j,k,l);
    }
    if (efb)   // e_f1 = e_f0 + dt R0
        for (k = 2; k < K+2; k++) for (j = 2; j < J+2; j++) setEf(j, k, Ekeep(j,k) + dt_fluid*Et0(j,k));
    boundary_conditions(un, parameters, id, p);

    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Jc, C1);   // same Jc: B unchanged
    if (efb)   // same stress (B frozen); e_f rate at the stage-1 state
    {
        for (k = 2; k < K+2; k++) for (j = 2; j < J+2; j++) for (l = 0; l < 3; l++) C1(j,k,l+1) += (*eDm)(j,k,l);
        efl_rhs(Pc, Pl);
    }
    for (int n = 0; n < n_fluid; n++)
    {
        l = fluid_comp[n];
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
                un(j,k,l) += 0.5*(C1(j,k,l) - C0(j,k,l));
    }
    if (efb)   // Heun, e_f = e_f0 + dt/2 (R0 + R1)
        for (k = 2; k < K+2; k++) for (j = 2; j < J+2; j++)
            setEf(j, k, Ekeep(j,k) + 0.5*dt_fluid*(Et0(j,k) + Pc(j,k) + Pl(j,k)));
    boundary_conditions(un, parameters, id, p);

    // ---- EM half step 2: v at t^{n+1} ----
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    freeze_fluid();
    if (efb) keep_ef();
    em_subcycle(dt_fluid);
    if (efb) { restore_ef(); efb_apply(); }
    if (efb && efb_log)   // EFB_LOG: per fluid step sums (x dA) of the exchange, its J.E part and R_RK
    {
        double a[3] = {efb_sQ, efb_sJE, efb_sR}, g[3];
        MPI_Allreduce(a, g, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        static double t_e = 0.0; t_e += dt_fluid;
        if (id == 0) printf("EFB t %.5f Q %+.9e JE %+.9e RRK %+.9e\n", t_e, g[0], g[1], g[2]);
        efb_sQ = efb_sJE = efb_sR = 0.0;
    }

    // NEW (2026-11-06): EDRIFT=1 (environment): total energy sum e_f + 1/2 <B, B'> per fluid step (B' = B at d_e = 0);
    // diagnostic only, it changes nothing in the solution
    if (getenv("EDRIFT"))
    {
        double ef = 0.0, eb = 0.0, g2[2], l2[2];
        const double h1 = dx_cell(2), h2 = dy_cell(2);
        for (k = 2; k < K+2; k++) for (j = 2; j < J+2; j++)
        {
            ef += ef_c(j,k)*h1*h2;
            if (de > 0.0) eb += 0.5*(B1_xf(j,k,0)*(*Bp1)(j,k,0) + B2_yf(j,k,0)*(*Bp2)(j,k,0) + un(j,k,6)*(*Bpz)(j,k,0))*h1*h2;
            else          eb += 0.5*(B1_xf(j,k,0)*B1_xf(j,k,0) + B2_yf(j,k,0)*B2_yf(j,k,0) + un(j,k,6)*un(j,k,6))*h1*h2;
        }
        l2[0] = ef; l2[1] = eb;
        MPI_Allreduce(l2, g2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        static double t_ed = 0.0; t_ed += dt_fluid;
        if (id == 0) printf("EDRIFT t %.5f Ef %.12e EB %.12e Etot %.12e\n", t_ed, g2[0], g2[1], g2[0] + g2[1]);
    }
}
