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
    const double c_whistler = (n_par > 7) ? parameters(7) : 0.0;   // whistler fraction in the EM upwinding speed
    const bool   walls  = ((int) parameters(1) == 1);   // reconnection: conducting walls in y
    const int    yw_J   = walls ? 3 : 0;                  // ghost-fill modes for boundary_conditions_ct
    const int    yw_Bx  = walls ? 1 : 0;
    const int    yw_By  = walls ? 2 : 0;
    // NEW (2026-10-31): electron skin depth (electron inertia), parameters(8), default 0 = MHD2d
    const double de     = (n_par > 8) ? parameters(8) : 0.0;
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
                  eta_H > 0.0 ? &E_hyper : 0, c_whistler);
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
                      eta_H > 0.0 ? &E_hyper : 0, c_whistler);
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
                    if (rate > 0.0)
                        dtp_em = std::min(dtp_em, cfl_em/rate);
                    // parabolic limit of the explicit 4th-order term (SSP-RK3 real-axis
                    // stability ~2.5; max eigenvalue of lap^2 on the 5-point lattice 64/h^4)
                    if (eta_H > 0.0)
                    {
                        double h = std::min(dx_cell(j), dy_cell(k));
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
            em_rhs(1.0/6.0);  em_stage(0.0,     1.0);       // u1 = u0 + L(u0)
            em_rhs(1.0/6.0);  em_stage(0.75,    0.25);      // u2 = 3/4 u0 + 1/4 (u1 + L(u1))
            em_rhs(2.0/3.0);  em_stage(1.0/3.0, 2.0/3.0);   // u3 = 1/3 u0 + 2/3 (u2 + L(u2))

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

    // ---- EM half step 1: v at t^n ----
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    freeze_fluid();
    em_subcycle(0.5*dt_fluid);

    // ---- fluid RK2, B frozen at t^{n+1/2} ----
    compute_J();
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Jc, C0);
    for (int n = 0; n < n_fluid; n++)
    {
        l = fluid_comp[n];
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
                un(j,k,l) += C0(j,k,l);
    }
    boundary_conditions(un, parameters, id, p);

    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Jc, C1);   // same Jc: B unchanged
    for (int n = 0; n < n_fluid; n++)
    {
        l = fluid_comp[n];
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
                un(j,k,l) += 0.5*(C1(j,k,l) - C0(j,k,l));
    }
    boundary_conditions(un, parameters, id, p);

    // ---- EM half step 2: v at t^{n+1} ----
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    freeze_fluid();
    em_subcycle(dt_fluid);
}