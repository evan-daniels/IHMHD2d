////////////////////////////////////////////////////////////////////////////////
//
// C_flux_em.cc — EM flux computation for the CT sub-stepper
//
// Computes:
//   (a) Upwinded Ez at x-faces and y-faces (for Bx/By CT update)
//   (b) Corner Ez from 4-point average (drives Faraday's law)
//   (c) Cx_em, Cy_em: Faraday increments for Bx and By face arrays
//   (d) Cz_em: FV increment for cell-centered Bz (using upwinded Ey, Ex)
//
// Upwinding speed: Alfven + whistler,  |B|/sqrt(rho) + di*|B|/(rho*dx)
// (the Alfven term supplies dissipation when di = 0; without it the
//  forward-Euler EM update is centered and unconditionally unstable)
//
// ywall = 1 (reconnection): perfectly conducting walls at the y-boundaries.
//   Tangential E vanishes on a conductor, so
//     - Ez = 0 at the wall corners (k = 1 bottom, k = K+1 top): By on the wall
//       faces then never changes, keeping normal B = 0 and div B ~ 0 next to the wall
//     - Ex = 0 on the wall y-faces: zero y-flux of Bz through the wall
// ywall = 0: periodic in y (Orszag-Tang), unchanged behavior.
//
// s_B  scales the upwind dissipation in the CT Ez (Bx, By update),
// s_Bz scales the upwind dissipation in the Ey, -Ex fluxes (Bz update).
// Both are 1 by default (parameters(4), parameters(5) in the input file).
// c_whistler (parameters(7), default 0) adds that fraction of the whistler
// speed di*|B|/(rho*dx) to the Alfven upwinding speed at every face.
//
// Hyperresistivity, E += -eta_H lap(J) (parameters(6), default 0): hyper_em() below
// builds J on the staggered lattice where each E component lives and returns
// -eta_H lap(J) there; C_flux_em adds it before the conducting-wall zeroing.
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"
#include <cmath>

using namespace std;

void CENTPACK::C_flux_em(
    const doublearray2d& rho_N, const doublearray2d& rho_S,
    const doublearray2d& rho_E, const doublearray2d& rho_W,
    const doublearray3d& B_N, const doublearray3d& B_S,
    const doublearray3d& B_E, const doublearray3d& B_W,
    const doublearray3d& E_N, const doublearray3d& E_S,
    const doublearray3d& E_E, const doublearray3d& E_W,
    const doublearray1d& lambda_em, const doublearray1d& mu_em,
    const doublearray1d& dx_interface, const doublearray1d& dy_interface,
    const double& di,
    doublearray2d& Cx_em, doublearray2d& Cy_em, doublearray2d& Cz_em,
    const int& ywall, const double& s_B, const double& s_Bz,
    // ORIGINAL (pre-2026-11-06):
    // const doublearray3d* E_hyper, const double& c_whistler)
    // NEW (2026-11-06): c_whistler is ignored (removed option, kept in the signature); c_wi = c_whistler_inertia
    // (used only when de > 0); Eedge, if given, receives the edge E* (Ez nodes, Ey x-faces, Ex y-faces)
    const doublearray3d* E_hyper, const double& c_whistler,
    const double& c_wi, const double& de, doublearray3d* Eedge)

{

    // c_whistler (parameters(7), default 0): fraction of the whistler speed in the LxF upwinding
    //   1.0 = old behaviour (needed with forward-Euler EM)
    //   0.0 = Alfven-only upwinding (relies on RK3 EM sub-cycle for whistler stability)

    long j, k, J, K;
    J = B_N.getIndex1Size() - 4;
    K = B_N.getIndex2Size() - 4;

    double Bl, Br, rl, rr, al, ar, ax;
    double Bb, Bt, rb, rt, ab, at_v, ay;

    // Internal face/corner arrays
    doublearray2d Ez_xf(J+4, K+4);    // upwinded Ez at x-faces
    doublearray2d Ez_yf(J+4, K+4);    // upwinded Ez at y-faces
    doublearray2d Fx_xf(J+4, K+4);    // Ey x-flux for Bz update
    doublearray2d Gy_yf(J+4, K+4);    // -Ex y-flux for Bz update
    doublearray2d Ez_corner(J+4, K+4); // corner Ez

    // ===== (a) Upwinded Ez and Ey at x-faces =====
    // x-face at (j+1/2, k): Left = E_E(j,k), Right = E_W(j+1,k)
    for (k = 1; k < K+3; k++)
    {
        for (j = 1; j < J+2; j++)
        {
            Bl = sqrt(B_E(j,k,0)*B_E(j,k,0) + B_E(j,k,1)*B_E(j,k,1) + B_E(j,k,2)*B_E(j,k,2));
            Br = sqrt(B_W(j+1,k,0)*B_W(j+1,k,0) + B_W(j+1,k,1)*B_W(j+1,k,1) + B_W(j+1,k,2)*B_W(j+1,k,2));
            rl = std::max(rho_E(j,k), 0.1);
            rr = std::max(rho_W(j+1,k), 0.1);
            // al = di * Bl / (rl * dx_interface(j)); ORIGINAL
            // ar = di * Br / (rr * dx_interface(j)); ORIGINAL
            // al = Bl / sqrt(rl) + di * Bl / (rl * dx_interface(j));
            // ar = Br / sqrt(rr) + di * Br / (rr * dx_interface(j));
            // ORIGINAL (pre-2026-11-06):
            // al = Bl / sqrt(rl) + c_whistler * di * Bl / (rl * dx_interface(j));
            // ar = Br / sqrt(rr) + c_whistler * di * Br / (rr * dx_interface(j));
            // NEW (2026-11-06): uncapped whistler term removed (identical to the original with c_whistler = 0)
            al = Bl / sqrt(rl);
            ar = Br / sqrt(rr);
            // NEW (2026-11-06): capped whistler upwind, d_e > 0 only: + c_wi d_i |B| khat / (rho (1 + d_e^2 khat^2)), khat = pi/h
            if (de > 0.0 && c_wi > 0.0)
            {
                const double kh = M_PI/dx_interface(j), cap = 1.0/(1.0 + de*de*kh*kh);
                al += c_wi*di*Bl*kh*cap/rl;
                ar += c_wi*di*Br*kh*cap/rr;
            }


            ax = std::max(al, ar);

            // Ez: dBy/dt = +dEz/dx, so flux of By is -Ez
            // LxF gives: Ez = avg + 0.5*a*(By_R - By_L)
            Ez_xf(j,k) = 0.5*(E_E(j,k,2) + E_W(j+1,k,2))
                        + 0.5*s_B*ax*(B_W(j+1,k,1) - B_E(j,k,1));

            // Ey: dBz/dt = -dEy/dx + ..., so x-flux of Bz is F = Ey
            // LxF gives: F = avg - 0.5*a*(Bz_R - Bz_L)
            Fx_xf(j,k) = 0.5*(E_E(j,k,1) + E_W(j+1,k,1))
                        - 0.5*s_Bz*ax*(B_W(j+1,k,2) - B_E(j,k,2));
        }
    }

    // ===== (b) Upwinded Ez and -Ex at y-faces =====
    // y-face at (j, k+1/2): Bottom = E_N(j,k), Top = E_S(j,k+1)
    for (k = 1; k < K+2; k++)
    {
        for (j = 1; j < J+3; j++)
        {
            Bb = sqrt(B_N(j,k,0)*B_N(j,k,0) + B_N(j,k,1)*B_N(j,k,1) + B_N(j,k,2)*B_N(j,k,2));
            Bt = sqrt(B_S(j,k+1,0)*B_S(j,k+1,0) + B_S(j,k+1,1)*B_S(j,k+1,1) + B_S(j,k+1,2)*B_S(j,k+1,2));
            rb = std::max(rho_N(j,k), 0.1);
            rt = std::max(rho_S(j,k+1), 0.1);
            // ab = di * Bb / (rb * dy_interface(k)); ORIGINAL
            // at_v = di * Bt / (rt * dy_interface(k)); ORIGINAL
            // ab   = Bb / sqrt(rb) + di * Bb / (rb * dy_interface(k));
            // at_v = Bt / sqrt(rt) + di * Bt / (rt * dy_interface(k));
            // ORIGINAL (pre-2026-11-06):
            // ab   = Bb / sqrt(rb) + c_whistler * di * Bb / (rb * dy_interface(k));
            // at_v = Bt / sqrt(rt) + c_whistler * di * Bt / (rt * dy_interface(k));
            // NEW (2026-11-06): uncapped whistler term removed (identical to the original with c_whistler = 0)
            ab   = Bb / sqrt(rb);
            at_v = Bt / sqrt(rt);
            // NEW (2026-11-06): capped whistler upwind, d_e > 0 only
            if (de > 0.0 && c_wi > 0.0)
            {
                const double kh = M_PI/dy_interface(k), cap = 1.0/(1.0 + de*de*kh*kh);
                ab   += c_wi*di*Bb*kh*cap/rb;
                at_v += c_wi*di*Bt*kh*cap/rt;
            }


            ay = std::max(ab, at_v);

            // Ez: dBx/dt = -dEz/dy, so flux of Bx is +Ez
            // LxF gives: Ez = avg - 0.5*a*(Bx_T - Bx_B)
            Ez_yf(j,k) = 0.5*(E_N(j,k,2) + E_S(j,k+1,2))
                        - 0.5*s_B*ay*(B_S(j,k+1,0) - B_N(j,k,0));

            // -Ex: dBz/dt = ... + dEx/dy, so y-flux of Bz is G = -Ex
            // LxF gives: G = avg(-Ex) - 0.5*a*(Bz_T - Bz_B)
            Gy_yf(j,k) = -0.5*(E_N(j,k,0) + E_S(j,k+1,0))
                         - 0.5*s_Bz*ay*(B_S(j,k+1,2) - B_N(j,k,2));
        }
    }

    // Hyperresistive Ey (x-faces) and Ex (y-faces): F = Ey, G = -Ex
    if (E_hyper)
    {
        for (k = 1; k < K+3; k++)
            for (j = 1; j < J+2; j++)
                Fx_xf(j,k) += (*E_hyper)(j,k,1);
        for (k = 1; k < K+2; k++)
            for (j = 1; j < J+3; j++)
                Gy_yf(j,k) -= (*E_hyper)(j,k,0);
    }

    // Conducting walls: Ex = 0 on the wall y-faces (bottom k = 1, top k = K+1)
    // -> no flux of Bz through the wall
    if (ywall == 1)
    {
        for (j = 1; j < J+3; j++)
        {
            Gy_yf(j,1)   = 0.0;
            Gy_yf(j,K+1) = 0.0;
        }
    }

    // ===== (c) Corner Ez: 4-point average of surrounding face values =====
    for (k = 1; k < K+2; k++)
        for (j = 1; j < J+2; j++)
            Ez_corner(j,k) = 0.25*(Ez_xf(j,k) + Ez_xf(j,k+1)
                                  + Ez_yf(j,k) + Ez_yf(j+1,k));

    // Hyperresistive Ez at the corners
    if (E_hyper)
        for (k = 1; k < K+2; k++)
            for (j = 1; j < J+2; j++)
                Ez_corner(j,k) += (*E_hyper)(j,k,2);

    // Conducting walls: Ez = 0 at the wall corners (bottom k = 1, top k = K+1)
    // -> By on the wall faces stays 0, Bx next to the wall sees the conductor
    if (ywall == 1)
    {
        for (j = 1; j < J+2; j++)
        {
            Ez_corner(j,1)   = 0.0;
            Ez_corner(j,K+1) = 0.0;
        }
    }

    // NEW (2026-11-06): edge-field export for Formulation B (energy flux S = E* x B from the exact edge E*)
    if (Eedge)
        for (k = 0; k < K+4; k++)
            for (j = 0; j < J+4; j++)
            {
                bool okz = (j >= 1 && j <= J+1 && k >= 1 && k <= K+1);
                bool oky = (j >= 1 && j <= J+1 && k >= 1 && k <= K+2);
                bool okx = (j >= 1 && j <= J+2 && k >= 1 && k <= K+1);
                (*Eedge)(j,k,0) = okz ? Ez_corner(j,k) : 0.0;
                (*Eedge)(j,k,1) = oky ? Fx_xf(j,k) : 0.0;
                (*Eedge)(j,k,2) = okx ? -Gy_yf(j,k) : 0.0;
            }

    // ===== (d) Cx_em: Bx face increment, dBx/dt = -dEz/dy =====
    for (k = 2; k < K+2; k++)
        for (j = 2; j < J+2; j++)
            Cx_em(j,k) = -mu_em(k) * (Ez_corner(j,k) - Ez_corner(j,k-1));

    // ===== (e) Cy_em: By face increment, dBy/dt = +dEz/dx =====
    for (k = 2; k < K+2; k++)
        for (j = 2; j < J+2; j++)
            Cy_em(j,k) = lambda_em(j) * (Ez_corner(j,k) - Ez_corner(j-1,k));

    // ===== (f) Cz_em: Bz cell-center increment =====
    // dBz/dt + d(Ey)/dx + d(-Ex)/dy = 0
    // Bz_new = Bz_old - lambda*(Fx_right - Fx_left) - mu*(Gy_top - Gy_bottom)
    for (k = 2; k < K+2; k++)
        for (j = 2; j < J+2; j++)
            Cz_em(j,k) = -lambda_em(j) * (Fx_xf(j,k) - Fx_xf(j-1,k))
                         - mu_em(k) * (Gy_yf(j,k) - Gy_yf(j,k-1));
}

////////////////////////////////////////////////////////////////////////////////
//
// hyper_em -- hyperresistive electric field  E_H = -eta_H lap(J)  on the CT lattice
//
//   Jz at corners (j+1/2, k+1/2) = d(By_face)/dx - d(Bx_face)/dy   -> E_hyper(j,k,2)
//   Jy at x-faces (j+1/2, k)     = -d(Bz)/dx                         -> E_hyper(j,k,1)
//   Jx at y-faces (j, k+1/2)     =  d(Bz)/dy                         -> E_hyper(j,k,0)
//
// Each J is differenced from the ghost-filled B, so its wall values follow the
// existing B wall ghosts: Jz = 0 at the wall corners (By = 0 on the wall faces,
// Bx even), Jx = 0 on the wall faces (Bz even), Jy even across the wall -- the
// same parity as the cell-centred Jc (boundary_conditions_ct mode 3). x is periodic
// through the B ghosts. lap is the compact 5-point stencil on the same lattice
// (uniform mesh assumed, as in mesh.cc).
//
// Returns the rank-local dissipation  eta_H * sum J.(-lap J) dA  over owned points,
// which by summation by parts equals eta_H * sum |grad_h J|^2 dA >= 0.
//
////////////////////////////////////////////////////////////////////////////////

double CENTPACK::hyper_em(const doublearray3d& B1_xf, const doublearray3d& B2_yf,
    const doublearray3d& un, const doublearray1d& dx_interface,
    const doublearray1d& dy_interface, const double& eta_H, doublearray3d& E_hyper)
{
    long j, k, J, K;
    J = un.getIndex1Size() - 4;
    K = un.getIndex2Size() - 4;

    const double dx = dx_interface(2), dy = dy_interface(2);
    const double idx2 = 1.0/(dx*dx), idy2 = 1.0/(dy*dy);

    doublearray2d Jz_c(J+4,K+4), Jy_xf(J+4,K+4), Jx_yf(J+4,K+4);

    for (k = 0; k < K+3; k++)
        for (j = 0; j < J+3; j++)
            Jz_c(j,k) = (B2_yf(j+1,k,0) - B2_yf(j,k,0))/dx - (B1_xf(j,k+1,0) - B1_xf(j,k,0))/dy;
    for (k = 0; k < K+4; k++)
        for (j = 0; j < J+3; j++)
            Jy_xf(j,k) = -(un(j+1,k,6) - un(j,k,6))/dx;
    for (k = 0; k < K+3; k++)
        for (j = 0; j < J+4; j++)
            Jx_yf(j,k) = (un(j,k+1,6) - un(j,k,6))/dy;

    auto lap = [&](const doublearray2d& f, long jj, long kk)
    {
        return (f(jj+1,kk) - 2.0*f(jj,kk) + f(jj-1,kk))*idx2
             + (f(jj,kk+1) - 2.0*f(jj,kk) + f(jj,kk-1))*idy2;
    };

    for (k = 1; k < K+2; k++)
        for (j = 1; j < J+2; j++)
            E_hyper(j,k,2) = -eta_H*lap(Jz_c, j, k);
    for (k = 1; k < K+3; k++)
        for (j = 1; j < J+2; j++)
            E_hyper(j,k,1) = -eta_H*lap(Jy_xf, j, k);
    for (k = 1; k < K+2; k++)
        for (j = 1; j < J+3; j++)
            E_hyper(j,k,0) = -eta_H*lap(Jx_yf, j, k);

    // owned points: corners / faces (j, k), j, k = 2..J+1, 2..K+1
    // (the bottom wall corner/face k = 1 carries J = 0 and is left out)
    double d = 0.0;
    for (k = 2; k < K+2; k++)
        for (j = 2; j < J+2; j++)
            d += Jz_c(j,k)*E_hyper(j,k,2) + Jy_xf(j,k)*E_hyper(j,k,1) + Jx_yf(j,k)*E_hyper(j,k,0);
    return d*dx*dy;
}
