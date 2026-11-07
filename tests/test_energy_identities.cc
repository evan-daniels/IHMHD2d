////////////////////////////////////////////////////////////////////////////////
//
// test_energy_identities.cc -- unit tests (a), (b), (4) of tools/design/inertia_energy.md
// (scratch only). Uses the real boundary_conditions_ct, helmholtz_apply_walls and
// helmholtz_solve_walls; CT curl stencils copied from C_flux_em (d)-(f).
//
// (a) symmetry of H = 1 - de^2 lap_h with wall images, <a,Hb> - <Ha,b> and the same
//     for H^-1 (the FFT solve), with fields supported (i) everywhere, (ii) only in the
//     3 rows next to each wall, (iii) only in the interior half; plus the per-cell
//     Green residual a L b - b L a - div G and the wall flux G_wall from the ghosts.
// (b) staggered product rule  B.(C E) = div_h S + avg(J.E)  per cell (cell assignment (3)),
//     and the adjoint <B, C E> = <C^T B, E>.
// (4) the full local field-energy law with the inertia flux, bdot = H^-1(-C E):
//     wdot_direct + div S + avg(J.E) + 1/2 de^2 div G(bdot, b) = 0 per cell.
//
// Run: mpirun -np P ./test_energy_identities   (P divides NX)
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"
static double g_worst = 0.0;   // NEW (2026-11-06): largest relative residual (rank 0), for the exit status
#include <random>
#include <vector>
#include <cstdint>

using namespace std;
using namespace CENTPACK;

static double gmax(double v) { double g; MPI_Allreduce(&v, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return g; }
static double gsum(double v) { double g; MPI_Allreduce(&v, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return g; }

// rank-independent random number for (field tag, global j, k)
static double rnd(uint64_t tag, long jg, long k)
{
    uint64_t s = tag*1000003ULL + (uint64_t) jg*7919ULL + (uint64_t) k*104729ULL + 0x9E3779B97F4A7C15ULL;
    mt19937_64 g(s); g.discard(3);
    return uniform_real_distribution<double>(-1.0, 1.0)(g);
}

static long NX = 128, K = 64, J, id_, p_;
static double dx = 0.05, dy = 0.05;

struct F { doublearray3d a; F() : a(J+4, K+4, 1) { for (long j = 0; j < J+4; j++) for (long k = 0; k < K+4; k++) a(j,k,0) = 0.0; }
           double& operator()(long j, long k) { return a(j,k,0); } };

// fill unknowns: kind 1 (cell rows 2..K+1), kind 2 (face rows 1..K+1, walls 0 if walls)
// support: 0 all, 1 wall band (3 rows next to each wall), 2 interior half
static void fill(F& f, uint64_t tag, int yw, int support)
{
    long k0 = (yw == 2) ? 1 : 2, k1 = K+1;
    for (long j = 2; j < J+2; j++)
        for (long k = k0; k <= k1; k++)
        {
            long jg = id_*J + j - 2;
            bool in = true;
            long dw = std::min(k - k0, k1 - k);           // distance from the outermost row
            if (support == 1) in = (dw < 3 + (yw == 2 ? 1 : 0));
            if (support == 2) in = (k > K/4 && k < 3*K/4);
            f(j,k) = in ? rnd(tag, jg, k) : 0.0;
            if (yw == 2 && (k == 1 || k == K+1)) f(j,k) = 0.0;
        }
    boundary_conditions_ct(f.a, (int) id_, (int) p_, yw);
}

// sum over the unknowns of a lattice (Y lattice with walls: rows 2..K; else 2..K+1 / 1..K for periodic Y)
static double dot(F& a, F& b, int yw, int rows = 0 /*0 all, 1 wall rows, 2 interior*/)
{
    long k0 = 2, k1 = (yw == 2) ? K : K+1;
    double s = 0.0;
    for (long j = 2; j < J+2; j++)
        for (long k = k0; k <= k1; k++)
        {
            bool w = (k <= k0+2 || k >= k1-2);
            if (rows == 1 && !w) continue;
            if (rows == 2 && w) continue;
            s += a(j,k)*b(j,k);
        }
    return gsum(s*dx*dy);
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int id, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    id_ = id; p_ = p; J = NX/p;
    const double des[2] = {0.2, 0.4};          // de/dx = 4, 8
    const char* sup[3] = {"all", "wall band", "interior"};

    // ======================= (a) symmetry of H and H^-1 =======================
    if (id == 0) printf("(a) symmetry  rel = |<a,Hb> - <Ha,b>| / (|a| |Hb|),  NX=%ld K=%ld dx=%g, P=%d\n", NX, K, dx, p);
    for (int walls = 1; walls >= 0; walls--)
    for (int c = 0; c < 3; c++)                                // lattice: Bx (X), By (Y), Bz (C)
    for (int t = 0; t < 2; t++)
    {
        const double de = des[t];
        int yw  = walls ? (c == 1 ? 2 : 1) : 0;
        int hz  = walls ? (c == 1 ? HZ_DIR : HZ_NEU) : HZ_PER;
        const char* nm[3] = {"Bx/X ", "By/Y ", "Bz/C "};
        for (int s = 0; s < 3; s++)
        {
            F a, b, Ha, Hb, Ia, Ib;
            fill(a, 11 + 100*c + 10*t, yw, s);
            fill(b, 13 + 100*c + 10*t, yw, s);
            helmholtz_apply_walls(a.a, 0, Ha.a, 0, hz, de, dx, dy);
            helmholtz_apply_walls(b.a, 0, Hb.a, 0, hz, de, dx, dy);
            helmholtz_solve_walls(a.a, 0, Ia.a, 0, hz, de, dx, dy, id, p);
            helmholtz_solve_walls(b.a, 0, Ib.a, 0, hz, de, dx, dy, id, p);
            double nrm = sqrt(dot(a,a,yw)*dot(Hb,Hb,yw));
            double rH  = fabs(dot(a,Hb,yw) - dot(Ha,b,yw))/nrm;
            double nI  = sqrt(dot(a,a,yw)*dot(Ib,Ib,yw));
            double rI  = fabs(dot(a,Ib,yw) - dot(Ia,b,yw))/nI;

            // per-cell Green residual a L b - b L a - div G, G_{k+1/2} = (a_k b_{k+1} - b_k a_{k+1})/h,
            // and the wall flux G at the wall half-points from the ghosts
            double rg_w = 0.0, rg_i = 0.0, lmax = 0.0, gw = 0.0;
            long k0 = 2, k1 = (yw == 2) ? K : K+1;
            auto Gx = [&](long j, long k) { return (a(j,k)*b(j+1,k) - b(j,k)*a(j+1,k))/dx; };
            auto Gy = [&](long j, long k) { return (a(j,k)*b(j,k+1) - b(j,k)*a(j,k+1))/dy; };
            for (long j = 2; j < J+2; j++)
            {
                for (long k = k0; k <= k1; k++)
                {
                    double La = (a(j+1,k)-2*a(j,k)+a(j-1,k))/(dx*dx) + (a(j,k+1)-2*a(j,k)+a(j,k-1))/(dy*dy);
                    double Lb = (b(j+1,k)-2*b(j,k)+b(j-1,k))/(dx*dx) + (b(j,k+1)-2*b(j,k)+b(j,k-1))/(dy*dy);
                    double lhs = a(j,k)*Lb - b(j,k)*La;
                    double dG = (Gx(j,k) - Gx(j-1,k))/dx + (Gy(j,k) - Gy(j,k-1))/dy;
                    lmax = std::max(lmax, fabs(a(j,k)*Lb) + fabs(b(j,k)*La));
                    bool w = (k <= k0+2 || k >= k1-2);
                    (w ? rg_w : rg_i) = std::max(w ? rg_w : rg_i, fabs(lhs - dG));
                }
                if (walls)
                {
                    // wall half-points: below row k0 (between k0-1 and k0) and above row k1
                    gw = std::max(gw, std::max(fabs(Gy(j, k0-1)), fabs(Gy(j, k1))));
                }
            }
            rg_w = gmax(rg_w)/gmax(lmax); rg_i = gmax(rg_i)/gmax(lmax); gw = gmax(gw)/gmax(lmax)*dy;
            if (id == 0)
                printf("  %s %s de/dx=%g %-9s  H: %.1e  H^-1: %.1e   Green/cell wall rows %.1e interior %.1e   G_wall %.1e\n",
                       walls ? "walls" : "per  ", nm[c], de/dx, sup[s], rH, rI, rg_w, rg_i, gw);
            if (id == 0) g_worst = std::max(g_worst, std::max(std::max(rH, rI), std::max(std::max(rg_w, rg_i), gw)));   // NEW (2026-11-06)
        }
    }

    // ============ (b) product rule and (4) local field-energy law, walls ============
    if (id == 0) printf("\n(b) product rule B.(CE) = div S + avg(J.E), and (4) the local law with the inertia flux\n");
    for (int walls = 1; walls >= 0; walls--)
    for (int t = 0; t < 2; t++)
    {
        const double de = des[t];
        const int yX = walls ? 1 : 0, yY = walls ? 2 : 0, yC = walls ? 1 : 0, yN = walls ? 2 : 0;
        const int hX = walls ? HZ_NEU : HZ_PER, hY = walls ? HZ_DIR : HZ_PER, hC = walls ? HZ_NEU : HZ_PER;
        F Bx, By, Bz, Ez, Ey, Ex;
        fill(Bx, 201 + t, yX, 0); fill(By, 202 + t, yY, 0); fill(Bz, 203 + t, yC, 0);
        fill(Ez, 204 + t, yN, 0);                    // nodes, k = 1..K+1, 0 on wall nodes (C_flux_em)
        fill(Ey, 205 + t, yX, 0);                    // x-faces
        fill(Ex, 206 + t, yY, 0);                    // y-faces, 0 on wall faces (C_flux_em)

        // C E (so that Bdot' = -C E, cf. Cx_em = -dt CE_x etc.)
        F CEx, CEy, CEz;
        for (long j = 2; j < J+2; j++)
        {
            for (long k = 2; k <= K+1; k++)
            {
                CEx(j,k) = (Ez(j,k) - Ez(j,k-1))/dy;
                CEz(j,k) = (Ey(j,k) - Ey(j-1,k))/dx - (Ex(j,k) - Ex(j,k-1))/dy;
            }
            for (long k = 1; k <= K+1; k++)
                CEy(j,k) = -(Ez(j,k) - Ez(j-1,k))/dx;
        }
        boundary_conditions_ct(CEx.a, id, p, yX); boundary_conditions_ct(CEy.a, id, p, yY); boundary_conditions_ct(CEz.a, id, p, yC);

        // lattice fluxes and node/face J.E
        auto F1y = [&](long j, long k) { return  Ez(j,k)*0.5*(Bx(j,k) + Bx(j,k+1)); };     // X lattice, y
        auto F2x = [&](long j, long k) { return -Ez(j,k)*0.5*(By(j,k) + By(j+1,k)); };     // Y lattice, x
        auto F3x = [&](long j, long k) { return  Ey(j,k)*0.5*(Bz(j,k) + Bz(j+1,k)); };     // C, x
        auto F3y = [&](long j, long k) { return -Ex(j,k)*0.5*(Bz(j,k) + Bz(j,k+1)); };     // C, y
        auto qa  = [&](long j, long k) { return -(Bx(j,k+1) - Bx(j,k))/dy*Ez(j,k); };      // node
        auto qb  = [&](long j, long k) { return  (By(j+1,k) - By(j,k))/dx*Ez(j,k); };      // node
        auto qc  = [&](long j, long k) { return -(Bz(j+1,k) - Bz(j,k))/dx*Ey(j,k); };      // x-face, Jy Ey
        auto qd  = [&](long j, long k) { return  (Bz(j,k+1) - Bz(j,k))/dy*Ex(j,k); };      // y-face, Jx Ex

        // X lattice point (j,k): div and avg(JE) ; Y lattice point (j,k) ; C cell
        auto PX  = [&](long j, long k) { return Bx(j,k)*CEx(j,k); };
        auto PY  = [&](long j, long k) { return (walls && (k == 1 || k == K+1)) ? 0.0 : By(j,k)*CEy(j,k); };
        auto dX  = [&](long j, long k) { return (F1y(j,k) - F1y(j,k-1))/dy; };
        auto dY  = [&](long j, long k) { return (walls && (k == 1 || k == K+1)) ? 0.0 : (F2x(j,k) - F2x(j-1,k))/dx; };
        auto jX  = [&](long j, long k) { return 0.5*(qa(j,k) + qa(j,k-1)); };
        auto jY  = [&](long j, long k) { return (walls && (k == 1 || k == K+1)) ? 0.0 : 0.5*(qb(j,k) + qb(j-1,k)); };

        // (4): bdot = H^-1 (-C E), direct wdot = 1/2 (bdot.Hb + b.H bdot) on each lattice
        F mX, mY, mZ, dbx, dby, dbz, HBx, HBy, HBz, Hdx, Hdy, Hdz;
        for (long j = 0; j < J+4; j++) for (long k = 0; k < K+4; k++) { mX(j,k) = -CEx(j,k); mY(j,k) = -CEy(j,k); mZ(j,k) = -CEz(j,k); }
        helmholtz_solve_walls(mX.a, 0, dbx.a, 0, hX, de, dx, dy, id, p);
        helmholtz_solve_walls(mY.a, 0, dby.a, 0, hY, de, dx, dy, id, p);
        helmholtz_solve_walls(mZ.a, 0, dbz.a, 0, hC, de, dx, dy, id, p);
        boundary_conditions_ct(dbx.a, id, p, yX); boundary_conditions_ct(dby.a, id, p, yY); boundary_conditions_ct(dbz.a, id, p, yC);
        helmholtz_apply_walls(Bx.a, 0, HBx.a, 0, hX, de, dx, dy);  helmholtz_apply_walls(dbx.a, 0, Hdx.a, 0, hX, de, dx, dy);
        helmholtz_apply_walls(By.a, 0, HBy.a, 0, hY, de, dx, dy);  helmholtz_apply_walls(dby.a, 0, Hdy.a, 0, hY, de, dx, dy);
        helmholtz_apply_walls(Bz.a, 0, HBz.a, 0, hC, de, dx, dy);  helmholtz_apply_walls(dbz.a, 0, Hdz.a, 0, hC, de, dx, dy);
        boundary_conditions_ct(HBx.a, id, p, yX); boundary_conditions_ct(Hdx.a, id, p, yX);   // seam ghosts for the cell average
        boundary_conditions_ct(HBy.a, id, p, yY); boundary_conditions_ct(Hdy.a, id, p, yY);
        boundary_conditions_ct(HBz.a, id, p, yC); boundary_conditions_ct(Hdz.a, id, p, yC);
        auto wd  = [&](F& b, F& d, F& Hb, F& Hd, long j, long k) { return 0.5*(d(j,k)*Hb(j,k) + b(j,k)*Hd(j,k)); };
        // inertia flux G(bdot, b) on a lattice: div G = bdot L b - b L bdot
        auto dG  = [&](F& f, F& g, long j, long k) {
            auto gx = [&](long jj, long kk) { return (f(jj,kk)*g(jj+1,kk) - g(jj,kk)*f(jj+1,kk))/dx; };
            auto gy = [&](long jj, long kk) { return (f(jj,kk)*g(jj,kk+1) - g(jj,kk)*f(jj,kk+1))/dy; };
            return (gx(j,k) - gx(j-1,k))/dx + (gy(j,k) - gy(j,k-1))/dy; };
        const double h2 = 0.5*de*de;

        double rb_w = 0, rb_i = 0, sc_b = 0, r4_w = 0, r4_i = 0, sc_4 = 0;
        double sBCE = 0, sJE = 0, sW = 0, sD = 0;
        for (long j = 2; j < J+2; j++)
            for (long k = 2; k <= K+1; k++)
            {
                double lhs = 0.5*(PX(j-1,k) + PX(j,k)) + 0.5*(PY(j,k-1) + PY(j,k)) + Bz(j,k)*CEz(j,k);
                double dS  = 0.5*(dX(j-1,k) + dX(j,k)) + 0.5*(dY(j,k-1) + dY(j,k))
                           + (F3x(j,k) - F3x(j-1,k))/dx + (F3y(j,k) - F3y(j,k-1))/dy;
                double JE  = 0.5*(jX(j-1,k) + jX(j,k)) + 0.5*(jY(j,k-1) + jY(j,k))
                           + 0.5*(qc(j,k) + qc(j-1,k)) + 0.5*(qd(j,k) + qd(j,k-1));
                bool w = (k <= 4 || k >= K-1);
                double r = fabs(lhs - dS - JE);
                (w ? rb_w : rb_i) = std::max(w ? rb_w : rb_i, r);
                sc_b = std::max(sc_b, fabs(lhs) + fabs(dS) + fabs(JE));
                sBCE += lhs; sJE += JE;

                // (4)
                bool yw0 = walls && (k == K+1);             // Y face K+1 is a wall face (b = bdot = 0)
                double wdir = 0.5*(wd(Bx,dbx,HBx,Hdx,j-1,k) + wd(Bx,dbx,HBx,Hdx,j,k))
                            + 0.5*((walls && k-1 == 1 ? 0.0 : wd(By,dby,HBy,Hdy,j,k-1)) + (yw0 ? 0.0 : wd(By,dby,HBy,Hdy,j,k)))
                            + wd(Bz,dbz,HBz,Hdz,j,k);
                double dGi = 0.5*(dG(dbx,Bx,j-1,k) + dG(dbx,Bx,j,k))
                           + 0.5*((walls && k-1 == 1 ? 0.0 : dG(dby,By,j,k-1)) + (yw0 ? 0.0 : dG(dby,By,j,k)))
                           + dG(dbz,Bz,j,k);
                double r4 = fabs(wdir + dS + JE + h2*dGi);
                (w ? r4_w : r4_i) = std::max(w ? r4_w : r4_i, r4);
                sc_4 = std::max(sc_4, fabs(wdir) + fabs(dS) + fabs(JE) + fabs(h2*dGi));
                sW += wdir; sD += h2*dGi;
            }
        // adjoint <B, CE> vs <C^T B, E> as lattice sums
        double sL = 0, sR = 0;
        for (long j = 2; j < J+2; j++)
        {
            for (long k = 2; k <= K+1; k++) sL += Bx(j,k)*CEx(j,k) + Bz(j,k)*CEz(j,k) + qc(j,k);
            for (long k = (walls ? 2 : 1); k <= (walls ? K : K); k++) sL += By(j,k)*CEy(j,k);
            for (long k = 1; k <= (walls ? K+1 : K); k++) sR += qa(j,k) + qb(j,k) + qd(j,k);
        }
        // periodic y: Y and node rows 1..K are the K distinct rows; walls: Y unknowns 2..K, nodes 1..K+1
        sL = gsum(sL); sR = gsum(sR);
        double sR_full = sR + 0.0;  // qc already moved to sL side: <B,CE> = sum(J.E) means sL_noqc = sR + sum qc
        (void) sR_full;
        rb_w = gmax(rb_w)/gmax(sc_b); rb_i = gmax(rb_i)/gmax(sc_b);
        r4_w = gmax(r4_w)/gmax(sc_4); r4_i = gmax(r4_i)/gmax(sc_4);
        sBCE = gsum(sBCE); sJE = gsum(sJE); sW = gsum(sW); sD = gsum(sD);
        if (id == 0)
        {
            printf("  %s de/dx=%g  (b) per cell: wall rows %.1e interior %.1e;  global |sum B.CE - sum J.E|/sum|.| = %.1e\n",
                   walls ? "walls" : "per  ", de/dx, rb_w, rb_i, fabs(sBCE - sJE)/(fabs(sBCE) + fabs(sJE)));
            printf("  %s de/dx=%g  (4) per cell: wall rows %.1e interior %.1e;  global sum(div G) = %.1e (rel. to sum|wdot| %.1e)\n",
                   walls ? "walls" : "per  ", de/dx, r4_w, r4_i, sD, fabs(sW));
            // NEW (2026-11-06): (b) and (4) residuals; sum(div G) relative to sum|wdot|
            g_worst = std::max(g_worst, std::max(std::max(rb_w, rb_i), fabs(sBCE - sJE)/(fabs(sBCE) + fabs(sJE))));
            g_worst = std::max(g_worst, std::max(std::max(r4_w, r4_i), fabs(sD)/fabs(sW)));
        }
    }
    // NEW (2026-11-06): pass/fail for make test (every metric is a relative round-off residual)
    const double tol = 1e-12;
    if (id == 0) printf("%s: worst relative residual %.1e (tolerance %.0e) -> %s\n", "test_energy_identities", g_worst, tol, g_worst < tol ? "PASS" : "FAIL");
    MPI_Finalize();
    // ORIGINAL (pre-2026-11-06):
    // return 0;
    return (id == 0 && !(g_worst < tol)) ? 1 : 0;
}
