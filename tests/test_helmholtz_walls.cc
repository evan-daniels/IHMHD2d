////////////////////////////////////////////////////////////////////////////////
//
// test_helmholtz_walls.cc -- unit test for helmholtz_solve_walls (INERTIA2d)
//
// Random parity-matched B' on the CT lattice: B'x, B'y = discrete curl of a random
// corner stream function psi (so D B' = 0 to round-off), psi = 0 on the wall corner
// rows (so B'y = 0 on the wall faces); B'z random at cell centres.
// Solve (1 - de^2 lap_h) B = B' per component, fill ghosts with boundary_conditions_ct,
// and report
//   residual = max|(1 - de^2 lap_h) B - B'| / max|B'|          per component
//   max|D B| * dx / max|B|  separately on the wall rows, the periodic-x seam columns
//                           and the interior
// for conducting walls in y (Bx: HZ_NEU, By: HZ_DIR, Bz: HZ_NEU) and for periodic y.
//
// Run: mpirun -np P ./test_helmholtz_walls      (P must divide NX)
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"
static double g_worst = 0.0;   // NEW (2026-11-06): largest relative residual (rank 0), for the exit status
#include <random>
#include <vector>

using namespace std;
using namespace CENTPACK;

static double gmax(double v)
{
    double g;
    MPI_Allreduce(&v, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return g;
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int id, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    const long NX = 128, K = 64;
    const long J = NX/p;
    const double dx = 0.1, dy = 0.1;
    const double des[3] = {0.05, 0.2, 1.0};

    for (int walls = 1; walls >= 0; walls--)
    for (int t = 0; t < 3; t++)
    {
        const double de = des[t];
        const int kx = walls ? HZ_NEU : HZ_PER, ky = walls ? HZ_DIR : HZ_PER, kz = walls ? HZ_NEU : HZ_PER;
        const int yx = walls ? 1 : 0, yy = walls ? 2 : 0, yz = walls ? 1 : 0;

        // global random corner stream function psi(jg, k), k = 1..K+1 (corner row k at y_{k+1/2})
        mt19937_64 rng(12345 + 7*t + 100*walls);
        uniform_real_distribution<double> U(-1.0, 1.0);
        vector<double> psi(NX*(K+2), 0.0), bz(NX*(K+4), 0.0);
        for (long jg = 0; jg < NX; jg++)
        {
            for (long k = 1; k <= K+1; k++) psi[jg*(K+2) + k] = U(rng);
            for (long k = 2; k <= K+1; k++) bz[jg*(K+4) + k] = U(rng);
            if (walls) { psi[jg*(K+2) + 1] = 0.0; psi[jg*(K+2) + K+1] = 0.0; }
            else psi[jg*(K+2) + 1] = psi[jg*(K+2) + K+1];        // periodic in y
        }
        auto PS = [&](long j, long k) { long jg = ((id*J + j - 2) % NX + NX) % NX; return psi[jg*(K+2) + k]; };

        doublearray3d Q1(J+4,K+4,1), Q2(J+4,K+4,1), Qz(J+4,K+4,1), B1(J+4,K+4,1), B2(J+4,K+4,1), Bz(J+4,K+4,1), R(J+4,K+4,1);
        for (long j = 0; j < J+4; j++)
            for (long k = 0; k < K+4; k++)
            { Q1(j,k,0) = Q2(j,k,0) = Qz(j,k,0) = B1(j,k,0) = B2(j,k,0) = Bz(j,k,0) = 0.0; }

        for (long j = 2; j < J+2; j++)
        {
            for (long k = 2; k <= K+1; k++)
            {
                Q1(j,k,0) = (PS(j,k) - PS(j,k-1))/dy;                     // B'x on x-face (j+1/2, k)
                Qz(j,k,0) = bz[(id*J + j - 2)*(K+4) + k];
            }
            for (long k = 1; k <= K+1; k++)
                Q2(j,k,0) = -(PS(j,k) - PS(j-1,k))/dx;                     // B'y on y-face (j, k+1/2)
        }
        boundary_conditions_ct(Q1, id, p, yx);
        boundary_conditions_ct(Q2, id, p, yy);
        boundary_conditions_ct(Qz, id, p, yz);

        helmholtz_solve_walls(Q1, 0, B1, 0, kx, de, dx, dy, id, p);
        helmholtz_solve_walls(Q2, 0, B2, 0, ky, de, dx, dy, id, p);
        helmholtz_solve_walls(Qz, 0, Bz, 0, kz, de, dx, dy, id, p);
        boundary_conditions_ct(B1, id, p, yx);
        boundary_conditions_ct(B2, id, p, yy);
        boundary_conditions_ct(Bz, id, p, yz);

        // residuals
        double res[3], qm[3] = {0, 0, 0}, bm = 0.0;
        doublearray3d* QQ[3] = {&Q1, &Q2, &Qz};
        doublearray3d* BB[3] = {&B1, &B2, &Bz};
        int kinds[3] = {kx, ky, kz};
        for (int c = 0; c < 3; c++)
        {
            helmholtz_apply_walls(*BB[c], 0, R, 0, kinds[c], de, dx, dy);
            double r = 0.0;
            long k1 = (kinds[c] == HZ_DIR) ? K : K+1;
            for (long j = 2; j < J+2; j++)
                for (long k = 2; k <= k1; k++)
                {
                    r = max(r, fabs(R(j,k,0) - (*QQ[c])(j,k,0)));
                    qm[c] = max(qm[c], fabs((*QQ[c])(j,k,0)));
                    if (c < 2) bm = max(bm, fabs((*BB[c])(j,k,0)));
                }
            res[c] = gmax(r)/gmax(qm[c]);
        }
        bm = gmax(bm);

        // div B and div B' by region
        double dw = 0, ds = 0, di_ = 0, dq = 0;
        for (long j = 2; j < J+2; j++)
            for (long k = 2; k <= K+1; k++)
            {
                double d  = fabs((B1(j,k,0) - B1(j-1,k,0))/dx + (B2(j,k,0) - B2(j,k-1,0))/dy)*dx/bm;
                double d2 = fabs((Q1(j,k,0) - Q1(j-1,k,0))/dx + (Q2(j,k,0) - Q2(j,k-1,0))/dy)*dx/bm;
                long jg = id*J + j - 2;
                bool seam = (jg == 0 || jg == NX-1), wall = (k == 2 || k == K+1);
                if (wall) dw = max(dw, d);
                if (seam) ds = max(ds, d);
                if (!wall && !seam) di_ = max(di_, d);
                dq = max(dq, d2);
            }
        dw = gmax(dw); ds = gmax(ds); di_ = gmax(di_); dq = gmax(dq);

        // wall faces of By must stay exactly 0
        double bw = 0.0;
        if (walls)
            for (long j = 2; j < J+2; j++) bw = max(bw, max(fabs(B2(j,1,0)), fabs(B2(j,K+1,0))));
        bw = gmax(bw);

        if (id == 0)
            printf("%s p=%d de=%.2f (de/dx=%.1f)  residual Bx %.1e By %.1e Bz %.1e | max|D B|dx/|B|: walls %.1e seam %.1e interior %.1e (D B' %.1e) | By on wall %.1e\n",
                   walls ? "walls   " : "periodic", p, de, de/dx, res[0], res[1], res[2], dw, ds, di_, dq, bw);
        // NEW (2026-11-06): residuals, div B / div B' (relative) and By on the wall (exactly 0)
        if (id == 0) g_worst = std::max({g_worst, res[0], res[1], res[2], dw, ds, di_, dq, bw});
    }

    // NEW (2026-11-06): pass/fail for make test (every metric is a relative round-off residual)
    const double tol = 1e-12;
    if (id == 0) printf("%s: worst relative residual %.1e (tolerance %.0e) -> %s\n", "test_helmholtz_walls", g_worst, tol, g_worst < tol ? "PASS" : "FAIL");
    MPI_Finalize();
    // ORIGINAL (pre-2026-11-06):
    // return 0;
    return (id == 0 && !(g_worst < tol)) ? 1 : 0;
}
