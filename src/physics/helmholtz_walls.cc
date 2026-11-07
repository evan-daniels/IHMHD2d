////////////////////////////////////////////////////////////////////////////////
//
// helmholtz_walls.cc -- (1 - de^2 lap_h) b = q for one field component on the
// CT lattice, periodic in x, periodic or conducting walls in y (INERTIA2d).
//
// Adapted from ~/Research/QHD src/2d/samples/MHD2d/helmholtz_fft.cc (doubly
// periodic r2c version): same gather-to-rank-0 / serial FFTW / scatter layout
// and plan cache, with a real-to-real transform chosen by the wall parity.
//
// lap_h is the compact 5-point Laplacian on the component's own lattice. With
// the ghost parities of boundary_conditions_ct the field is the image of a
// field periodic over 2*Ly, so one separable r2r transform diagonalises lap_h:
//
//   kind  lattice / wall rule                  x         y            Ny
//   HZ_PER  any, periodic in y                 R2HC      R2HC         K
//   HZ_NEU  cell-centred in y, even mirror     R2HC      REDFT10      K      (Bx on x-faces, Bz)
//           (ywall = 1)                                  (DCT-II)
//   HZ_DIR  y-faces, 0 on wall faces 1, K+1,   R2HC      RODFT00      K-1    (By, faces 2..K)
//           odd mirror (ywall = 2)                       (DST-I)
//
//   b_hat = q_hat / (1 + de^2 [ (4/dx^2) sin^2(pi m/Nx) + (4/dy^2) sin^2(theta_n) ])
//   theta_n = pi m_y/K (periodic, m_y = halfcomplex frequency), pi n/(2K) (DCT-II, n = 0..K-1),
//             pi (n+1)/(2K) (DST-I, n = 0..K-2)
//
// Requires a uniform mesh (checked by the caller). Uses serial FFTW only
// (no fftw3_mpi). Ghost cells of b are NOT filled; call boundary_conditions_ct.
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"
#include <fftw3.h>
#include <vector>

using namespace std;

namespace {
struct HzPlan { long Nx = -1, Ny = -1; double* in = 0; double* out = 0; fftw_plan fwd, bwd; };
HzPlan hz_plans[3];
}

void CENTPACK::helmholtz_solve_walls(const doublearray3d& q, const int& lq, doublearray3d& b, const int& lb,
                                     const int& kind, const double& de, const double& dx, const double& dy,
                                     const int& id, const int& p)
{
    long J = q.getIndex1Size() - 4;
    long K = q.getIndex2Size() - 4;
    long Nx = J*p;
    long Ny = (kind == HZ_DIR) ? K - 1 : K;
    long k0 = 2;                                   // first unknown row (cells 2..K+1, or faces 2..K)
    long j, k;

    vector<double> loc(J*Ny), glb;
    for (j = 2; j < J+2; j++)
        for (k = 0; k < Ny; k++)
            loc[(j-2)*Ny + k] = q(j, k0 + k, lq);

    if (id == 0) glb.resize(Nx*Ny);
    // rank r owns global x-indices r*J .. r*J+J-1, so the gather is a concatenation
    MPI_Gather(loc.data(), J*Ny, MPI_DOUBLE, glb.data(), J*Ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (id == 0)
    {
        HzPlan& P = hz_plans[kind];
        if (P.Nx != Nx || P.Ny != Ny)
        {
            if (P.in)
            {
                fftw_destroy_plan(P.fwd);
                fftw_destroy_plan(P.bwd);
                fftw_free(P.in);
                fftw_free(P.out);
            }
            P.in  = fftw_alloc_real(Nx*Ny);
            P.out = fftw_alloc_real(Nx*Ny);
            fftw_r2r_kind fy = (kind == HZ_PER) ? FFTW_R2HC : (kind == HZ_NEU) ? FFTW_REDFT10 : FFTW_RODFT00;
            fftw_r2r_kind by = (kind == HZ_PER) ? FFTW_HC2R : (kind == HZ_NEU) ? FFTW_REDFT01 : FFTW_RODFT00;
            P.fwd = fftw_plan_r2r_2d(Nx, Ny, P.in, P.out, FFTW_R2HC, fy, FFTW_ESTIMATE);
            P.bwd = fftw_plan_r2r_2d(Nx, Ny, P.out, P.in, FFTW_HC2R, by, FFTW_ESTIMATE);
            P.Nx = Nx;
            P.Ny = Ny;
        }

        for (long i = 0; i < Nx*Ny; i++) P.in[i] = glb[i];
        fftw_execute(P.fwd);

        const double pi = acos(-1.0);
        const double norm = 1.0/((double) Nx * (kind == HZ_PER ? (double) K : 2.0*K));
        const double cx = 4.0*de*de/(dx*dx), cy = 4.0*de*de/(dy*dy);
        for (long i = 0; i < Nx; i++)
        {
            long m = (i <= Nx/2) ? i : Nx - i;
            double sx = sin(pi*m/Nx);
            for (long n = 0; n < Ny; n++)
            {
                double th;
                if (kind == HZ_PER) th = pi*((n <= K/2) ? n : K - n)/K;
                else if (kind == HZ_NEU) th = pi*n/(2.0*K);
                else th = pi*(n + 1)/(2.0*K);
                double sy = sin(th);
                P.out[i*Ny + n] *= norm/(1.0 + cx*sx*sx + cy*sy*sy);
            }
        }

        fftw_execute(P.bwd);
        for (long i = 0; i < Nx*Ny; i++) glb[i] = P.in[i];
    }

    MPI_Scatter(glb.data(), J*Ny, MPI_DOUBLE, loc.data(), J*Ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (j = 2; j < J+2; j++)
        for (k = 0; k < Ny; k++)
            b(j, k0 + k, lb) = loc[(j-2)*Ny + k];
}

// q = (1 - de^2 lap_h) b on the unknowns of `kind` (b ghosts must be filled);
// for HZ_DIR the wall faces 1 and K+1 are set to 0
void CENTPACK::helmholtz_apply_walls(const doublearray3d& b, const int& lb, doublearray3d& q, const int& lq,
                                     const int& kind, const double& de, const double& dx, const double& dy)
{
    long J = b.getIndex1Size() - 4;
    long K = b.getIndex2Size() - 4;
    long k1 = (kind == HZ_DIR) ? K : K + 1;

    for (long j = 2; j < J+2; j++)
    {
        for (long k = 2; k <= k1; k++)
        {
            double lap = (b(j+1,k,lb) - 2.0*b(j,k,lb) + b(j-1,k,lb))/(dx*dx)
                       + (b(j,k+1,lb) - 2.0*b(j,k,lb) + b(j,k-1,lb))/(dy*dy);
            q(j,k,lq) = b(j,k,lb) - de*de*lap;
        }
        if (kind == HZ_DIR)
        {
            q(j,1,lq)   = 0.0;
            q(j,K+1,lq) = 0.0;
        }
    }
}
