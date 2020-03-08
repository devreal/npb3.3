#ifndef HAVE_ADI_H
#define HAVE_ADI_H

#include "NArray.h"

#include "x_solve.h"
#include "y_solve.h"
#include "z_solve.h"
#include "add.h"
#include "rhs.h"

template<typename ValueT, typename SizeT>
void adi(ValueT* rho_i,
         ValueT* us,
         ValueT* vs,
         ValueT* ws,
         ValueT* qs,
         ValueT* square,
         ValueT* rhs,
         ValueT* forcing,
         ValueT* u,
         SizeT nx, SizeT nxmax, SizeT ny, SizeT nz)
{
  compute_rhs(rho_i, us, vs, ws, qs, square, rhs,
              forcing, u, nx, nxmax, ny, nz);

  x_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz);
  // TODO: this taskwait is needed for some reason, maybe dependencies are broken?
#pragma omp taskwait

  y_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz);

  // wait for all tasks to complete
#pragma omp taskwait

  z_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz);

  // wait for all tasks to complete
#pragma omp taskwait

  add(u, rhs, nx, nxmax, ny, nz);
}

/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  adi(rho_i, us, vs, ws, qs, square, rhs,
     $                forcing, u, nx, nxmax, ny, nz)

      integer nx, nxmax, ny, nz
      double precision rho_i  (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 us     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 vs     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 ws     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 qs     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 square (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 rhs    (5,0:nxmax-1,0:ny-1,0:nz-1),
     $                 forcing(5,0:nxmax-1,0:ny-1,0:nz-1),
     $                 u      (5,0:nxmax-1,0:ny-1,0:nz-1)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      call compute_rhs(rho_i, us, vs, ws, qs, square, rhs,
     $                 forcing, u, nx, nxmax, ny, nz)

      call x_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz)

      call y_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz)

      call z_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz)

      call add(u, rhs, nx, nxmax, ny, nz)

      return
      end

*/

#endif // HAVE_ADI_H
