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

  y_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz);

  // wait for all tasks to complete
#pragma omp taskwait

  z_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz);

  // wait for all tasks to complete
#pragma omp taskwait

  add(u, rhs, nx, nxmax, ny, nz);
}

#endif // HAVE_ADI_H
