#ifndef HAVE_ADI_H
#define HAVE_ADI_H

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

#pragma oss task wait depend(out:u[0])
{
  compute_rhs(rho_i, us, vs, ws, qs, square, rhs,
              forcing, u, nx, nxmax, ny, nz);

  x_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz);

  y_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz);
}

  // task for with out(u[0]) inside this call
  z_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz);

  // task for with out(u[0]) inside this call
  add(u, rhs, nx, nxmax, ny, nz);
}

#endif // HAVE_ADI_H
