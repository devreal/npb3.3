
#include "npbparams_cc.h"

#include "rhs.h"

template<typename ValueT, typename SizeT>
void compute_rhs_parallel(
  ValueT* rho_i,
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
    compute_rhs(rho_i, us,
                vs, ws,
                qs, square,
                rhs, forcing, u,
                nx, nxmax, ny, nz);
#pragma oss taskwait
}


/**
 * Instantiate template
 */
template void compute_rhs_parallel(double*, double*, double*, double*, double*,
                                  double*, double*, double*, double*,
                                  size_type, size_type, size_type, size_type);
