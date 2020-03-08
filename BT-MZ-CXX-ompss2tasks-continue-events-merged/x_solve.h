#ifndef HAVE_X_SOLVE_H
#define HAVE_X_SOLVE_H


template<typename ValueT, typename SizeT>
void x_solve_impl(
  ValueT* rho_i_ptr,
  ValueT* qs_ptr,
  ValueT* square_ptr,
  ValueT* u_ptr,
  ValueT* rhs_ptr,
  SizeT k,
  SizeT nx,
  SizeT nxmax,
  SizeT ny,
  SizeT nz);

/**
c---------------------------------------------------------------------
c
c     Performs line solves in X direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix,
c     and then performing back substitution to solve for the unknow
c     vectors of each line.
c
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c
c---------------------------------------------------------------------
*/

template<typename SizeT, typename ValueT>
void x_solve(
  ValueT* rho_i_ptr,
  ValueT* qs_ptr,
  ValueT* square_ptr,
  ValueT* u_ptr,
  ValueT* rhs_ptr,
  SizeT nx,
  SizeT nxmax,
  SizeT ny,
  SizeT nz)
{

/*
c---------------------------------------------------------------------
c     This function computes the left hand side in the xi-direction
c---------------------------------------------------------------------

*/

  /*
c---------------------------------------------------------------------
c     determine a (labeled f) and n jacobians
c---------------------------------------------------------------------
  */
      //do k = 1, nz-2
  for (SizeT k = 1; k <= nz-2; k++) {
#pragma oss task depend(out: rhs_ptr[k])
{
    x_solve_impl(rho_i_ptr, qs_ptr, square_ptr, u_ptr, rhs_ptr, k, nx, nxmax, ny, nz);
} // pragma omp task
  }

  return;
}

#endif // HAVE_X_SOLVE_H
