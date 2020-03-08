#ifndef HAVE_Y_SOLVE_H
#define HAVE_Y_SOLVE_H

template<typename ValueT, typename SizeT>
void y_solve_impl(
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

template<typename SizeT, typename ValueT>
void y_solve(
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
c     Performs line solves in Y direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix,
c     and then performing back substitution to solve for the unknow
c     vectors of each line.
c
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c---------------------------------------------------------------------
*/

/*
c---------------------------------------------------------------------
c     This function computes the left hand side for the three y-factors
c---------------------------------------------------------------------
*/

/*
c---------------------------------------------------------------------
c     Compute the indices for storing the tri-diagonal matrix;
c     determine a (labeled f) and n jacobians for cell c
c---------------------------------------------------------------------
*/

  for (SizeT k = 1; k <= nz-2; ++k) {
#pragma oss task depend(out: rhs_ptr[k])
{
    y_solve_impl(rho_i_ptr, qs_ptr, square_ptr, u_ptr, rhs_ptr, k, nx, nxmax, ny, nz);
} // pragma omp task
  } // enddo

  return;

}

#endif // HAVE_Y_SOLVE_H
