#ifndef HAVE_Z_SOLVE_H
#define HAVE_Z_SOLVE_H

#include <omp.h>


template<typename ValueT, typename SizeT>
void z_solve_impl(
  ValueT* rho_i_ptr,
  ValueT* qs_ptr,
  ValueT* square_ptr,
  ValueT* u_ptr,
  ValueT* rhs_ptr,
  SizeT j,
  SizeT nx,
  SizeT nxmax,
  SizeT ny,
  SizeT nz);

template<typename SizeT, typename ValueT>
void z_solve(
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

//---------------------------------------------------------------------
//     Performs line solves in Z direction by first factoring
//     the block-tridiagonal matrix into an upper triangular matrix,
//     and then performing back substitution to solve for the unknow
//     vectors of each line.
//
//     Make sure we treat elements zero to cell_size in the direction
//     of the sweep.
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//---------------------------------------------------------------------

      //if (timeron) call timer_start(t_zsolve)

//---------------------------------------------------------------------
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     This function computes the left hand side for the three z-factors
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     Compute the indices for storing the block-diagonal matrix;
//     determine c (labeled f) and s jacobians
//---------------------------------------------------------------------

  // do j = 1, ny-2
  // use taskloop without dependencies, there is no overlap possible with other places
#pragma oss task for out(u_ptr[0])
  for (SizeT j = 1; j <= ny-2; ++j) {
    z_solve_impl(rho_i_ptr, qs_ptr, square_ptr, u_ptr, rhs_ptr, j, nx, nxmax, ny, nz);
  } // enddo

  return;

}

#endif // HAVE_Z_SOLVE_H
