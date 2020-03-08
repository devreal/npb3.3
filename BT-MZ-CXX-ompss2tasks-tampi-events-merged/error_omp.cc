#ifndef HAVE_ERROR_H
#define HAVE_ERROR_H

#include <cmath>

#include "npbparams_cc.h"

void error_norm_impl(
  double* rms_ptr,
  double* u_ptr,
  double* rms_loc_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);


void rhs_norm_impl(
  double* rms_ptr,
  double* rhs_ptr,
  double* rms_loc_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine error_norm(rms, u, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     this function computes the norm of the difference between the
c     computed solution and the exact solution
c---------------------------------------------------------------------
*/
void error_norm(
  double* rms,
  double* u,
  size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  //do m = 1, 5
  for (int m = 0; m < 5; ++m) {
      rms[m] = 0.0;
  }

  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
#pragma oss task
{
    double rms_loc[5];
    //do m=1,5
    for (int m = 0; m < 5; ++m) {
      rms_loc[m]=0.0;
    }
    error_norm_impl(rms, u, rms_loc, k, nx, nxmax, ny, nz);

    //do m=1,5
    for (int m = 0; m < 5; ++m) {
  #pragma oss atomic
      rms[m]+=rms_loc[m];
    } // enddo
}
  } // enddo

#pragma oss taskwait

  //do m = 1, 5
  for (int m = 0; m < 5; ++m) {
    rms[m] = rms[m] / ((double)(nz-2)*(double)(ny-2)*(double)(nx-2));
    rms[m] = std::sqrt(rms[m]);
  } // enddo

  return;

}


void rhs_norm(
  double* rms,
  double* rhs,
  size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  //do m = 1, 5
  for (int m = 0; m < 5; ++m) {
      rms[m] = 0.0;
  } // enddo

  // do k = 1, nz-2
  for (size_type k = 1; k <= nz-2; ++k) {
#pragma oss task
{
    double rms_loc[5];
    // do m=1,5
    for (int m = 0; m < 5; ++m) {
        rms_loc[m]=0.0;
    } // enddo
    rhs_norm_impl(rms, rhs, rms_loc, k, nx, nxmax, ny, nz);

    // do m=1,5,1
    for (int m = 0; m < 5; ++m) {
  // !$OMP ATOMIC
  #pragma oss atomic
      rms[m]+=rms_loc[m];
    } // enddo
}
  } // enddo

#pragma oss taskwait

  //do m = 1, 5
  for (int m = 0; m < 5; ++m) {
    rms[m] = rms[m] / ((double)(nz-2)*(double)(ny-2)*(double)(nx-2));
    rms[m] = std::sqrt(rms[m]);
  } // enddo


  return;
}

#endif // HAVE_ERROR_H
