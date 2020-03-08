
#include "NArray.h"

#include "header.h"

#include "exact_solution.h"

void error_norm_impl(
  double* rms_ptr,
  double* u_ptr,
  double* rms_loc_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  StaticArrayViewT1<double, size_type, 5, 1> rms(rms_ptr);
  ArrayViewT4<double, size_type> u(u_ptr, 5, nxmax, ny, nz);
  StaticArrayViewT1<double, size_type, 5, 1> rms_loc{rms_loc_ptr};
  StaticArrayT1<double, size_type, 5, 1> u_exact;

  double zeta = dble(k) * dnzm1;
  // do j = 0, ny-1
  for (int j = 0; j <= ny-1; ++j) {
    double eta = dble(j) * dnym1;
    //do i = 0, nx-1
    for (int i = 0; i <= nx-1; ++i) {
      double xi = dble(i) * dnxm1;
      exact_solution(xi, eta, zeta, u_exact.begin());

      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
        double add = u(m,i,j,k)-u_exact(m);
        rms_loc(m) = rms_loc(m) + add*add;
      } // enddo
    } // enddo
  } // enddo
}


void rhs_norm_impl(
  double* rms_ptr,
  double* rhs_ptr,
  double* rms_loc_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  StaticArrayViewT1<double, size_type, 5, 1> rms(rms_ptr);
  ArrayViewT4<double, size_type> rhs(rhs_ptr, 5, nxmax, ny, nz);
  StaticArrayViewT1<double, size_type, 5, 1> rms_loc{rms_loc_ptr};

  // do j = 1, ny-2
  for (size_type j = 1; j <= ny-2; ++j) {
    // do i = 1, nx-2
    for (size_type i = 1; i <= nx-2; ++i) {
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
        double add = rhs(m,i,j,k);
        rms_loc(m) = rms_loc(m) + add*add;
      } // enddo
    } //enddo
  } // enddo

  return;
}

