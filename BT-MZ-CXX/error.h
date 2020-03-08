
#include "NArray.h"

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
template<typename ValueT, typename SizeT>
void error_norm(
  ValueT* rms_ptr,
  ValueT* u_ptr,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz)
{
  StaticArrayViewT1<ValueT, SizeT, 5, 1> rms(rms_ptr);
  ArrayViewT4<ValueT, SizeT> u(u_ptr, 5, nxmax, ny, nz);
  //do m = 1, 5
  for (int m = 1; m <= 5; ++m) {
      rms(m) = 0.0;
  }

#pragma omp parallel
{
  StaticArrayT1<ValueT, SizeT, 5, 1> rms_loc;
  StaticArrayT1<ValueT, SizeT, 5, 1> u_exact;
  //do m=1,5
  for (int m = 1; m <= 5; ++m) {
    rms_loc(m)=0.0;
  }
//!$OMP DO
#pragma omp for nowait
  //do k = 0, nz-1
  for (SizeT k = 0; k <= nz-1; ++k) {
    ValueT zeta = dble(k) * dnzm1;
    // do j = 0, ny-1
    for (int j = 0; j <= ny-1; ++j) {
      ValueT eta = dble(j) * dnym1;
      //do i = 0, nx-1
      for (int i = 0; i <= nx-1; ++i) {
        ValueT xi = dble(i) * dnxm1;
        exact_solution(xi, eta, zeta, u_exact.begin());

        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          ValueT add = u(m,i,j,k)-u_exact(m);
          rms_loc(m) = rms_loc(m) + add*add;
        } // enddo
      } // enddo
    } // enddo
  } // enddo
  //do m=1,5
  for (int m = 1; m <= 5; ++m) {
#pragma omp atomic
    rms(m)+=rms_loc(m);
  } // enddo
} // omp parallel

  //do m = 1, 5
  for (int m = 1; m <= 5; ++m) {
    rms(m) = rms(m) / (dble(nz-2)*dble(ny-2)*dble(nx-2));
    rms(m) = dsqrt(rms(m));
  } // enddo

  return;

}



template<typename ValueT, typename SizeT>
void rhs_norm(
  ValueT* rms_ptr,
  ValueT* rhs_ptr,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz)
{
  StaticArrayViewT1<ValueT, SizeT, 5, 1> rms(rms_ptr);
  ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr, 5, nxmax, ny, nz);
  //do m = 1, 5
  for (int m = 1; m <= 5; ++m) {
      rms(m) = 0.0;
  } // enddo

#pragma omp parallel
{
  StaticArrayT1<ValueT, SizeT, 5, 1> rms_loc;
  // do m=1,5
  for (int m = 1; m <= 5; ++m) {
      rms_loc(m)=0.0;
  } // enddo
#pragma omp for
  // do k = 1, nz-2
  for (SizeT k = 1; k <= nz-2; ++k) {
    // do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      // do i = 1, nx-2
      for (SizeT i = 1; i <= nx-2; ++i) {
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          ValueT add = rhs(m,i,j,k);
          rms_loc(m) = rms_loc(m) + add*add;
        } // enddo
      } //enddo
    } // enddo
  } // enddo
  // do m=1,5,1
  for (int m = 1; m <= 5; ++m) {
#pragma omp atomic
    rms(m)+=rms_loc(m);
  } // enddo
} // omp parallel

  //do m = 1, 5
  for (int m = 1; m <= 5; ++m) {
    rms(m) = rms(m) / (dble(nz-2)*dble(ny-2)*dble(nx-2));
    rms(m) = dsqrt(rms(m));
  } // enddo

  return;
}

