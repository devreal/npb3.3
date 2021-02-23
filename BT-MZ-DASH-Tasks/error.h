
#include <libdash.h>

#include "NArray.h"

template<typename ValueT, typename SizeT>
void error_norm(
  ValueT* rms_ptr,
  ValueT* u_ptr,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz)
{
  StaticArrayViewT1<ValueT, SizeT, 5, 1> rms(rms_ptr);
  StaticArrayT1<dash::tasks::Combinator<ValueT>, SizeT, 5, 1> rms_loc;
  ArrayViewT4<ValueT, SizeT> u(u_ptr, 5, nxmax, ny, nz);
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

      include 'header.h'

      integer nx, nxmax, ny, nz
      double precision u(5,0:nxmax-1,0:ny-1,0:nz-1)

      integer i, j, k, m
      double precision xi, eta, zeta, u_exact(5), rms(5), add
      double precision rms_loc(5)

*/

  //do m = 1, 5
  //for (int m = 1; m <= 5; ++m) {
  //    rms(m) = 0.0;
  //}

//!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rms_loc,add,m,u_exact,xi,i,eta,
//!$OMP& j,zeta,k)
//!$OMP&  SHARED(dnxm1,nx,dnym1,ny,dnzm1,nz)
//#pragma omp parallel
{
  //do m=1,5
  //for (int m = 1; m <= 5; ++m) {
  //  rms_loc(m)=0.0;
  //}
//!$OMP DO
//#pragma omp for nowait
  //do k = 0, nz-1
  //for (SizeT k = 0; k <= nz-1; ++k) {
  dash::tasks::taskloop(0, nz-1+1, [&](SizeT k_begin, SizeT k_end) {
    StaticArrayT1<ValueT, SizeT, 5, 1> u_exact;
    for (SizeT k = k_begin; k < k_end; ++k) {
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
            rms_loc(m).local() = rms_loc(m).local() + add*add;
          } // enddo
        } // enddo
      } // enddo
    } // enddo
  });
// !$OMP END DO nowait
  //do m=1,5
//  for (int m = 1; m <= 5; ++m) {
//!$OMP ATOMIC
//#pragma omp atomic
//    rms(m)+=rms_loc(m);
//  } // enddo
} // omp parallel
//!$OMP END PARALLEL

  dash::tasks::complete(true);

  for (int m = 1; m <= 5; ++m) {
    rms(m) = rms_loc(m).reduce(dash::plus<ValueT>{});
  }

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
  StaticArrayT1<dash::tasks::Combinator<ValueT>, SizeT, 5, 1> rms_loc;
  ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr, 5, nxmax, ny, nz);
/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine rhs_norm(rms, rhs, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      include 'header.h'

      integer nx, nxmax, ny, nz
      double precision rhs(5,0:nxmax-1,0:ny-1,0:nz-1)

      integer i, j, k, m
      double precision rms(5), add
      double precision rms_loc(5)
*/
  //do m = 1, 5
  //for (int m = 1; m <= 5; ++m) {
  //    rms(m) = 0.0;
  //} // enddo

//!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rms_loc,add,m,i,j,k)
//!$OMP&  SHARED(nx,ny,nz)
//#pragma omp parallel
{
  // do m=1,5
  //for (int m = 1; m <= 5; ++m) {
  //    rms_loc(m)=0.0;
  //} // enddo
//!$OMP DO
//#pragma omp for
  // do k = 1, nz-2
  //for (SizeT k = 1; k <= nz-2; ++k) {
  dash::tasks::taskloop(0, nz-2+1, [&](SizeT k_begin, SizeT k_end) {
    StaticArrayT1<ValueT, SizeT, 5, 1> u_exact;
    for (SizeT k = k_begin; k < k_end; ++k) {
    // do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      // do i = 1, nx-2
      for (SizeT i = 1; i <= nx-2; ++i) {
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          ValueT add = rhs(m,i,j,k);
          rms_loc(m).local() = rms_loc(m).local() + add*add;
        } // enddo
      } //enddo
    } // enddo
  } // enddo
  });
//!$OMP END DO nowait
} // omp parallel
//!$OMP END PARALLEL

  dash::tasks::complete(true);

  for (int m = 1; m <= 5; ++m) {
    rms(m) = rms_loc(m).reduce(dash::plus<ValueT>{});
  }

  //do m = 1, 5
  for (int m = 1; m <= 5; ++m) {
    rms(m) = rms(m) / (dble(nz-2)*dble(ny-2)*dble(nx-2));
    rms(m) = dsqrt(rms(m));
  } // enddo

  return;
}

