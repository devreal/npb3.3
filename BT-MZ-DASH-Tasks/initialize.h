#ifndef HAVE_INITIALIZE_H
#define HAVE_INITIALIZE_H

#include "NArray.h"

#include "exact_solution.h"

template<typename ValueT>
void initialize(ValueT* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<ValueT, size_type> u{u_ptr, 5, nxmax, ny, nz};
//#pragma omp parallel shared(dnxm1,nx,dnym1,ny,dnzm1,nz)
{
  size_type i, j, k;
  ValueT xi, eta, zeta;

  // do k = 0, nz-1
//#pragma omp for
  //for (size_type k = 0; k <= nz-1; ++k) {
  dash::tasks::taskloop(0, nz-1+1, [&](auto k_begin, auto k_end) {
    for (size_type k = k_begin; k < k_end; ++k) {
      //do j = 0, ny-1
      for (size_type j = 0; j <= ny-1; ++j) {
        // do i = 0, nx-1
        for (size_type i = 0; i <= nx-1; ++i) {
          // do m = 1, 5
          for (size_type m = 1; m <= 5; ++m) {
            u(m,i,j,k) = 1.0;
          } // end do
        } // end do
      } //end do
    }
  }); // end do

  dash::tasks::complete(true);

//#pragma omp for
  //do k = 0, nz-1
  //for (size_type k = 0; k <= nz-1; ++k) {
  dash::tasks::taskloop(0, nz-1+1, [&](auto k_begin, auto k_end) {
    StaticArrayT3<ValueT, size_type, 5, 3, 2, 1, 1, 1> Pface{};
    for (size_type k = k_begin; k < k_end; ++k) {
      ValueT zeta = dble(k) * dnzm1;
      // do j = 0, ny-1
      for (size_type j = 0; j <= ny-1; ++j) {
        ValueT eta = dble(j) * dnym1;
        // do i = 0, nx-1
        for (size_type i = 0; i <= nx-1; ++i) {
          ValueT xi = dble(i) * dnxm1;

          // do ix = 1, 2
          for  (size_type ix = 1; ix <= 2; ++ix) {
            exact_solution(dble(ix-1), eta, zeta, &Pface(1,1,ix));
          } // enddo

          // do iy = 1, 2
          for  (size_type iy = 1; iy <= 2; ++iy) {
            exact_solution(xi, dble(iy-1) , zeta, &Pface(1,2,iy));
          } // enddo

          // do iz = 1, 2
          for  (size_type iz = 1; iz <= 2; ++iz) {
            exact_solution(xi, eta, dble(iz-1), &Pface(1,3,iz));
          } // enddo

          // do m = 1, 5
          for (int m = 1; m <= 5; ++m) {
            ValueT Pxi   = xi   * Pface(m,1,2) +
                      (1.0-xi)   * Pface(m,1,1);
            ValueT Peta  = eta  * Pface(m,2,2) +
                      (1.0-eta)  * Pface(m,2,1);
            ValueT Pzeta = zeta * Pface(m,3,2) +
                      (1.0-zeta) * Pface(m,3,1);

            u(m,i,j,k) = Pxi + Peta + Pzeta -
                      Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
                      Pxi*Peta*Pzeta;

          } // enddo
        } // enddo
      } // enddo
    }// enddo
  });

  dash::tasks::complete(true);


//---------------------------------------------------------------------
//     now store the exact values on the boundaries
//---------------------------------------------------------------------

  // NOTE: no synchronization is needed for the tasks initializing the faces

//---------------------------------------------------------------------
//     west face
//---------------------------------------------------------------------
  i = 0;
  xi = 0.0;
//#pragma omp for
  //do k = 0, nz-1
  //for (size_type k = 0; k <= nz-1; ++k) {
  dash::tasks::taskloop(0, nz-1+1, [&, i, xi](auto k_begin, auto k_end) {
    StaticArrayT1<ValueT, size_type, 5, 1> temp{};
    for (size_type k = k_begin; k < k_end; ++k) {
      ValueT zeta = dble(k) * dnzm1;
      //do j = 0, ny-1
      for (size_type j = 0; j <= ny-1; ++j) {
        ValueT eta = dble(j) * dnym1;
        exact_solution(xi, eta, zeta, temp.begin());
        //do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
            u(m,i,j,k) = temp(m);
        } // enddo
      } // enddo
    } // enddo
  });

  dash::tasks::complete(true);

//---------------------------------------------------------------------
//     east face
//---------------------------------------------------------------------

  i = nx-1;
  xi = 1.0;
//#pragma omp for
  //do k = 0, nz-1
  //for (size_type k = 0; k <= nz-1; ++k) {
  dash::tasks::taskloop(0, nz-1+1, [&, i, xi](auto k_begin, auto k_end) {
    StaticArrayT1<ValueT, size_type, 5, 1> temp{};
    for (size_type k = k_begin; k < k_end; ++k) {
      ValueT zeta = dble(k) * dnzm1;
      //do j = 0, ny-1
      for (size_type j = 0; j <= ny-1; ++j) {
        ValueT eta = dble(j) * dnym1;
        exact_solution(xi, eta, zeta, temp.begin());
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = temp(m);
        } // enddo
      } // enddo
    } // enddo
  });

  dash::tasks::complete(true);

//---------------------------------------------------------------------
//     south face
//---------------------------------------------------------------------
  j = 0;
  eta = 0.0;
//#pragma omp for
  //do k = 0, nz-1
  //for (size_type k = 0; k <= nz-1; ++k) {
  dash::tasks::taskloop(0, nz-1+1, [&, j, eta](auto k_begin, auto k_end) {
    StaticArrayT1<ValueT, size_type, 5, 1> temp{};
    for (size_type k = k_begin; k < k_end; ++k) {
      ValueT zeta = dble(k) * dnzm1;
      // do i = 0, nx-1
      for (size_type i = 0; i <= nx-1; ++i) {
        ValueT xi = dble(i) * dnxm1;
        exact_solution(xi, eta, zeta, temp.begin());
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
            u(m,i,j,k) = temp(m);
        } // enddo
      } // enddo
    } // enddo
  });
  dash::tasks::complete(true);


//---------------------------------------------------------------------
//     north face
//---------------------------------------------------------------------
  j = ny-1;
  eta = 1.0;
//#pragma omp for
  //do k = 0, nz-1
  //for (size_type k = 0; k <= nz-1; ++k) {
  dash::tasks::taskloop(0, nz-1+1, [&, j, eta](auto k_begin, auto k_end) {
    StaticArrayT1<ValueT, size_type, 5, 1> temp{};
    for (size_type k = k_begin; k < k_end; ++k) {
      ValueT zeta = dble(k) * dnzm1;
      // do i = 0, nx-1
      for (size_type i = 0; i <= nx-1; ++i) {
        ValueT xi = dble(i) * dnxm1;
        exact_solution(xi, eta, zeta, temp.begin());
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = temp(m);
        }// enddo
      }// enddo
    } // enddo
  });
  dash::tasks::complete(true);


//---------------------------------------------------------------------
//     bottom face
//---------------------------------------------------------------------
  k = 0;
  zeta = 0.0;
//#pragma omp for
  // do j = 0, ny-1
  //for (size_type j = 0; j <= ny-1; ++j) {
  dash::tasks::taskloop(0, ny-1+1, [&, k, zeta](auto j_begin, auto j_end) {
    StaticArrayT1<ValueT, size_type, 5, 1> temp{};
    for (size_type j = j_begin; j < j_end; ++j) {
      ValueT eta = dble(j) * dnym1;
      // do i =0, nx-1
      for (size_type i = 0; i <= nx-1; ++i) {
        ValueT xi = dble(i) *dnxm1;
        exact_solution(xi, eta, zeta, temp.begin());
        //do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
            u(m,i,j,k) = temp(m);
        } // enddo
      } // enddo
    } // enddo
  });
  dash::tasks::complete(true);


//---------------------------------------------------------------------
//     top face
//---------------------------------------------------------------------
  k = nz-1;
  zeta = 1.0;
//#pragma omp for
  // do j = 0, ny-1
  //for (size_type j = 0; j <= ny-1; ++j) {
  dash::tasks::taskloop(0, ny-1+1, [&, k, zeta](auto j_begin, auto j_end) {
    StaticArrayT1<ValueT, size_type, 5, 1> temp{};
    for (size_type j = j_begin; j < j_end; ++j) {
      ValueT eta = dble(j) * dnym1;
      //do i =0, nx-1
      for (size_type i = 0; i <= nx-1; ++i) {
        ValueT xi = dble(i) * dnxm1;
        exact_solution(xi, eta, zeta, temp.begin());
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = temp(m);
        } // enddo
      } // enddo
    } // enddo
  });

} // omp parallel


  dash::tasks::complete(true);

  return;
}
/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  initialize(u, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     This subroutine initializes the field variable u using
c     tri-linear transfinite interpolation of the boundary values
c---------------------------------------------------------------------

      include 'header.h'

      integer nx, nxmax, ny, nz
      double precision u(5,0:nxmax-1,0:ny-1,0:nz-1)
      integer i, j, k, m, ix, iy, iz
      double precision  xi, eta, zeta, Pface(5,3,2), Pxi, Peta,
     >     Pzeta, temp(5)

c---------------------------------------------------------------------
c  Later (in compute_rhs) we compute 1/u for every element. A few of
c  the corner elements are not used, but it convenient (and faster)
c  to compute the whole thing with a simple loop. Make sure those
c  values are nonzero by initializing the whole thing here.
c---------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(temp,Pzeta,Peta,Pxi,m,iz,iy,
!$OMP& Pface,ix,xi,i,eta,j,zeta,k)
!$OMP&  SHARED(dnxm1,nx,dnym1,ny,dnzm1,nz)
!$OMP DO SCHEDULE(STATIC)
      do k = 0, nz-1
         do j = 0, ny-1
            do i = 0, nx-1
               do m = 1, 5
                  u(m,i,j,k) = 1.0
               end do
            end do
         end do
      end do
!$OMP END DO nowait
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c     first store the "interpolated" values everywhere on the zone
c---------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      do k = 0, nz-1
         zeta = dble(k) * dnzm1
         do j = 0, ny-1
            eta = dble(j) * dnym1
            do i = 0, nx-1
               xi = dble(i) * dnxm1

               do ix = 1, 2
                  call exact_solution(dble(ix-1), eta, zeta,
     >                    Pface(1,1,ix))
               enddo

               do iy = 1, 2
                  call exact_solution(xi, dble(iy-1) , zeta,
     >                    Pface(1,2,iy))
               enddo

               do iz = 1, 2
                  call exact_solution(xi, eta, dble(iz-1),
     >                    Pface(1,3,iz))
               enddo

               do m = 1, 5
                  Pxi   = xi   * Pface(m,1,2) +
     >                    (1.0d0-xi)   * Pface(m,1,1)
                  Peta  = eta  * Pface(m,2,2) +
     >                    (1.0d0-eta)  * Pface(m,2,1)
                  Pzeta = zeta * Pface(m,3,2) +
     >                    (1.0d0-zeta) * Pface(m,3,1)

                  u(m,i,j,k) = Pxi + Peta + Pzeta -
     >                    Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
     >                    Pxi*Peta*Pzeta

               enddo
            enddo
         enddo
      enddo
!$OMP END DO nowait

c---------------------------------------------------------------------
c     now store the exact values on the boundaries
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     west face
c---------------------------------------------------------------------
      i = 0
      xi = 0.0d0
!$OMP DO SCHEDULE(STATIC)
      do k = 0, nz-1
         zeta = dble(k) * dnzm1
         do j = 0, ny-1
            eta = dble(j) * dnym1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo
!$OMP END DO nowait

c---------------------------------------------------------------------
c     east face
c---------------------------------------------------------------------

      i = nx-1
      xi = 1.0d0
!$OMP DO SCHEDULE(STATIC)
      do k = 0, nz-1
         zeta = dble(k) * dnzm1
         do j = 0, ny-1
            eta = dble(j) * dnym1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo
!$OMP END DO nowait

c---------------------------------------------------------------------
c     south face
c---------------------------------------------------------------------
      j = 0
      eta = 0.0d0
!$OMP DO SCHEDULE(STATIC)
      do k = 0, nz-1
         zeta = dble(k) * dnzm1
         do i = 0, nx-1
            xi = dble(i) * dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo
!$OMP END DO nowait


c---------------------------------------------------------------------
c     north face
c---------------------------------------------------------------------
      j = ny-1
      eta = 1.0d0
!$OMP DO SCHEDULE(STATIC)
      do k = 0, nz-1
         zeta = dble(k) * dnzm1
         do i = 0, nx-1
            xi = dble(i) * dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo
!$OMP END DO

c---------------------------------------------------------------------
c     bottom face
c---------------------------------------------------------------------
      k = 0
      zeta = 0.0d0
!$OMP DO SCHEDULE(STATIC)
      do j = 0, ny-1
         eta = dble(j) * dnym1
         do i =0, nx-1
            xi = dble(i) *dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo
!$OMP END DO nowait

c---------------------------------------------------------------------
c     top face
c---------------------------------------------------------------------
      k = nz-1
      zeta = 1.0d0
!$OMP DO SCHEDULE(STATIC)
      do j = 0, ny-1
         eta = dble(j) * dnym1
         do i =0, nx-1
            xi = dble(i) * dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo
!$OMP END DO nowait
!$OMP END PARALLEL

      return
      end

*/



//---------------------------------------------------------------------
//---------------------------------------------------------------------

template<typename ValueT, typename SizeT>
void lhsinit(ValueT* lhs_ptr,
             SizeT   size)
{
  NArrayView<ValueT, 4, NArrayShape<4, SizeT, NARRAY_MEMORDER_COL, 1, 1, 1, 0>>
    lhs(lhs_ptr, 5, 5, 3, size+1);
  SizeT i = size;

//---------------------------------------------------------------------
//     zero the whole left hand side for starters
//---------------------------------------------------------------------
  //do m = 1, 5
  DO(n, 1, 5, {
      //do n = 1, 5
    DO(m, 1, 5, {
        lhs(m,n,1,0) = 0.0;
        lhs(m,n,2,0) = 0.0;
        lhs(m,n,3,0) = 0.0;
        lhs(m,n,1,i) = 0.0;
        lhs(m,n,2,i) = 0.0;
        lhs(m,n,3,i) = 0.0;
    }); // enddo
  }); // enddo

//---------------------------------------------------------------------
//     next, set all diagonal values to 1. This is overkill, but convenient
//---------------------------------------------------------------------
  //do m = 1, 5
  DO(m, 1, 5, {
    lhs(m,m,2,0) = 1.0;
    lhs(m,m,2,i) = 1.0;
  }); // enddo

}

/*
      subroutine lhsinit(lhs, size)
      implicit none
      integer size
      double precision lhs(5,5,3,0:size)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer i, m, n

      i = size
c---------------------------------------------------------------------
c     zero the whole left hand side for starters
c---------------------------------------------------------------------
      do m = 1, 5
         do n = 1, 5
            lhs(m,n,1,0) = 0.0d0
            lhs(m,n,2,0) = 0.0d0
            lhs(m,n,3,0) = 0.0d0
            lhs(m,n,1,i) = 0.0d0
            lhs(m,n,2,i) = 0.0d0
            lhs(m,n,3,i) = 0.0d0
         enddo
      enddo

c---------------------------------------------------------------------
c     next, set all diagonal values to 1. This is overkill, but convenient
c---------------------------------------------------------------------
      do m = 1, 5
         lhs(m,m,2,0) = 1.0d0
         lhs(m,m,2,i) = 1.0d0
      enddo

      return
      end
*/

#endif // HAVE_INITIALIZE_H
