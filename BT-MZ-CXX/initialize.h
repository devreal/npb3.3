#ifndef HAVE_INITIALIZE_H
#define HAVE_INITIALIZE_H

#include "NArray.h"

#include "exact_solution.h"
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
*/
template<typename ValueT>
void initialize(ValueT* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  ArrayViewT4<ValueT, size_type> u{u_ptr, 5, nxmax, ny, nz};
/*
c---------------------------------------------------------------------
c  Later (in compute_rhs) we compute 1/u for every element. A few of
c  the corner elements are not used, but it convenient (and faster)
c  to compute the whole thing with a simple loop. Make sure those
c  values are nonzero by initializing the whole thing here.
c---------------------------------------------------------------------
 */
#pragma omp parallel shared(dnxm1,nx,dnym1,ny,dnzm1,nz, u)
{
  size_type i, j, k;
  ValueT xi, eta, zeta;

  //StaticArrayT3<ValueT, size_type, 5, 3, 2, 1, 1, 1> Pface{};
  //StaticArrayT1<ValueT, size_type, 5, 1> temp{};

  // do k = 0, nz-1
#pragma omp for
  for (size_type k = 0; k <= nz-1; ++k) {
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
  } // end do

#pragma omp for
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
    ValueT zeta = dble(k) * dnzm1;
    // do j = 0, ny-1
    for (size_type j = 0; j <= ny-1; ++j) {
      ValueT eta = dble(j) * dnym1;
      // do i = 0, nx-1
      for (size_type i = 0; i <= nx-1; ++i) {
        StaticArrayT3<ValueT, size_type, 5, 3, 2, 1, 1, 1> Pface{};
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
  } // enddo


//---------------------------------------------------------------------
//     now store the exact values on the boundaries
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     west face
//---------------------------------------------------------------------
  i = 0;
  xi = 0.0;
#pragma omp for
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
    ValueT zeta = dble(k) * dnzm1;
    //do j = 0, ny-1
    for (size_type j = 0; j <= ny-1; ++j) {
      ValueT eta = dble(j) * dnym1;
      StaticArrayT1<ValueT, size_type, 5, 1> temp{};
      exact_solution(xi, eta, zeta, temp.begin());
      //do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = temp(m);
      } // enddo
    } // enddo

  } // enddo

//---------------------------------------------------------------------
//     east face
//---------------------------------------------------------------------

  i = nx-1;
  xi = 1.0;
#pragma omp for
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
    ValueT zeta = dble(k) * dnzm1;
    //do j = 0, ny-1
    for (size_type j = 0; j <= ny-1; ++j) {
      ValueT eta = dble(j) * dnym1;
      StaticArrayT1<ValueT, size_type, 5, 1> temp{};
      exact_solution(xi, eta, zeta, temp.begin());
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
        u(m,i,j,k) = temp(m);
      } // enddo
    } // enddo
  } // enddo

//---------------------------------------------------------------------
//     south face
//---------------------------------------------------------------------
  j = 0;
  eta = 0.0;
#pragma omp for
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
    ValueT zeta = dble(k) * dnzm1;
    // do i = 0, nx-1
    for (size_type i = 0; i <= nx-1; ++i) {
      ValueT xi = dble(i) * dnxm1;
      StaticArrayT1<ValueT, size_type, 5, 1> temp{};
      exact_solution(xi, eta, zeta, temp.begin());
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = temp(m);
      } // enddo
    } // enddo
  } // enddo


//---------------------------------------------------------------------
//     north face
//---------------------------------------------------------------------
  j = ny-1;
  eta = 1.0;
#pragma omp for
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
      ValueT zeta = dble(k) * dnzm1;
      StaticArrayT1<ValueT, size_type, 5, 1> temp{};
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

//---------------------------------------------------------------------
//     bottom face
//---------------------------------------------------------------------
  k = 0;
  zeta = 0.0;
#pragma omp for
  // do j = 0, ny-1
  for (size_type j = 0; j <= ny-1; ++j) {
    ValueT eta = dble(j) * dnym1;
    StaticArrayT1<ValueT, size_type, 5, 1> temp{};
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

//---------------------------------------------------------------------
//     top face
//---------------------------------------------------------------------
  k = nz-1;
  zeta = 1.0;
#pragma omp for
  // do j = 0, ny-1
  for (size_type j = 0; j <= ny-1; ++j) {
    ValueT eta = dble(j) * dnym1;
    StaticArrayT1<ValueT, size_type, 5, 1> temp{};
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

} // omp parallel

}


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
  //DO(n, 1, 5, {
  for (int n = 1; n <= 5; ++n) {
      //do n = 1, 5
    //DO(m, 1, 5, {
    for (int m = 1; m <= 5; ++m) {
        lhs(m,n,1,0) = 0.0;
        lhs(m,n,2,0) = 0.0;
        lhs(m,n,3,0) = 0.0;
        lhs(m,n,1,i) = 0.0;
        lhs(m,n,2,i) = 0.0;
        lhs(m,n,3,i) = 0.0;
    } // ); // enddo
  } // ); // enddo

//---------------------------------------------------------------------
//     next, set all diagonal values to 1. This is overkill, but convenient
//---------------------------------------------------------------------
  //do m = 1, 5
  //DO(m, 1, 5, {
  for (int m = 1; m <= 5; ++m) {
    lhs(m,m,2,0) = 1.0;
    lhs(m,m,2,i) = 1.0;
  } // ); // enddo

}


#endif // HAVE_INITIALIZE_H
