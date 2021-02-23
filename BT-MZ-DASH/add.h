
#include "NArray.h"

template<typename ValueT, typename SizeT>
void add(
  ValueT *__restrict u_ptr,
  ValueT *__restrict rhs_ptr,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz)
{

  ArrayViewT4<ValueT, SizeT> u(u_ptr, 5, nxmax, ny, nz);
  ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr, 5, nxmax, ny, nz);

//---------------------------------------------------------------------
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     addition of update to the vector u
//---------------------------------------------------------------------

  //if (timeron) call timer_start(t_add)
//!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m,i,j,k)
//!$OMP&  SHARED(nx,ny,nz)
//      do     k = 1, nz-2
//         do     j = 1, ny-2
//            do     i = 1, nx-2
//               do    m = 1, 5
#pragma omp parallel for
  for (SizeT k = 1; k <= nz-2; ++k) {
    for (SizeT j = 1; j <= ny-2; ++j) {
      for (SizeT i = 1; i <= nx-2; ++i) {
        for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = u(m,i,j,k) + rhs(m,i,j,k);
        } // enddo
      } // enddo
    } // enddo
  } // enddo

// !$OMP END PARALLEL DO
      //if (timeron) call timer_stop(t_add)

      return;
}
/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  add(u, rhs, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     addition of update to the vector u
c---------------------------------------------------------------------

      include 'header.h'

      integer nx, nxmax, ny, nz
      double precision rhs(5,0:nxmax-1,0:ny-1,0:nz-1),
     $                 u  (5,0:nxmax-1,0:ny-1,0:nz-1)

      integer i, j, k, m

      if (timeron) call timer_start(t_add)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m,i,j,k)
!$OMP&  SHARED(nx,ny,nz)
      do     k = 1, nz-2
         do     j = 1, ny-2
            do     i = 1, nx-2
               do    m = 1, 5
                  u(m,i,j,k) = u(m,i,j,k) + rhs(m,i,j,k)
               enddo
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO
      if (timeron) call timer_stop(t_add)

      return
      end
*/
