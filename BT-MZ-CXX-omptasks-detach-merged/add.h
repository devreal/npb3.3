
#include "NArray.h"

template<typename ValueT, typename SizeT>
void add(
  ValueT* u_ptr,
  ValueT* rhs_ptr,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz)
{

  ArrayViewT4<ValueT, SizeT> u(u_ptr, 5, nxmax, ny, nz);
  ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr, 5, nxmax, ny, nz);

//---------------------------------------------------------------------
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     addition of update to the vector u
//---------------------------------------------------------------------
#pragma omp taskloop
  for (SizeT k = 1; k <= nz-2; ++k) {
    for (SizeT j = 1; j <= ny-2; ++j) {
      for (SizeT i = 1; i <= nx-2; ++i) {
        for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = u(m,i,j,k) + rhs(m,i,j,k);
        } // enddo
      } // enddo
    } // enddo
  } // enddo

      return;
}

