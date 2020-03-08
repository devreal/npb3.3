
#include "NArray.h"

#include "npbparams_cc.h"

template<typename ValueT, typename SizeT>
void add_impl(
  ValueT* u_ptr,
  ValueT* rhs_ptr,
  SizeT k,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz)
{

  ArrayViewT4<ValueT, SizeT> u(u_ptr, 5, nxmax, ny, nz);
  ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr, 5, nxmax, ny, nz);

  for (SizeT j = 1; j <= ny-2; ++j) {
    for (SizeT i = 1; i <= nx-2; ++i) {
      for (int m = 1; m <= 5; ++m) {
        u(m,i,j,k) = u(m,i,j,k) + rhs(m,i,j,k);
      } // enddo
    } // enddo
  } // enddo

      return;
}

/* instantiate function template */
template void add_impl(
  double*, double*,
  size_type, size_type, size_type, size_type, size_type);


