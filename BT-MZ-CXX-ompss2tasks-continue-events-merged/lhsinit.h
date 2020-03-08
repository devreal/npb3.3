#ifndef HAVE_LHSINIT_H
#define HAVE_LHSINIT_H

#include "NArray.h"
#include "header.h"

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

#endif // HAVE_LHSINIT_H
