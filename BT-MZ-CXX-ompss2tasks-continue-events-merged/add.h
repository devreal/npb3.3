#ifndef HAVE_ADD_H
#define HAVE_ADD_H


template<typename ValueT, typename SizeT>
void add_impl(
  ValueT* u_ptr,
  ValueT* rhs_ptr,
  SizeT k,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz);

template<typename ValueT, typename SizeT>
void add(
  ValueT* u_ptr,
  ValueT* rhs_ptr,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz)
{

//---------------------------------------------------------------------
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     addition of update to the vector u
//---------------------------------------------------------------------

#pragma oss task for out(u_ptr[0])
  for (SizeT k = 1; k <= nz-2; ++k) {
    add_impl(u_ptr, rhs_ptr, k, nx, nxmax, ny, nz);
  }

      return;
}

#endif // HAVE_ADD_H

