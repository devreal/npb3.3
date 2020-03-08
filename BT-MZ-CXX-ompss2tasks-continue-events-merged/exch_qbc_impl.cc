
#include <iostream>
#include <mpi.h>

#include "NArray.h"
#include "header.h"

#include "exch_qbc.h"

/* 3D column-major Array, 0-index-based */
using QBCArrayT = NArrayView<double, 3, NArrayShape<3, size_type, NARRAY_MEMORDER_COL, 1, 1, 1>>;
/* 4D column-major Array, 0-index-based except on the fastest dimension */
using UArrayT = NArrayView<double, 4, NArrayShape<4, size_type, NARRAY_MEMORDER_COL, 1, 0, 0, 0>>;


void
copy_y_face_impl(
  double* u_ptr,
  double* qbc_ptr,
  int k,
  int nx,
  int nxmax,
  int ny,
  int nz,
  int jloc,
  int dir)
{
  UArrayT u(u_ptr, 5, nxmax, ny, nz);
  QBCArrayT qbc(qbc_ptr, 5, nx-2, nz-2);

  int j = jloc;
  if (dir == DIR_IN) {
    for (int i = 1; i <= nx-2; i++) {
      for (int m = 1; m <= 5; ++m) {
        u(m,i,j,k) = qbc(m,i,k);
      } // end do
    } // end do
  } else if (dir == DIR_OUT) {
    for (int i = 1; i <= nx-2; i++) {
      for (int m = 1; m <= 5; ++m) {
        qbc(m,i,k) = u(m,i,j,k);
      } // end do
    } // end do
  } else {
    std::cout << "Erroneous data designation: " << dir << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 6);
  } // endif

  return;
}

void
copy_x_face_impl(
  double* u_ptr,
  double* qbc_ptr,
  int k,
  int nx,
  int nxmax,
  int ny,
  int nz,
  int iloc,
  int dir)
{
  UArrayT u(u_ptr, 5, nxmax, ny, nz);
  QBCArrayT qbc(qbc_ptr, 5, ny-2, nz-2);

  int i = iloc;
  if (dir == DIR_IN) {
    for (int j = 1; j <= ny-2; ++j) {
      for (int m = 1; m <= 5; ++m) {
        u(m,i,j,k) = qbc(m,j,k);
      } // end do
    } // end do
  } else if (dir == DIR_OUT) {
    for (int j = 1; j <= ny-2; ++j) {
      for (int m = 1; m <= 5; ++m) {
        qbc(m,j,k) = u(m,i,j,k);
      } // end do
    } // end do
  } else {
    std::cout << "Erroneous data designation: " << dir << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 6);
  } // endif

  return;
}
