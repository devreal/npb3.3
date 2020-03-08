
#include "npbparams_cc.h"
#include "initialize.h"

void initialize(double* u_ptr, size_type nx, size_type nxmax, size_type ny, size_type nz)
{
  // do k = 0, nz-1
#pragma oss task for
  for (size_type k = 0; k <= nz-1; ++k) {
    initialize_1(k, u_ptr, nx, nxmax, ny, nz);
  } // end do
#pragma oss taskwait

#pragma oss task for 
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
    initialize_2(k, u_ptr, nx, nxmax, ny, nz);
  } // enddo

#pragma oss taskwait


//---------------------------------------------------------------------
//     now store the exact values on the boundaries
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     west face
//---------------------------------------------------------------------
#pragma oss task for
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
    initialize_3(k, u_ptr, nx, nxmax, ny, nz);
  } // enddo

#pragma oss taskwait

//---------------------------------------------------------------------
//     east face
//---------------------------------------------------------------------

#pragma oss task for
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
    initialize_4(k, u_ptr, nx, nxmax, ny, nz);
  } // enddo

#pragma oss taskwait

//---------------------------------------------------------------------
//     south face
//---------------------------------------------------------------------
#pragma oss task for
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
    initialize_5(k, u_ptr, nx, nxmax, ny, nz);
  } // enddo

#pragma oss taskwait


//---------------------------------------------------------------------
//     north face
//---------------------------------------------------------------------
#pragma oss task for
  //do k = 0, nz-1
  for (size_type k = 0; k <= nz-1; ++k) {
    initialize_6(k, u_ptr, nx, nxmax, ny, nz);
  } // enddo

#pragma oss taskwait

//---------------------------------------------------------------------
//     bottom face
//---------------------------------------------------------------------
#pragma oss task for
  // do j = 0, ny-1
  for (size_type j = 0; j <= ny-1; ++j) {
    initialize_7(j, u_ptr, nx, nxmax, ny, nz);
  } // enddo

#pragma oss taskwait

//---------------------------------------------------------------------
//     top face
//---------------------------------------------------------------------
#pragma oss task for
  // do j = 0, ny-1
  for (size_type j = 0; j <= ny-1; ++j) {
    initialize_8(j, u_ptr, nx, nxmax, ny, nz);
  } // enddo

#pragma oss taskwait

}

