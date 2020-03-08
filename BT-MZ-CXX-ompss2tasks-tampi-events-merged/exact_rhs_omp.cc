
#include <iostream>
#include <omp.h>

#include "npbparams_cc.h"


void exact_rhs_1(
  double* forcing_ptr,
  size_type k,
  size_type nx,
  size_type nxmax,
  size_type ny,
  size_type nz);

void exact_rhs_2(
  double* forcing_ptr,
  size_type k,
  size_type nx,
  size_type nxmax,
  size_type ny,
  size_type nz);

void exact_rhs_3(
  double* forcing_ptr,
  size_type k,
  size_type nx,
  size_type nxmax,
  size_type ny,
  size_type nz);

void exact_rhs_4(
  double* forcing_ptr,
  size_type j,
  size_type nx,
  size_type nxmax,
  size_type ny,
  size_type nz);

void exact_rhs_5(
  double* forcing_ptr,
  size_type k,
  size_type nx,
  size_type nxmax,
  size_type ny,
  size_type nz);

/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine exact_rhs(forcing, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     compute the right hand side based on exact solution
c---------------------------------------------------------------------
*/
void exact_rhs(
  double* forcing_ptr,
  size_type nx,
  size_type nxmax,
  size_type ny,
  size_type nz)
{

/*
c---------------------------------------------------------------------
c     initialize
c---------------------------------------------------------------------
*/
#pragma oss task for
  for (int k= 0; k <= nz-1; ++k) {
    exact_rhs_1(forcing_ptr, k, nx, nxmax, ny, nz);
  } // enddo
#pragma oss taskwait


/*
c---------------------------------------------------------------------
c     xi-direction flux differences
c---------------------------------------------------------------------
*/
#pragma oss task for
  for (int k = 1; k <= nz-2; ++k) {
    exact_rhs_2(forcing_ptr, k, nx, nxmax, ny, nz);
  } //enddo
#pragma oss taskwait


/*
c---------------------------------------------------------------------
c     eta-direction flux differences
c---------------------------------------------------------------------
*/
#pragma oss task for
  for (int k = 1; k <= nz-2; ++k) {
    exact_rhs_3(forcing_ptr, k, nx, nxmax, ny, nz);
  } // enddo
#pragma oss taskwait


/*
c---------------------------------------------------------------------
c     zeta-direction flux differences
c---------------------------------------------------------------------
*/
#pragma oss task for
  for (int j=1; j <= ny-2; ++j) {
    exact_rhs_4(forcing_ptr, j, nx, nxmax, ny, nz);
  } // enddo
#pragma oss taskwait


/*
c---------------------------------------------------------------------
c     now change the sign of the forcing function,
c---------------------------------------------------------------------
*/
//!$OMP DO SCHEDULE(STATIC)
#pragma oss task for
  for (int k = 1; k <= nz-2; ++k) {
    exact_rhs_5(forcing_ptr, k, nx, nxmax, ny, nz);
  } // enddo
#pragma oss taskwait

  return;
}
