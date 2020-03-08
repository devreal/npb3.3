#ifndef HAVE_RHS_H
#define HAVE_RHS_H

#include "npbparams_cc.h"

/*
 * Outlined function declarations.
 */

void compute_rhs_1(
  double* rho_i_ptr,
  double* us_ptr,
  double* vs_ptr,
  double* ws_ptr,
  double* qs_ptr,
  double* square_ptr,
  double* u_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_2(
  double* rhs_ptr,
  double* forcing_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_3(
  double* rho_i_ptr,
  double* us_ptr,
  double* vs_ptr,
  double* ws_ptr,
  double* qs_ptr,
  double* square_ptr,
  double* rhs_ptr,
  double* u_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_4(
  double* rho_i_ptr,
  double* us_ptr,
  double* vs_ptr,
  double* ws_ptr,
  double* qs_ptr,
  double* square_ptr,
  double* rhs_ptr,
  double* u_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_5(
  double* rho_i_ptr,
  double* us_ptr,
  double* vs_ptr,
  double* ws_ptr,
  double* qs_ptr,
  double* square_ptr,
  double* rhs_ptr,
  double* u_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_6(
  double* rhs_ptr,
  double* u_ptr,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_7(
  double* rhs_ptr,
  double* u_ptr,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_8(
  double* rhs_ptr,
  double* u_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_9(
  double* rhs_ptr,
  double* u_ptr,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_10(
  double* rhs_ptr,
  double* u_ptr,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void compute_rhs_11(
  double* rhs_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);





template<typename ValueT, typename SizeT>
void compute_rhs(
  ValueT* rho_i_ptr,
  ValueT* us_ptr,
  ValueT* vs_ptr,
  ValueT* ws_ptr,
  ValueT* qs_ptr,
  ValueT* square_ptr,
  ValueT* rhs_ptr,
  ValueT* forcing_ptr,
  ValueT* u_ptr,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz)
{

  for (SizeT k = 0; k <= nz-1; ++k) {
#pragma oss task depend(out: square_ptr[k])
{
    compute_rhs_1(rho_i_ptr, us_ptr, vs_ptr, ws_ptr, qs_ptr,
                  square_ptr,
                  u_ptr, k, nx, nxmax, ny, nz);
} // pragma omp task
  } // enddo

//---------------------------------------------------------------------
// copy the exact forcing term to the right hand side;  because
// this forcing term is known, we can store it on the whole zone
// including the boundary
//---------------------------------------------------------------------

  for (SizeT k = 0; k <= nz-1; ++k) {
#pragma oss task depend(out: rhs_ptr[k])
{
    compute_rhs_2(rhs_ptr, forcing_ptr, k, nx, nxmax, ny, nz);
} // pragma omp task
  } // enddo


//---------------------------------------------------------------------
//     compute xi-direction fluxes
//---------------------------------------------------------------------
  for (SizeT k = 1; k <= nz-2; ++k) {
#pragma oss task depend(out: rhs_ptr[k]) depend(in:square_ptr[k])
{
    compute_rhs_3(rho_i_ptr, us_ptr, vs_ptr, ws_ptr, qs_ptr,
                  square_ptr, rhs_ptr,
                  u_ptr, k, nx, nxmax, ny, nz);
} // pragma omp task
  } // enddo

  for (SizeT k = 1; k <= nz-2; ++k) {
#pragma oss task depend(out: rhs_ptr[k]) depend(in:square_ptr[k])
{
    compute_rhs_4(rho_i_ptr, us_ptr, vs_ptr, ws_ptr, qs_ptr,
                  square_ptr, rhs_ptr,
                  u_ptr, k, nx, nxmax, ny, nz);
} // pragma omp parallel
  } // enddo

  for (SizeT k = 1; k <= nz-2; ++k) {
#pragma oss task depend(out: rhs_ptr[k]) depend(in:square_ptr[k], square_ptr[k-1], square_ptr[k+1])
{
    compute_rhs_5(rho_i_ptr, us_ptr, vs_ptr, ws_ptr, qs_ptr,
                  square_ptr, rhs_ptr,
                  u_ptr, k, nx, nxmax, ny, nz);
} // pragma omp task
  } // enddo

//---------------------------------------------------------------------
//     add fourth order zeta-direction dissipation
//---------------------------------------------------------------------
  SizeT k = 1;
#pragma oss task depend(out: rhs_ptr[k])
{
    compute_rhs_6(rhs_ptr, u_ptr, nx, nxmax, ny, nz);
} // pragma omp task

  k = 2;
#pragma oss task depend(out: rhs_ptr[k])
{
    compute_rhs_7(rhs_ptr, u_ptr, nx, nxmax, ny, nz);
} // pragma omp task

  for (SizeT k = 3; k <= nz-4; ++k) {
#pragma oss task depend(out: rhs_ptr[k])
{
    compute_rhs_8(rhs_ptr, u_ptr, k, nx, nxmax, ny, nz);
} // pragma omp task
  } // enddo

  k = nz-3;
#pragma oss task depend(out: rhs_ptr[k])
{
    compute_rhs_9(rhs_ptr, u_ptr, nx, nxmax, ny, nz);
} // pragma omp task

  k = nz-2;
#pragma oss task depend(out: rhs_ptr[k])
{
    compute_rhs_10(rhs_ptr, u_ptr, nx, nxmax, ny, nz);
} // pragma omp task

  for (SizeT k = 1; k <= nz-2; ++k) {
#pragma oss task depend(out: rhs_ptr[k])
{
    compute_rhs_11(rhs_ptr, k, nx, nxmax, ny, nz);
} // pragma omp task
  }


  return;
}



template<typename ValueT, typename SizeT>
extern
void compute_rhs_parallel(
  ValueT* rho_i,
  ValueT* us,
  ValueT* vs,
  ValueT* ws,
  ValueT* qs,
  ValueT* square,
  ValueT* rhs,
  ValueT* forcing,
  ValueT* u,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz);


#endif // HAVE_RHS_H
