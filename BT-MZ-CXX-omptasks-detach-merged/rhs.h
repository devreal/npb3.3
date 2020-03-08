#ifndef HAVE_RHS_H
#define HAVE_RHS_H

#include "NArray.h"

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

  ArrayViewT3<ValueT, SizeT> rho_i(rho_i_ptr, nxmax, ny, nz);
  ArrayViewT3<ValueT, SizeT> us(us_ptr, nxmax, ny, nz);
  ArrayViewT3<ValueT, SizeT> vs(vs_ptr, nxmax, ny, nz);
  ArrayViewT3<ValueT, SizeT> ws(ws_ptr, nxmax, ny, nz);
  ArrayViewT3<ValueT, SizeT> qs(qs_ptr, nxmax, ny, nz);
  ArrayViewT3<ValueT, SizeT> square(square_ptr, nxmax, ny, nz);
  ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr, 5, nxmax, ny, nz);
  ArrayViewT4<ValueT, SizeT> forcing(forcing_ptr, 5, nxmax, ny, nz);
  ArrayViewT4<ValueT, SizeT> u(u_ptr, 5, nxmax, ny, nz);

//---------------------------------------------------------------------
//     compute the reciprocal of density and the kinetic energy,
//---------------------------------------------------------------------
{
  // do k = 0, nz-1
  for (int k = 0; k <= nz-1; ++k) {
#pragma omp task depend(out: square_ptr[k]) firstprivate(rho_i_ptr, qs_ptr, square_ptr, u_ptr, rhs_ptr, nx, nxmax, ny, nz, k)
{
    ArrayViewT3<ValueT, SizeT> rho_i(rho_i_ptr, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> us(us_ptr, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> vs(vs_ptr, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> ws(ws_ptr, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> qs(qs_ptr, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> square(square_ptr, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr, 5, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> u(u_ptr, 5, nxmax, ny, nz);

    //do j = 0, ny-1
    for (SizeT j = 0; j <= ny-1; ++j) {
      // do i = 0, nx-1
      for (SizeT i = 0; i <= nx-1; ++i) {
          double rho_inv = 1.00/u(1,i,j,k);
          rho_i(i,j,k) = rho_inv;
          us(i,j,k) = u(2,i,j,k) * rho_inv;
          vs(i,j,k) = u(3,i,j,k) * rho_inv;
          ws(i,j,k) = u(4,i,j,k) * rho_inv;
          square(i,j,k)     = 0.50* (
                 u(2,i,j,k)*u(2,i,j,k) +
                 u(3,i,j,k)*u(3,i,j,k) +
                 u(4,i,j,k)*u(4,i,j,k) ) * rho_inv;
          qs(i,j,k) = square(i,j,k) * rho_inv;
      } // enddo
    } // enddo
} // pragma omp task
  } // enddo

//---------------------------------------------------------------------
// copy the exact forcing term to the right hand side;  because
// this forcing term is known, we can store it on the whole zone
// including the boundary
//---------------------------------------------------------------------

  // do k = 0, nz-1
  for (int k = 0; k <= nz-1; ++k) {
#pragma omp task depend(out: rhs_ptr[k]) firstprivate(rho_i_ptr, qs_ptr, square_ptr, u_ptr, rhs_ptr, nx, nxmax, ny, nz, k)
{
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr, 5, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> forcing(forcing_ptr, 5, nxmax, ny, nz);
    //do j = 0, ny-1
    for (SizeT j = 0; j <= ny-1; ++j) {
      // do i = 0, nx-1
      for (SizeT i = 0; i <= nx-1; ++i) {
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = forcing(m,i,j,k);
        } // enddo
      } // enddo
    } // enddo
} // pragma omp task
  } // enddo

//---------------------------------------------------------------------
//     compute xi-direction fluxes
//---------------------------------------------------------------------
  // do k = 1, nz-2
  for (int k = 1; k <= nz-2; ++k) {
#pragma omp task depend(out: rhs_ptr[k]) depend(in:square_ptr[k]) firstprivate(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz, k)
{
    //do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      // do i = 1, nx-2
      for (SizeT i = 1; i <= nx-2; ++i) {
        ValueT uijk = us(i,j,k);
        ValueT up1  = us(i+1,j,k);
        ValueT um1  = us(i-1,j,k);

        rhs(1,i,j,k) = rhs(1,i,j,k) + dx1tx1 *
                   (u(1,i+1,j,k) - 2.0*u(1,i,j,k) +
                   u(1,i-1,j,k)) -
                   tx2 * (u(2,i+1,j,k) - u(2,i-1,j,k));

        rhs(2,i,j,k) = rhs(2,i,j,k) + dx2tx1 *
                   (u(2,i+1,j,k) - 2.0*u(2,i,j,k) +
                   u(2,i-1,j,k)) +
                   xxcon2*con43 * (up1 - 2.0*uijk + um1) -
                   tx2 * (u(2,i+1,j,k)*up1 -
                   u(2,i-1,j,k)*um1 +
                   (u(5,i+1,j,k)- square(i+1,j,k)-
                   u(5,i-1,j,k)+ square(i-1,j,k))*
                   c2);

        rhs(3,i,j,k) = rhs(3,i,j,k) + dx3tx1 *
                   (u(3,i+1,j,k) - 2.0*u(3,i,j,k) +
                   u(3,i-1,j,k)) +
                   xxcon2 * (vs(i+1,j,k) - 2.0*vs(i,j,k) +
                   vs(i-1,j,k)) -
                   tx2 * (u(3,i+1,j,k)*up1 -
                   u(3,i-1,j,k)*um1);

        rhs(4,i,j,k) = rhs(4,i,j,k) + dx4tx1 *
                   (u(4,i+1,j,k) - 2.00*u(4,i,j,k) +
                   u(4,i-1,j,k)) +
                   xxcon2 * (ws(i+1,j,k) - 2.0*ws(i,j,k) +
                   ws(i-1,j,k)) -
                   tx2 * (u(4,i+1,j,k)*up1 -
                   u(4,i-1,j,k)*um1);

        rhs(5,i,j,k) = rhs(5,i,j,k) + dx5tx1 *
                   (u(5,i+1,j,k) - 2.0*u(5,i,j,k) +
                   u(5,i-1,j,k)) +
                   xxcon3 * (qs(i+1,j,k) - 2.0*qs(i,j,k) +
                   qs(i-1,j,k)) +
                   xxcon4 * (up1*up1 -       2.0*uijk*uijk +
                   um1*um1) +
                   xxcon5 * (u(5,i+1,j,k)*rho_i(i+1,j,k) -
                   2.0*u(5,i,j,k)*rho_i(i,j,k) +
                   u(5,i-1,j,k)*rho_i(i-1,j,k)) -
                   tx2 * ( (c1*u(5,i+1,j,k) -
                   c2*square(i+1,j,k))*up1 -
                   (c1*u(5,i-1,j,k) -
                   c2*square(i-1,j,k))*um1 );
      } // enddo
    } // enddo

//---------------------------------------------------------------------
//     add fourth order xi-direction dissipation
//---------------------------------------------------------------------

    //do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      SizeT i = 1;
      // do m = 1, 5
      for (SizeT m = 1; m <= 5; ++m) {
        rhs(m,i,j,k) = rhs(m,i,j,k)- dssp *
                  ( 5.0*u(m,i,j,k) - 4.0*u(m,i+1,j,k) +
                  u(m,i+2,j,k));
      } // enddo

      i = 2;
      //do m = 1, 5
      for (SizeT m = 1; m <= 5; ++m) {
        rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                  (-4.0*u(m,i-1,j,k) + 6.0*u(m,i,j,k) -
                  4.0*u(m,i+1,j,k) + u(m,i+2,j,k));
      } // enddo
    } // enddo

    // do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      // do i = 3,nx-4
      for (SizeT i = 3; i <= nx-4; ++i) {
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                    (  u(m,i-2,j,k) - 4.0*u(m,i-1,j,k) +
                    6.0*u(m,i,j,k) - 4.0*u(m,i+1,j,k) +
                    u(m,i+2,j,k) );
        } // enddo
      } // enddo
    } // enddo

    //do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      SizeT i = nx-3;
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
        rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                  ( u(m,i-2,j,k) - 4.0*u(m,i-1,j,k) +
                  6.0*u(m,i,j,k) - 4.0*u(m,i+1,j,k) );
      } // enddo

      i = nx-2;
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
        rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                  ( u(m,i-2,j,k) - 4.0*u(m,i-1,j,k) +
                  5.0*u(m,i,j,k) );
      } // enddo
    } // enddo
} // pragma omp task
  } // enddo

//---------------------------------------------------------------------
//     compute eta-direction fluxes
//---------------------------------------------------------------------
  // do k = 1, nz-2
  for (int k = 1; k <= nz-2; ++k) {
#pragma omp task depend(out: rhs_ptr[k]) depend(in:square_ptr[k]) firstprivate(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz, k)
{
    // do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      // do i = 1, nx-2
      for (SizeT i = 1; i <= nx-2; ++i) {
          ValueT vijk = vs(i,j,k);
          ValueT vp1  = vs(i,j+1,k);
          ValueT vm1  = vs(i,j-1,k);
          rhs(1,i,j,k) = rhs(1,i,j,k) + dy1ty1 *
                  (u(1,i,j+1,k) - 2.0*u(1,i,j,k) +
                  u(1,i,j-1,k)) -
                  ty2 * (u(3,i,j+1,k) - u(3,i,j-1,k));
          rhs(2,i,j,k) = rhs(2,i,j,k) + dy2ty1 *
                  (u(2,i,j+1,k) - 2.0*u(2,i,j,k) +
                  u(2,i,j-1,k)) +
                  yycon2 * (us(i,j+1,k) - 2.0*us(i,j,k) +
                  us(i,j-1,k)) -
                  ty2 * (u(2,i,j+1,k)*vp1 -
                  u(2,i,j-1,k)*vm1);
          rhs(3,i,j,k) = rhs(3,i,j,k) + dy3ty1 *
                  (u(3,i,j+1,k) - 2.0*u(3,i,j,k) +
                  u(3,i,j-1,k)) +
                  yycon2*con43 * (vp1 - 2.0*vijk + vm1) -
                  ty2 * (u(3,i,j+1,k)*vp1 -
                  u(3,i,j-1,k)*vm1 +
                  (u(5,i,j+1,k) - square(i,j+1,k) -
                  u(5,i,j-1,k) + square(i,j-1,k))
                  *c2);
          rhs(4,i,j,k) = rhs(4,i,j,k) + dy4ty1 *
                  (u(4,i,j+1,k) - 2.0*u(4,i,j,k) +
                  u(4,i,j-1,k)) +
                  yycon2 * (ws(i,j+1,k) - 2.0*ws(i,j,k) +
                  ws(i,j-1,k)) -
                  ty2 * (u(4,i,j+1,k)*vp1 -
                  u(4,i,j-1,k)*vm1);
          rhs(5,i,j,k) = rhs(5,i,j,k) + dy5ty1 *
                  (u(5,i,j+1,k) - 2.0*u(5,i,j,k) +
                  u(5,i,j-1,k)) +
                  yycon3 * (qs(i,j+1,k) - 2.0*qs(i,j,k) +
                  qs(i,j-1,k)) +
                  yycon4 * (vp1*vp1       - 2.0*vijk*vijk +
                  vm1*vm1) +
                  yycon5 * (u(5,i,j+1,k)*rho_i(i,j+1,k) -
                  2.0*u(5,i,j,k)*rho_i(i,j,k) +
                  u(5,i,j-1,k)*rho_i(i,j-1,k)) -
                  ty2 * ((c1*u(5,i,j+1,k) -
                  c2*square(i,j+1,k)) * vp1 -
                  (c1*u(5,i,j-1,k) -
                  c2*square(i,j-1,k)) * vm1);
      } // enddo
    } // enddo

//---------------------------------------------------------------------
//     add fourth order eta-direction dissipation
//---------------------------------------------------------------------
    SizeT j = 1;
    // do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
        rhs(m,i,j,k) = rhs(m,i,j,k)- dssp *
                  ( 5.0*u(m,i,j,k) - 4.0*u(m,i,j+1,k) +
                    u(m,i,j+2,k));
      } // enddo
    } // enddo

    j = 2;
    // do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      //do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
        rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                  (-4.0*u(m,i,j-1,k) + 6.0*u(m,i,j,k) -
                    4.0*u(m,i,j+1,k) + u(m,i,j+2,k));
      } // enddo
    } // enddo

    // do j = 3, ny-4
    for (SizeT j = 3; j <= ny-4; ++j) {
      // do i = 1,nx-2
      for (SizeT i = 1; i <= nx-2; ++i) {
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                    (  u(m,i,j-2,k) - 4.0*u(m,i,j-1,k) +
                     6.0*u(m,i,j,k) - 4.0*u(m,i,j+1,k) +
                      u(m,i,j+2,k) );
        } // enddo
      } // enddo
    } // enddo

    j = ny-3;
    // do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                    ( u(m,i,j-2,k) - 4.0*u(m,i,j-1,k) +
                    6.0*u(m,i,j,k) - 4.0*u(m,i,j+1,k) );
      } // enddo
    } // enddo

    j = ny-2;
    // do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                    ( u(m,i,j-2,k) - 4.0*u(m,i,j-1,k) +
                    5.0*u(m,i,j,k) );
      } // enddo
    } // enddo
} // pragma omp parallel
  } // enddo

//---------------------------------------------------------------------
//     compute zeta-direction fluxes
//---------------------------------------------------------------------
  // do k = 1, nz-2
  for (int k = 1; k <= nz-2; ++k) {
#pragma omp task depend(out: rhs_ptr[k]) depend(in:square_ptr[k], square_ptr[k-1], square_ptr[k+1]) firstprivate(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz, k)
{
    // do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      // do i = 1, nx-2
      for (SizeT i = 1; i <= nx-2; ++i) {
        ValueT wijk = ws(i,j,k);
        ValueT wp1  = ws(i,j,k+1);
        ValueT wm1  = ws(i,j,k-1);

        rhs(1,i,j,k) = rhs(1,i,j,k) + dz1tz1 *
                  (u(1,i,j,k+1) - 2.0*u(1,i,j,k) +
                  u(1,i,j,k-1)) -
                  tz2 * (u(4,i,j,k+1) - u(4,i,j,k-1));
        rhs(2,i,j,k) = rhs(2,i,j,k) + dz2tz1 *
                  (u(2,i,j,k+1) - 2.0*u(2,i,j,k) +
                  u(2,i,j,k-1)) +
                  zzcon2 * (us(i,j,k+1) - 2.0*us(i,j,k) +
                  us(i,j,k-1)) -
                  tz2 * (u(2,i,j,k+1)*wp1 -
                  u(2,i,j,k-1)*wm1);
        rhs(3,i,j,k) = rhs(3,i,j,k) + dz3tz1 *
                  (u(3,i,j,k+1) - 2.0*u(3,i,j,k) +
                  u(3,i,j,k-1)) +
                  zzcon2 * (vs(i,j,k+1) - 2.0*vs(i,j,k) +
                  vs(i,j,k-1)) -
                  tz2 * (u(3,i,j,k+1)*wp1 -
                  u(3,i,j,k-1)*wm1);
        rhs(4,i,j,k) = rhs(4,i,j,k) + dz4tz1 *
                  (u(4,i,j,k+1) - 2.0*u(4,i,j,k) +
                  u(4,i,j,k-1)) +
                  zzcon2*con43 * (wp1 - 2.0*wijk + wm1) -
                  tz2 * (u(4,i,j,k+1)*wp1 -
                  u(4,i,j,k-1)*wm1 +
                  (u(5,i,j,k+1) - square(i,j,k+1) -
                  u(5,i,j,k-1) + square(i,j,k-1))
                  *c2);
        rhs(5,i,j,k) = rhs(5,i,j,k) + dz5tz1 *
                  (u(5,i,j,k+1) - 2.0*u(5,i,j,k) +
                  u(5,i,j,k-1)) +
                  zzcon3 * (qs(i,j,k+1) - 2.0*qs(i,j,k) +
                  qs(i,j,k-1)) +
                  zzcon4 * (wp1*wp1 - 2.0*wijk*wijk +
                  wm1*wm1) +
                  zzcon5 * (u(5,i,j,k+1)*rho_i(i,j,k+1) -
                  2.0*u(5,i,j,k)*rho_i(i,j,k) +
                  u(5,i,j,k-1)*rho_i(i,j,k-1)) -
                  tz2 * ( (c1*u(5,i,j,k+1) -
                  c2*square(i,j,k+1))*wp1 -
                  (c1*u(5,i,j,k-1) -
                  c2*square(i,j,k-1))*wm1);
      } // enddo
    } // enddo
} // pragma omp task
  } // enddo

//---------------------------------------------------------------------
//     add fourth order zeta-direction dissipation
//---------------------------------------------------------------------
  int k = 1;
#pragma omp task depend(out: rhs_ptr[k]) firstprivate(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz, k)
{
  for (SizeT j = 1; j <= ny-2; ++j) {
    // do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = rhs(m,i,j,k)- dssp *
                    ( 5.0*u(m,i,j,k) - 4.0*u(m,i,j,k+1) +
                    u(m,i,j,k+2));
      } // enddo
    } // enddo
  } // enddo
} // pragma omp task

  k = 2;
#pragma omp task depend(out: rhs_ptr[k]) firstprivate(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz, k)
{
  // do j = 1, ny-2
  for (SizeT j = 1; j <= ny-2; ++j) {
    // do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
        rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                    (-4.0*u(m,i,j,k-1) + 6.0*u(m,i,j,k) -
                      4.0*u(m,i,j,k+1) + u(m,i,j,k+2));
      } // enddo
    } // enddo
  } // enddo
} // pragma omp task

  for (int k = 3; k <= nz-4; ++k) {
#pragma omp task depend(out: rhs_ptr[k]) firstprivate(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz, k)
{
    // do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      //do i = 1,nx-2
      for (SizeT i = 1; i <= nx-2; ++i) {
        // do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                    (  u(m,i,j,k-2) - 4.0*u(m,i,j,k-1) +
                     6.0*u(m,i,j,k) - 4.0*u(m,i,j,k+1) +
                     u(m,i,j,k+2) );
        } // enddo
      } // enddo
    } // enddo
} // pragma omp task
  } // enddo

  k = nz-3;
#pragma omp task depend(out: rhs_ptr[k]) firstprivate(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz, k)
{
  //do j = 1, ny-2
  for (SizeT j = 1; j <= ny-2; ++j) {
    // do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                    ( u(m,i,j,k-2) - 4.0*u(m,i,j,k-1) +
                    6.0*u(m,i,j,k) - 4.0*u(m,i,j,k+1) );
      } // enddo
    } // enddo
  } // enddo
} // pragma omp task

  k = nz-2;
#pragma omp task depend(out: rhs_ptr[k]) firstprivate(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz, k)
{
  // do j = 1, ny-2
  for (SizeT j = 1; j <= ny-2; ++j) {
    // do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      // do m = 1, 5
      for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
                     ( u(m,i,j,k-2) - 4.0*u(m,i,j,k-1) +
                     5.0*u(m,i,j,k) );
      } // enddo
    } // enddo
  } // enddo
} // pragma omp task

  for (int k = 1; k <= nz-2; ++k) {
#pragma omp task depend(out: rhs_ptr[k]) firstprivate(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz, k)
{
    // do j = 1, ny-2
    for (SizeT j = 1; j <= ny-2; ++j) {
      // do i = 1, nx-2
      for (SizeT i = 1; i <= nx-2; ++i) {
        //do m = 1, 5
        for (int m = 1; m <= 5; ++m) {
          rhs(m,i,j,k) = rhs(m,i,j,k) * dt;
        } // enddo
      } // enddo
    } // enddo
} // pragma omp task
  } // enddo
} // omp parallel

  return;
}

#endif // HAVE_RHS_H
