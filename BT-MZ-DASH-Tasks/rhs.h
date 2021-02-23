#ifndef HAVE_RHS_H
#define HAVE_RHS_H

#include "NArray.h"

template<typename ValueT, typename SizeT>
void compute_rhs(
  ValueT *__restrict rho_i_ptr,
  ValueT *__restrict us_ptr,
  ValueT *__restrict vs_ptr,
  ValueT *__restrict ws_ptr,
  ValueT *__restrict qs_ptr,
  ValueT *__restrict square_ptr,
  ValueT *__restrict rhs_ptr,
  ValueT *__restrict forcing_ptr,
  ValueT *__restrict u_ptr,
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
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//---------------------------------------------------------------------

/*
      include 'header.h'

      integer nx, nxmax, ny, nz
      double precision rho_i  (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 us     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 vs     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 ws     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 qs     (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 square (  0:nxmax-1,0:ny-1,0:nz-1),
     $                 rhs    (5,0:nxmax-1,0:ny-1,0:nz-1),
     $                 forcing(5,0:nxmax-1,0:ny-1,0:nz-1),
     $                 u      (5,0:nxmax-1,0:ny-1,0:nz-1)

      integer i, j, k, m
      double precision rho_inv, uijk, up1, um1, vijk, vp1, vm1,
     >     wijk, wp1, wm1
*/
  // if (timeron) call timer_start(t_rhs)
//!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,i,j,k,rho_inv,
//!$OMP&  um1,up1,uijk,vm1,vp1,vijk,wm1,wp1,wijk)
//!$OMP& SHARED(dssp,c1,c2,con43,dt,nx,ny,nz,xxcon2,yycon2,zzcon2,
//!$OMP&  xxcon5,xxcon4,xxcon3,dx5tx1,dx4tx1,dx3tx1,dx2tx1,tx2,dx1tx1,
//!$OMP&  yycon5,yycon4,yycon3,dy5ty1,dy4ty1,dy3ty1,dy2ty1,ty2,dy1ty1,
//!$OMP&  zzcon5,zzcon4,zzcon3,dz5tz1,dz4tz1,dz3tz1,dz2tz1,tz2,dz1tz1)
//---------------------------------------------------------------------
//     compute the reciprocal of density and the kinetic energy,
//---------------------------------------------------------------------
//!$OMP DO SCHEDULE(STATIC)
//#pragma omp parallel
{
//#pragma omp for schedule(static) nowait
  // do k = 0, nz-1
  //for (SizeT k = 0; k <= nz-1; ++k) {
  dash::tasks::taskloop(0, nz-1+1, dash::tasks::chunk_size(1),
    [=](auto k, auto) {

    ValueT *__restrict rho_i_ptr_ = rho_i_ptr;
    ValueT *__restrict us_ptr_    = us_ptr;
    ValueT *__restrict vs_ptr_    = vs_ptr;
    ValueT *__restrict ws_ptr_    = ws_ptr;
    ValueT *__restrict qs_ptr_    = qs_ptr;
    ValueT *__restrict square_ptr_= square_ptr;
    ValueT *__restrict u_ptr_     = u_ptr;

    ArrayViewT3<ValueT, SizeT> rho_i(rho_i_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> us(us_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> vs(vs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> ws(ws_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> qs(qs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> square(square_ptr_, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> u(u_ptr_, 5, nxmax, ny, nz);

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
  },
  [&](auto k, auto, auto deps){
    deps = dash::tasks::in(u(1, 1, 1, k));
    deps = dash::tasks::out(square(1, 1, k));
  }); // enddo
//!$OMP END DO nowait

//---------------------------------------------------------------------
// copy the exact forcing term to the right hand side;  because
// this forcing term is known, we can store it on the whole zone
// including the boundary
//---------------------------------------------------------------------

//!$OMP DO SCHEDULE(STATIC)
//#pragma omp for schedule(static)
  // do k = 0, nz-1
  //for (SizeT k = 0; k <= nz-1; ++k) {
  dash::tasks::taskloop(0, nz-1+1, dash::tasks::chunk_size(1),
    [=](auto k, auto) {
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);
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
  },
  [&](auto k, auto, auto deps){
    deps = dash::tasks::out(rhs(1, 1, 1, k));
  }); // enddo
// !$OMP END DO

// !$OMP MASTER
  // if (timeron) call timer_start(t_rhsx)
//!$OMP END MASTER
//---------------------------------------------------------------------
//     compute xi-direction fluxes
//---------------------------------------------------------------------
//!$OMP DO SCHEDULE(STATIC)
  // do k = 1, nz-2
//#pragma omp for schedule(static) nowait
  //for (SizeT k = 1; k <= nz-2; ++k) {
  dash::tasks::taskloop(1, nz-2+1, dash::tasks::chunk_size(1),
    [=](auto k, auto) {
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);

    ValueT *__restrict rho_i_ptr_ = rho_i_ptr;
    ValueT *__restrict us_ptr_    = us_ptr;
    ValueT *__restrict vs_ptr_    = vs_ptr;
    ValueT *__restrict ws_ptr_    = ws_ptr;
    ValueT *__restrict qs_ptr_    = qs_ptr;
    ValueT *__restrict square_ptr_= square_ptr;
    ValueT *__restrict u_ptr_     = u_ptr;

    ArrayViewT3<ValueT, SizeT> rho_i(rho_i_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> us(us_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> vs(vs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> ws(ws_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> qs(qs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> square(square_ptr_, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> u(u_ptr_, 5, nxmax, ny, nz);
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
  },
  [&](auto k, auto, auto deps){
    deps = dash::tasks::out(rhs(1, 1, 1, k));
    deps = dash::tasks::in(u(1, 1, 1, k));
    deps = dash::tasks::in(square(1, 1, k));
  }); // enddo
//!$OMP END DO nowait
//!$OMP MASTER
      //if (timeron) call timer_stop(t_rhsx)

      // if (timeron) call timer_start(t_rhsy)
//!$OMP END MASTER
//---------------------------------------------------------------------
//     compute eta-direction fluxes
//---------------------------------------------------------------------
//!$OMP DO SCHEDULE(STATIC)
//#pragma omp for schedule(static) nowait
  // do k = 1, nz-2
  //for (SizeT k = 1; k <= nz-2; ++k) {
  dash::tasks::taskloop(1, nz-2+1, dash::tasks::chunk_size(1),
    [=](auto k, auto) {
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);

    ValueT *__restrict rho_i_ptr_ = rho_i_ptr;
    ValueT *__restrict us_ptr_    = us_ptr;
    ValueT *__restrict vs_ptr_    = vs_ptr;
    ValueT *__restrict ws_ptr_    = ws_ptr;
    ValueT *__restrict qs_ptr_    = qs_ptr;
    ValueT *__restrict square_ptr_= square_ptr;
    ValueT *__restrict u_ptr_     = u_ptr;

    ArrayViewT3<ValueT, SizeT> rho_i(rho_i_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> us(us_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> vs(vs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> ws(ws_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> qs(qs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> square(square_ptr_, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> u(u_ptr_, 5, nxmax, ny, nz);

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
  },
  [&](auto k, auto, auto deps){
    deps = dash::tasks::out(rhs(1, 1, 1, k));
    deps = dash::tasks::in(u(1, 1, 1, k));
    deps = dash::tasks::in(square(1, 1, k));
  }); // enddo
//!$OMP END DO nowait
//!$OMP MASTER
  //if (timeron) call timer_stop(t_rhsy)

  //if (timeron) call timer_start(t_rhsz)
//!$OMP END MASTER
//---------------------------------------------------------------------
//     compute zeta-direction fluxes
//---------------------------------------------------------------------
//!$OMP DO SCHEDULE(STATIC)
//#pragma omp for schedule(static)
  // do k = 1, nz-2
  // for (SizeT k = 1; k <= nz-2; ++k) {
  dash::tasks::taskloop(1, nz-2+1, dash::tasks::chunk_size(1),
    [=](auto k, auto) {
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);

    ValueT *__restrict rho_i_ptr_ = rho_i_ptr;
    ValueT *__restrict us_ptr_    = us_ptr;
    ValueT *__restrict vs_ptr_    = vs_ptr;
    ValueT *__restrict ws_ptr_    = ws_ptr;
    ValueT *__restrict qs_ptr_    = qs_ptr;
    ValueT *__restrict square_ptr_= square_ptr;
    ValueT *__restrict u_ptr_     = u_ptr;

    ArrayViewT3<ValueT, SizeT> rho_i(rho_i_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> us(us_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> vs(vs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> ws(ws_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> qs(qs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> square(square_ptr_, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> u(u_ptr_, 5, nxmax, ny, nz);

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
  },
  [&](auto k, auto, auto deps){
    deps = dash::tasks::out(rhs(1, 1, 1, k));
    deps = dash::tasks::in(u(1, 1, 1, k));
    deps = dash::tasks::in(u(1, 1, 1, k-1));
    deps = dash::tasks::in(u(1, 1, 1, k+1));
    deps = dash::tasks::in(square(1, 1, k));
    deps = dash::tasks::in(square(1, 1, k-1));
    deps = dash::tasks::in(square(1, 1, k+1));
  }); // enddo
//!$OMP END DO

//---------------------------------------------------------------------
//     add fourth order zeta-direction dissipation
//---------------------------------------------------------------------
  SizeT k = 1;
// !$OMP DO SCHEDULE(STATIC)
  // do j = 1, ny-2
//#pragma omp for schedule (static) nowait
  dash::tasks::async([=]() {
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);
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
  },
  dash::tasks::out(rhs(1, 1, 1, k)),
  dash::tasks::in(u(1, 1, 1, k)));
// !$OMP END DO nowait

  k = 2;
//!$OMP DO SCHEDULE(STATIC)
//#pragma omp for schedule(static) nowait
  dash::tasks::async([=]() {
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);
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
  },
  dash::tasks::out(rhs(1, 1, 1, k)),
  dash::tasks::in(u(1, 1, 1, k-1)),
  dash::tasks::in(u(1, 1, 1, k)),
  dash::tasks::in(u(1, 1, 1, k+1)),
  dash::tasks::in(u(1, 1, 1, k+2)));

// !$OMP END DO nowait

//!$OMP DO SCHEDULE(STATIC)
  // do k = 3, nz-4
//#pragma omp for schedule(static) nowait
  //for (SizeT k = 3; k <= nz-4; ++k) {
  dash::tasks::taskloop(3, nz-4+1, dash::tasks::chunk_size(1),
    [=](auto k, auto) {
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);

    ValueT *__restrict rho_i_ptr_ = rho_i_ptr;
    ValueT *__restrict us_ptr_    = us_ptr;
    ValueT *__restrict vs_ptr_    = vs_ptr;
    ValueT *__restrict ws_ptr_    = ws_ptr;
    ValueT *__restrict qs_ptr_    = qs_ptr;
    ValueT *__restrict square_ptr_= square_ptr;
    ValueT *__restrict u_ptr_     = u_ptr;

    ArrayViewT3<ValueT, SizeT> rho_i(rho_i_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> us(us_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> vs(vs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> ws(ws_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> qs(qs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> square(square_ptr_, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> u(u_ptr_, 5, nxmax, ny, nz);

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
  },
  [&](auto k, auto, auto deps){
    deps = dash::tasks::out(rhs(1, 1, 1, k));
    deps = dash::tasks::in(u(1, 1, 1, k-2));
    deps = dash::tasks::in(u(1, 1, 1, k-1));
    deps = dash::tasks::in(u(1, 1, 1, k));
    deps = dash::tasks::in(u(1, 1, 1, k+1));
    deps = dash::tasks::in(u(1, 1, 1, k+2));
  }); // enddo
//!$OMP END DO nowait

  k = nz-3;
//!$OMP DO SCHEDULE(STATIC)
//#pragma omp for schedule(static) nowait
  dash::tasks::async([=]() {
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);
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
  },
  dash::tasks::out(rhs(1, 1, 1, k)),
  dash::tasks::in(u(1, 1, 1, k-2)),
  dash::tasks::in(u(1, 1, 1, k-1)),
  dash::tasks::in(u(1, 1, 1, k)),
  dash::tasks::in(u(1, 1, 1, k+1)));
//!$OMP END DO nowait

  k = nz-2;
//!$OMP DO SCHEDULE(STATIC)
//#pragma omp for schedule(static)
  dash::tasks::async([=](){
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);
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
  },
  dash::tasks::out(rhs(1, 1, 1, k)),
  dash::tasks::in(u(1, 1, 1, k-2)),
  dash::tasks::in(u(1, 1, 1, k-1)),
  dash::tasks::in(u(1, 1, 1, k)));

//!$OMP END DO
//!$OMP MASTER
  // if (timeron) call timer_stop(t_rhsz)
//!$OMP END MASTER

//!$OMP DO SCHEDULE(STATIC) nowait
  // do k = 1, nz-2
  //for (SizeT k = 1; k <= nz-2; ++k) {
  dash::tasks::taskloop(1, nz-2+1, dash::tasks::chunk_size(1),
    [=](auto k, auto) {
    ValueT *__restrict rhs_ptr_ = rhs_ptr;
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);
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
  },
  [&](auto k, auto, auto deps){
    deps = dash::tasks::out(rhs(1, 1, 1, k));
  }); // enddo
//!$OMP END DO nowait
} // omp parallel
//!$OMP END PARALLEL
  // if (timeron) call timer_stop(t_rhs)

  return;

}



#endif // HAVE_RHS_H
