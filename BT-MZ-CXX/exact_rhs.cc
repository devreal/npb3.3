
#include <iostream>

#include "NArray.h"
#include "header.h"

#include "exact_solution.h"

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
  ArrayViewT4<double, size_type> forcing(forcing_ptr, 5, nxmax, ny, nz);
#pragma omp parallel
{
  // thread_local
  StaticArrayT1<double, size_type, problem_size+1, 0> cuf{}/*(0:problem_size)*/;
  StaticArrayT1<double, size_type, problem_size+1, 0> q{}/*(0:problem_size)*/;
  // thread_local
  StaticArrayT2<double, size_type, problem_size+1, 5, 0, 1> ue{}/*(0:problem_size,5)*/;
  StaticArrayT2<double, size_type, problem_size+1, 5, 0, 1> buf{}/*(0:problem_size,5)*/;
  StaticArrayT1<double, int, 5, 1> dtemp;
/*
c---------------------------------------------------------------------
c     initialize
c---------------------------------------------------------------------
*/
#pragma omp for
  for (int k= 0; k <= nz-1; ++k) {
    for (int j = 0; j <= ny-1; ++j) {
      for (int i = 0; i <= nx-1; ++i) {
        for (int m = 1; m <= 5; ++m) {
          forcing(m,i,j,k) = 0.0;
        } // enddo
      } // enddo
    } // enddo
  } // enddo

/*
c---------------------------------------------------------------------
c     xi-direction flux differences
c---------------------------------------------------------------------
*/
#pragma omp for
  for (int k = 1; k <= nz-2; ++k) {
    double zeta = k * dnzm1;
    for (int j = 1; j <= ny-2; ++j) {
      double eta = j * dnym1;

      for (int i=0; i <= nx-1; ++i) {
          double xi = i * dnxm1;

          exact_solution(xi, eta, zeta, dtemp.begin());
          // doloop(1, 5, [&](auto m) {
          for (int m = 1; m <= 5; ++m) {
            ue(i,m) = dtemp(m);
          } // ); // enddo

          double dtpp = 1.0 / dtemp(1);

          // doloop(2, 5, [&](auto m) {
          for (int m = 2; m <= 5; ++m) {
            buf(i,m) = dtpp * dtemp(m);
          } // ); // enddo

          cuf(i)   = buf(i,2) * buf(i,2);
          buf(i,1) = cuf(i) + buf(i,3) * buf(i,3) +
                        buf(i,4) * buf(i,4);
          q(i) = 0.50*(buf(i,2)*ue(i,2) + buf(i,3)*ue(i,3) +
                        buf(i,4)*ue(i,4));

      } // enddo

      // doloop(1, nx-2, [&](auto i){
      for (size_type i = 1; i <= nx-2; ++i) {
          int im1 = i-1;
          int ip1 = i+1;

          forcing(1,i,j,k) = forcing(1,i,j,k) -
                 tx2*( ue(ip1,2)-ue(im1,2) )+
                 dx1tx1*(ue(ip1,1)-2.0*ue(i,1)+ue(im1,1));

          forcing(2,i,j,k) = forcing(2,i,j,k) - tx2 * (
                 (ue(ip1,2)*buf(ip1,2)+c2*(ue(ip1,5)-q(ip1)))-
                 (ue(im1,2)*buf(im1,2)+c2*(ue(im1,5)-q(im1))))+
                 xxcon1*(buf(ip1,2)-2.0*buf(i,2)+buf(im1,2))+
                 dx2tx1*( ue(ip1,2)-2.0* ue(i,2)+ue(im1,2));

          forcing(3,i,j,k) = forcing(3,i,j,k) - tx2 * (
                 ue(ip1,3)*buf(ip1,2)-ue(im1,3)*buf(im1,2))+
                 xxcon2*(buf(ip1,3)-2.0*buf(i,3)+buf(im1,3))+
                 dx3tx1*( ue(ip1,3)-2.0*ue(i,3) +ue(im1,3));

          forcing(4,i,j,k) = forcing(4,i,j,k) - tx2*(
                 ue(ip1,4)*buf(ip1,2)-ue(im1,4)*buf(im1,2))+
                 xxcon2*(buf(ip1,4)-2.0*buf(i,4)+buf(im1,4))+
                 dx4tx1*( ue(ip1,4)-2.0* ue(i,4)+ ue(im1,4));

          forcing(5,i,j,k) = forcing(5,i,j,k) - tx2*(
                 buf(ip1,2)*(c1*ue(ip1,5)-c2*q(ip1))-
                 buf(im1,2)*(c1*ue(im1,5)-c2*q(im1)))+
                 0.5*xxcon3*(buf(ip1,1)-2.0*buf(i,1)+
                 buf(im1,1))+
                 xxcon4*(cuf(ip1)-2.0*cuf(i)+cuf(im1))+
                 xxcon5*(buf(ip1,5)-2.0*buf(i,5)+buf(im1,5))+
                 dx5tx1*( ue(ip1,5)-2.0* ue(i,5)+ ue(im1,5));
      } // ); // enddo

/*
c---------------------------------------------------------------------
c     Fourth-order dissipation
c---------------------------------------------------------------------
*/
      // doloop(1, 5, [&](auto m) {
      for (int m = 1; m <= 5; ++m) {
          int i = 1;
          forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                    (5.0*ue(i,m) - 4.0*ue(i+1,m) +ue(i+2,m));
          i = 2;
          forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                    (-4.0*ue(i-1,m) + 6.0*ue(i,m) -
                    4.0*ue(i+1,m) +       ue(i+2,m));
      } // ); // enddo

      // doloop(1, 5, [&](auto m) {
      for (int m = 1; m <= 5; ++m) {
          // doloop(3, nx-4, [&](auto i){
          for (size_type i = 3; i <= nx-4; ++i) {
            forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
                    (ue(i-2,m) - 4.0*ue(i-1,m) +
                    6.0*ue(i,m) - 4.0*ue(i+1,m) + ue(i+2,m));
          } // ); // enddo
      } // ); // enddo

      //doloop(1, 5, [&](auto m) {
      for (int m = 1; m <= 5; ++m) {
          auto i = nx-3;
          forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                    (ue(i-2,m) - 4.0*ue(i-1,m) +
                    6.0*ue(i,m) - 4.0*ue(i+1,m));
          i = nx-2;
          forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                    (ue(i-2,m) - 4.0*ue(i-1,m) + 5.0*ue(i,m));
      } // ); // enddo

    } // enddo
  } //enddo

/*
c---------------------------------------------------------------------
c     eta-direction flux differences
c---------------------------------------------------------------------
*/
#pragma omp for
  for (int k = 1; k <= nz-2; ++k) {
    double zeta = k * dnzm1;
    // doloop(1, nx-2, [&](auto i){
    for (size_type i = 1; i <= nx-2; ++i) {
      double xi = i * dnxm1;

      // doloop(0, ny-1, [&](auto j){
      for (size_type j = 0; j <= ny-1; ++j) {
        double eta = j * dnym1;

        exact_solution(xi, eta, zeta, dtemp.begin());
        // doloop(1, 5, [&](auto m){
        for (int m = 1; m <= 5; ++m) {
          ue(j,m) = dtemp(m);
        } // ); // enddo

        double dtpp = 1.0/dtemp(1);

        // doloop(2, 5, [&](auto m){
        for (int m = 2; m <= 5; ++m) {
          buf(j,m) = dtpp * dtemp(m);
        } // ); // enddo

        cuf(j)   = buf(j,3) * buf(j,3);
        buf(j,1) = cuf(j) + buf(j,2) * buf(j,2) +
              buf(j,4) * buf(j,4);
        q(j) = 0.5*(buf(j,2)*ue(j,2) + buf(j,3)*ue(j,3) +
              buf(j,4)*ue(j,4));
      } // ); //enddo

      //doloop(1, ny-2, [&](auto j){
      for (size_type j = 1; j <= ny-2; ++j) {
        auto jm1 = j-1;
        auto jp1 = j+1;

        forcing(1,i,j,k) = forcing(1,i,j,k) -
              ty2*( ue(jp1,3)-ue(jm1,3) )+
              dy1ty1*(ue(jp1,1)-2.0*ue(j,1)+ue(jm1,1));

        forcing(2,i,j,k) = forcing(2,i,j,k) - ty2*(
              ue(jp1,2)*buf(jp1,3)-ue(jm1,2)*buf(jm1,3))+
              yycon2*(buf(jp1,2)-2.0*buf(j,2)+buf(jm1,2))+
              dy2ty1*( ue(jp1,2)-2.0* ue(j,2)+ ue(jm1,2));

        forcing(3,i,j,k) = forcing(3,i,j,k) - ty2*(
              (ue(jp1,3)*buf(jp1,3)+c2*(ue(jp1,5)-q(jp1)))-
              (ue(jm1,3)*buf(jm1,3)+c2*(ue(jm1,5)-q(jm1))))+
              yycon1*(buf(jp1,3)-2.0*buf(j,3)+buf(jm1,3))+
              dy3ty1*( ue(jp1,3)-2.0*ue(j,3) +ue(jm1,3));

        forcing(4,i,j,k) = forcing(4,i,j,k) - ty2*(
              ue(jp1,4)*buf(jp1,3)-ue(jm1,4)*buf(jm1,3))+
              yycon2*(buf(jp1,4)-2.0*buf(j,4)+buf(jm1,4))+
              dy4ty1*( ue(jp1,4)-2.0*ue(j,4)+ ue(jm1,4));

        forcing(5,i,j,k) = forcing(5,i,j,k) - ty2*(
              buf(jp1,3)*(c1*ue(jp1,5)-c2*q(jp1))-
              buf(jm1,3)*(c1*ue(jm1,5)-c2*q(jm1)))+
              0.5*yycon3*(buf(jp1,1)-2.0*buf(j,1)+
              buf(jm1,1))+
              yycon4*(cuf(jp1)-2.0*cuf(j)+cuf(jm1))+
              yycon5*(buf(jp1,5)-2.0*buf(j,5)+buf(jm1,5))+
              dy5ty1*(ue(jp1,5)-2.0*ue(j,5)+ue(jm1,5));
      } // ); // enddo

/*
c---------------------------------------------------------------------
c     Fourth-order dissipation
c---------------------------------------------------------------------
*/
      // doloop(1, 5, [&](auto m){
      for (int m = 1; m <= 5; ++m) {
        auto j = 1;
        forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                  (5.0*ue(j,m) - 4.0*ue(j+1,m) +ue(j+2,m));
        j = 2;
        forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                  (-4.0*ue(j-1,m) + 6.0*ue(j,m) -
                  4.0*ue(j+1,m) +       ue(j+2,m));
      } // ); // enddo

      // doloop(1, 5, [&](auto m){
      for (int m = 1; m <= 5; ++m) {
        // doloop(3, ny-4, [&](auto j){
        for (size_type j = 3; j <= ny-4; ++j) {
          forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
                  (ue(j-2,m) - 4.0*ue(j-1,m) +
                  6.0*ue(j,m) - 4.0*ue(j+1,m) + ue(j+2,m));
        } // ); //enddo
      } // ); // enddo

      // doloop(1, 5, [&](auto m){
      for (int m = 1; m <= 5; ++m) {
        auto j = ny-3;
        forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                  (ue(j-2,m) - 4.0*ue(j-1,m) +
                  6.0*ue(j,m) - 4.0*ue(j+1,m));
        j = ny-2;
        forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                  (ue(j-2,m) - 4.0*ue(j-1,m) + 5.0*ue(j,m));

      } // ); // enddo

    } // ); // enddo
  } // enddo

/*
c---------------------------------------------------------------------
c     zeta-direction flux differences
c---------------------------------------------------------------------
*/
#pragma omp for
  for (int j=1; j <= ny-2; ++j) {
    double eta = j * dnym1;
    //doloop(1, nx-2, [&](auto i){
    for (size_type i = 1; i <= nx-2; ++i) {
      double xi = i * dnxm1;

      // doloop(0, nz-1, [&](auto k) {
      for (size_type k = 0; k <= nz-1; ++k) {
        double zeta = k * dnzm1;

        exact_solution(xi, eta, zeta, dtemp.begin());
        // doloop(1, 5, [&](auto m){
        for (int m = 1; m <= 5; ++m) {
          ue(k,m) = dtemp(m);
        } // ); // enddo

        auto dtpp = 1.0/dtemp(1);

        // doloop(2, 5, [&](auto m){
        for (int m = 2; m <= 5; ++m) {
          buf(k,m) = dtpp * dtemp(m);
        } // ); // enddo

        cuf(k)   = buf(k,4) * buf(k,4);
        buf(k,1) = cuf(k) + buf(k,2) * buf(k,2) +
                buf(k,3) * buf(k,3);
        q(k) = 0.5*(buf(k,2)*ue(k,2) + buf(k,3)*ue(k,3) +
                buf(k,4)*ue(k,4));
      } // ); // enddo

      // doloop(1, nz-2, [&](auto k){
      for (size_type k = 1; k <= nz-2; ++k) {
        auto km1 = k-1;
        auto kp1 = k+1;

        forcing(1,i,j,k) = forcing(1,i,j,k) -
              tz2*( ue(kp1,4)-ue(km1,4) )+
              dz1tz1*(ue(kp1,1)-2.0*ue(k,1)+ue(km1,1));

        forcing(2,i,j,k) = forcing(2,i,j,k) - tz2 * (
              ue(kp1,2)*buf(kp1,4)-ue(km1,2)*buf(km1,4))+
              zzcon2*(buf(kp1,2)-2.0*buf(k,2)+buf(km1,2))+
              dz2tz1*( ue(kp1,2)-2.0* ue(k,2)+ ue(km1,2));

        forcing(3,i,j,k) = forcing(3,i,j,k) - tz2 * (
              ue(kp1,3)*buf(kp1,4)-ue(km1,3)*buf(km1,4))+
              zzcon2*(buf(kp1,3)-2.0*buf(k,3)+buf(km1,3))+
              dz3tz1*(ue(kp1,3)-2.0*ue(k,3)+ue(km1,3));

        forcing(4,i,j,k) = forcing(4,i,j,k) - tz2 * (
              (ue(kp1,4)*buf(kp1,4)+c2*(ue(kp1,5)-q(kp1)))-
              (ue(km1,4)*buf(km1,4)+c2*(ue(km1,5)-q(km1))))+
              zzcon1*(buf(kp1,4)-2.0*buf(k,4)+buf(km1,4))+
              dz4tz1*( ue(kp1,4)-2.0*ue(k,4) +ue(km1,4));

        forcing(5,i,j,k) = forcing(5,i,j,k) - tz2 * (
              buf(kp1,4)*(c1*ue(kp1,5)-c2*q(kp1))-
              buf(km1,4)*(c1*ue(km1,5)-c2*q(km1)))+
              0.5*zzcon3*(buf(kp1,1)-2.0*buf(k,1)
              +buf(km1,1))+
              zzcon4*(cuf(kp1)-2.0*cuf(k)+cuf(km1))+
              zzcon5*(buf(kp1,5)-2.0*buf(k,5)+buf(km1,5))+
              dz5tz1*( ue(kp1,5)-2.0*ue(k,5)+ ue(km1,5));
      } // ); // enddo

/*
c---------------------------------------------------------------------
c     Fourth-order dissipation
c---------------------------------------------------------------------
*/
      // doloop(1, 5, [&](auto m){
      for (int m = 1; m <= 5; ++m) {
        auto k = 1;
        forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                  (5.0*ue(k,m) - 4.0*ue(k+1,m) +ue(k+2,m));
        k = 2;
        forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                  (-4.0*ue(k-1,m) + 6.0*ue(k,m) -
                    4.0*ue(k+1,m) +     ue(k+2,m));
      } // ); // enddo

      // doloop(1, 5, [&](auto m){
      for (int m = 1; m <= 5; ++m) {
        // doloop(3, nz-4, [&](auto k){
        for (size_type k = 3; k <= nz-4; ++k) {
          forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
                  (ue(k-2,m) - 4.0*ue(k-1,m) +
                  6.0*ue(k,m) - 4.0*ue(k+1,m) + ue(k+2,m));
        } // ); // enddo
      } // ); // enddo

      // doloop(1, 5, [&](auto m){
      for (int m = 1; m <= 5; ++m) {
        auto k = nz-3;
        forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                  (ue(k-2,m) - 4.0*ue(k-1,m) +
                  6.0*ue(k,m) - 4.0*ue(k+1,m));
        k = nz-2;
        forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
                  (ue(k-2,m) - 4.0*ue(k-1,m) + 5.0*ue(k,m));
      } // ); // enddo

    } // ); // enddo
  } // enddo

/*
c---------------------------------------------------------------------
c     now change the sign of the forcing function,
c---------------------------------------------------------------------
*/
#pragma omp for
  for (int k = 1; k <= nz-2; ++k) {
    // doloop(1, ny-2, [&](auto j) {
    for (size_type j = 1; j <= ny-2; ++j) {
      // doloop(1, nx-2, [&](auto i) {
      for (size_type i = 1; i <= nx-2; ++i) {
        //doloop(1, 5, [&](auto m) {
        for (int m = 1; m <= 5; ++m) {
          forcing(m,i,j,k) = -1.0 * forcing(m,i,j,k);
        } // ); // enddo
      } // ); // enddo
    } // ); // enddo
  } // enddo
}

  return;
}
