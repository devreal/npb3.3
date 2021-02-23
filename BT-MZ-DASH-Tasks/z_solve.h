
#include "NArray.h"

#define TASKLOOP_FACTOR 2

template<typename SizeT, typename ValueT>
void z_solve(
  ValueT *__restrict rho_i_ptr,
  ValueT *__restrict qs_ptr,
  ValueT *__restrict square_ptr,
  ValueT *__restrict u_ptr,
  ValueT *__restrict rhs_ptr,
  SizeT nx,
  SizeT nxmax,
  SizeT ny,
  SizeT nz)
{

//#pragma omp parallel
{
  SizeT ksize = nz-1;

  if (timeron) timer_start(t_zsolve);

//---------------------------------------------------------------------
//     Performs line solves in Z direction by first factoring
//     the block-tridiagonal matrix into an upper triangular matrix,
//     and then performing back substitution to solve for the unknow
//     vectors of each line.
//
//     Make sure we treat elements zero to cell_size in the direction
//     of the sweep.
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//---------------------------------------------------------------------

      //if (timeron) call timer_start(t_zsolve)

//---------------------------------------------------------------------
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     This function computes the left hand side for the three z-factors
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     Compute the indices for storing the block-diagonal matrix;
//     determine c (labeled f) and s jacobians
//---------------------------------------------------------------------
  // TODO: no dependencies here, add them if needed!
//#pragma omp for
  // do j = 1, ny-2
  //for (SizeT j = 1; j <= ny-2; ++j) {
  dash::tasks::taskloop(1, ny-2+1, dash::tasks::num_chunks(dash::tasks::numthreads()*TASKLOOP_FACTOR),
    [=](auto j_begin, auto j_end) {
    // thread-local temporary variables, allocated on the stack
    StaticArrayT3<double, size_type, 5, 5,    problem_size+1, 1, 1, 0> fjac/*(5, 5, problem_size+1)*/;
    StaticArrayT3<double, size_type, 5, 5,    problem_size+1, 1, 1, 0> njac/*(5, 5, problem_size+1)*/;
    StaticArrayT4<double, size_type, 5, 5, 3, problem_size+1, 1, 1, 1, 0> lhs /*(5, 5, 3, problem_size+1)*/;
    StaticArrayT2<double, size_type, 5, problem_size+1, 1, 0> rtmp/*(5, problem_size+1)*/;

    // Hack necessary because compilers fails to properly transfer __restrict into lambda
    ValueT *__restrict rho_i_ptr_  = rho_i_ptr;
    ValueT *__restrict qs_ptr_     = qs_ptr;
    ValueT *__restrict square_ptr_ = square_ptr;
    ValueT *__restrict u_ptr_      = u_ptr;
    ValueT *__restrict rhs_ptr_    = rhs_ptr;
    ArrayViewT3<ValueT, SizeT> rho_i(rho_i_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> qs(qs_ptr_, nxmax, ny, nz);
    ArrayViewT3<ValueT, SizeT> square(square_ptr_, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> u(u_ptr_, 5, nxmax, ny, nz);
    ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr_, 5, nxmax, ny, nz);


    for (SizeT j = j_begin; j < j_end; ++j) {
    // do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      // do k = 0, ksize
      for (SizeT k = 0; k <= ksize; ++k) {

        ValueT tmp1, tmp2, tmp3;

        tmp1 = 1.0 / u(1,i,j,k);
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;

        fjac(1,1,k) = 0.0;
        fjac(1,2,k) = 0.0;
        fjac(1,3,k) = 0.0;
        fjac(1,4,k) = 1.0;
        fjac(1,5,k) = 0.0;

        fjac(2,1,k) = - ( u(2,i,j,k)*u(4,i,j,k) )
              * tmp2;
        fjac(2,2,k) = u(4,i,j,k) * tmp1;
        fjac(2,3,k) = 0.0;
        fjac(2,4,k) = u(2,i,j,k) * tmp1;
        fjac(2,5,k) = 0.0;

        fjac(3,1,k) = - ( u(3,i,j,k)*u(4,i,j,k) )
              * tmp2;
        fjac(3,2,k) = 0.0;
        fjac(3,3,k) = u(4,i,j,k) * tmp1;
        fjac(3,4,k) = u(3,i,j,k) * tmp1;
        fjac(3,5,k) = 0.0;

        fjac(4,1,k) = - (u(4,i,j,k)*u(4,i,j,k) * tmp2 )
              + c2 * qs(i,j,k);
        fjac(4,2,k) = - c2 *  u(2,i,j,k) * tmp1;
        fjac(4,3,k) = - c2 *  u(3,i,j,k) * tmp1;
        fjac(4,4,k) = ( 2.0 - c2 )
              *  u(4,i,j,k) * tmp1;
        fjac(4,5,k) = c2;

        fjac(5,1,k) = ( c2 * 2.00 * square(i,j,k)
              - c1 * u(5,i,j,k) )
              * u(4,i,j,k) * tmp2;
        fjac(5,2,k) = - c2 * ( u(2,i,j,k)*u(4,i,j,k) )
              * tmp2;
        fjac(5,3,k) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
              * tmp2;
        fjac(5,4,k) = c1 * ( u(5,i,j,k) * tmp1 )
              - c2
              * ( qs(i,j,k)
              + u(4,i,j,k)*u(4,i,j,k) * tmp2 );
        fjac(5,5,k) = c1 * u(4,i,j,k) * tmp1;

        njac(1,1,k) = 0.0;
        njac(1,2,k) = 0.0;
        njac(1,3,k) = 0.0;
        njac(1,4,k) = 0.0;
        njac(1,5,k) = 0.0;

        njac(2,1,k) = - c3c4 * tmp2 * u(2,i,j,k);
        njac(2,2,k) =   c3c4 * tmp1;
        njac(2,3,k) =   0.0;
        njac(2,4,k) =   0.0;
        njac(2,5,k) =   0.0;

        njac(3,1,k) = - c3c4 * tmp2 * u(3,i,j,k);
        njac(3,2,k) =   0.0;
        njac(3,3,k) =   c3c4 * tmp1;
        njac(3,4,k) =   0.0;
        njac(3,5,k) =   0.0;

        njac(4,1,k) = - con43 * c3c4 * tmp2 * u(4,i,j,k);
        njac(4,2,k) =   0.0;
        njac(4,3,k) =   0.0;
        njac(4,4,k) =   con43 * c3 * c4 * tmp1;
        njac(4,5,k) =   0.0;

        njac(5,1,k) = - (  c3c4
              - c1345 ) * tmp3 * (u(2,i,j,k)*u(2,i,j,k))
              - ( c3c4 - c1345 ) * tmp3 * (u(3,i,j,k)*u(3,i,j,k))
              - ( con43 * c3c4
              - c1345 ) * tmp3 * (u(4,i,j,k)*u(4,i,j,k))
              - c1345 * tmp2 * u(5,i,j,k);

        njac(5,2,k) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k);
        njac(5,3,k) = (  c3c4 - c1345 ) * tmp2 * u(3,i,j,k);
        njac(5,4,k) = ( con43 * c3c4
              - c1345 ) * tmp2 * u(4,i,j,k);
        njac(5,5,k) = ( c1345 )* tmp1;


      } // enddo

//---------------------------------------------------------------------
//     now jacobians set, so form left hand side in z direction
//---------------------------------------------------------------------
      lhsinit(lhs.begin(), ksize);
      //do k = 1, ksize-1
      for (SizeT k = 1; k <= ksize-1; ++k) {

        ValueT tmp1 = dt * tz1;
        ValueT tmp2 = dt * tz2;

        lhs(1,1,aa,k) = - tmp2 * fjac(1,1,k-1)
              - tmp1 * njac(1,1,k-1)
              - tmp1 * dz1;
        lhs(1,2,aa,k) = - tmp2 * fjac(1,2,k-1)
              - tmp1 * njac(1,2,k-1);
        lhs(1,3,aa,k) = - tmp2 * fjac(1,3,k-1)
              - tmp1 * njac(1,3,k-1);
        lhs(1,4,aa,k) = - tmp2 * fjac(1,4,k-1)
              - tmp1 * njac(1,4,k-1);
        lhs(1,5,aa,k) = - tmp2 * fjac(1,5,k-1)
              - tmp1 * njac(1,5,k-1);

        lhs(2,1,aa,k) = - tmp2 * fjac(2,1,k-1)
              - tmp1 * njac(2,1,k-1);
        lhs(2,2,aa,k) = - tmp2 * fjac(2,2,k-1)
              - tmp1 * njac(2,2,k-1)
              - tmp1 * dz2;
        lhs(2,3,aa,k) = - tmp2 * fjac(2,3,k-1)
              - tmp1 * njac(2,3,k-1);
        lhs(2,4,aa,k) = - tmp2 * fjac(2,4,k-1)
              - tmp1 * njac(2,4,k-1);
        lhs(2,5,aa,k) = - tmp2 * fjac(2,5,k-1)
              - tmp1 * njac(2,5,k-1);

        lhs(3,1,aa,k) = - tmp2 * fjac(3,1,k-1)
              - tmp1 * njac(3,1,k-1);
        lhs(3,2,aa,k) = - tmp2 * fjac(3,2,k-1)
              - tmp1 * njac(3,2,k-1);
        lhs(3,3,aa,k) = - tmp2 * fjac(3,3,k-1)
              - tmp1 * njac(3,3,k-1)
              - tmp1 * dz3;
        lhs(3,4,aa,k) = - tmp2 * fjac(3,4,k-1)
              - tmp1 * njac(3,4,k-1);
        lhs(3,5,aa,k) = - tmp2 * fjac(3,5,k-1)
              - tmp1 * njac(3,5,k-1);

        lhs(4,1,aa,k) = - tmp2 * fjac(4,1,k-1)
              - tmp1 * njac(4,1,k-1);
        lhs(4,2,aa,k) = - tmp2 * fjac(4,2,k-1)
              - tmp1 * njac(4,2,k-1);
        lhs(4,3,aa,k) = - tmp2 * fjac(4,3,k-1)
              - tmp1 * njac(4,3,k-1);
        lhs(4,4,aa,k) = - tmp2 * fjac(4,4,k-1)
              - tmp1 * njac(4,4,k-1)
              - tmp1 * dz4;
        lhs(4,5,aa,k) = - tmp2 * fjac(4,5,k-1)
              - tmp1 * njac(4,5,k-1);

        lhs(5,1,aa,k) = - tmp2 * fjac(5,1,k-1)
              - tmp1 * njac(5,1,k-1);
        lhs(5,2,aa,k) = - tmp2 * fjac(5,2,k-1)
              - tmp1 * njac(5,2,k-1);
        lhs(5,3,aa,k) = - tmp2 * fjac(5,3,k-1)
              - tmp1 * njac(5,3,k-1);
        lhs(5,4,aa,k) = - tmp2 * fjac(5,4,k-1)
              - tmp1 * njac(5,4,k-1);
        lhs(5,5,aa,k) = - tmp2 * fjac(5,5,k-1)
              - tmp1 * njac(5,5,k-1)
              - tmp1 * dz5;

        lhs(1,1,bb,k) = 1.0
              + tmp1 * 2.0 * njac(1,1,k)
              + tmp1 * 2.0 * dz1;
        lhs(1,2,bb,k) = tmp1 * 2.0 * njac(1,2,k);
        lhs(1,3,bb,k) = tmp1 * 2.0 * njac(1,3,k);
        lhs(1,4,bb,k) = tmp1 * 2.0 * njac(1,4,k);
        lhs(1,5,bb,k) = tmp1 * 2.0 * njac(1,5,k);

        lhs(2,1,bb,k) = tmp1 * 2.0 * njac(2,1,k);
        lhs(2,2,bb,k) = 1.0
              + tmp1 * 2.0 * njac(2,2,k)
              + tmp1 * 2.0 * dz2;
        lhs(2,3,bb,k) = tmp1 * 2.0 * njac(2,3,k);
        lhs(2,4,bb,k) = tmp1 * 2.0 * njac(2,4,k);
        lhs(2,5,bb,k) = tmp1 * 2.0 * njac(2,5,k);

        lhs(3,1,bb,k) = tmp1 * 2.0 * njac(3,1,k);
        lhs(3,2,bb,k) = tmp1 * 2.0 * njac(3,2,k);
        lhs(3,3,bb,k) = 1.0
              + tmp1 * 2.0 * njac(3,3,k)
              + tmp1 * 2.0 * dz3;
        lhs(3,4,bb,k) = tmp1 * 2.0 * njac(3,4,k);
        lhs(3,5,bb,k) = tmp1 * 2.0 * njac(3,5,k);

        lhs(4,1,bb,k) = tmp1 * 2.0 * njac(4,1,k);
        lhs(4,2,bb,k) = tmp1 * 2.0 * njac(4,2,k);
        lhs(4,3,bb,k) = tmp1 * 2.0 * njac(4,3,k);
        lhs(4,4,bb,k) = 1.0
              + tmp1 * 2.0 * njac(4,4,k)
              + tmp1 * 2.0 * dz4;
        lhs(4,5,bb,k) = tmp1 * 2.0 * njac(4,5,k);

        lhs(5,1,bb,k) = tmp1 * 2.0 * njac(5,1,k);
        lhs(5,2,bb,k) = tmp1 * 2.0 * njac(5,2,k);
        lhs(5,3,bb,k) = tmp1 * 2.0 * njac(5,3,k);
        lhs(5,4,bb,k) = tmp1 * 2.0 * njac(5,4,k);
        lhs(5,5,bb,k) = 1.0
              + tmp1 * 2.0 * njac(5,5,k)
              + tmp1 * 2.0 * dz5;

        lhs(1,1,cc,k) =  tmp2 * fjac(1,1,k+1)
              - tmp1 * njac(1,1,k+1)
              - tmp1 * dz1;
        lhs(1,2,cc,k) =  tmp2 * fjac(1,2,k+1)
              - tmp1 * njac(1,2,k+1);
        lhs(1,3,cc,k) =  tmp2 * fjac(1,3,k+1)
              - tmp1 * njac(1,3,k+1);
        lhs(1,4,cc,k) =  tmp2 * fjac(1,4,k+1)
              - tmp1 * njac(1,4,k+1);
        lhs(1,5,cc,k) =  tmp2 * fjac(1,5,k+1)
              - tmp1 * njac(1,5,k+1);

        lhs(2,1,cc,k) =  tmp2 * fjac(2,1,k+1)
              - tmp1 * njac(2,1,k+1);
        lhs(2,2,cc,k) =  tmp2 * fjac(2,2,k+1)
              - tmp1 * njac(2,2,k+1)
              - tmp1 * dz2;
        lhs(2,3,cc,k) =  tmp2 * fjac(2,3,k+1)
              - tmp1 * njac(2,3,k+1);
        lhs(2,4,cc,k) =  tmp2 * fjac(2,4,k+1)
              - tmp1 * njac(2,4,k+1);
        lhs(2,5,cc,k) =  tmp2 * fjac(2,5,k+1)
              - tmp1 * njac(2,5,k+1);

        lhs(3,1,cc,k) =  tmp2 * fjac(3,1,k+1)
              - tmp1 * njac(3,1,k+1);
        lhs(3,2,cc,k) =  tmp2 * fjac(3,2,k+1)
              - tmp1 * njac(3,2,k+1);
        lhs(3,3,cc,k) =  tmp2 * fjac(3,3,k+1)
              - tmp1 * njac(3,3,k+1)
              - tmp1 * dz3;
        lhs(3,4,cc,k) =  tmp2 * fjac(3,4,k+1)
              - tmp1 * njac(3,4,k+1);
        lhs(3,5,cc,k) =  tmp2 * fjac(3,5,k+1)
              - tmp1 * njac(3,5,k+1);

        lhs(4,1,cc,k) =  tmp2 * fjac(4,1,k+1)
              - tmp1 * njac(4,1,k+1);
        lhs(4,2,cc,k) =  tmp2 * fjac(4,2,k+1)
              - tmp1 * njac(4,2,k+1);
        lhs(4,3,cc,k) =  tmp2 * fjac(4,3,k+1)
              - tmp1 * njac(4,3,k+1);
        lhs(4,4,cc,k) =  tmp2 * fjac(4,4,k+1)
              - tmp1 * njac(4,4,k+1)
              - tmp1 * dz4;
        lhs(4,5,cc,k) =  tmp2 * fjac(4,5,k+1)
              - tmp1 * njac(4,5,k+1);

        lhs(5,1,cc,k) =  tmp2 * fjac(5,1,k+1)
              - tmp1 * njac(5,1,k+1);
        lhs(5,2,cc,k) =  tmp2 * fjac(5,2,k+1)
              - tmp1 * njac(5,2,k+1);
        lhs(5,3,cc,k) =  tmp2 * fjac(5,3,k+1)
              - tmp1 * njac(5,3,k+1);
        lhs(5,4,cc,k) =  tmp2 * fjac(5,4,k+1)
              - tmp1 * njac(5,4,k+1);
        lhs(5,5,cc,k) =  tmp2 * fjac(5,5,k+1)
              - tmp1 * njac(5,5,k+1)
              - tmp1 * dz5;

      } // enddo

      // do k = 0, ksize
      for (SizeT k = 0; k <= ksize; ++k) {
          rtmp(1,k) = rhs(1,i,j,k);
          rtmp(2,k) = rhs(2,i,j,k);
          rtmp(3,k) = rhs(3,i,j,k);
          rtmp(4,k) = rhs(4,i,j,k);
          rtmp(5,k) = rhs(5,i,j,k);
      } // enddo

//---------------------------------------------------------------------
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     performs gaussian elimination on this cell.
//
//     assumes that unpacking routines for non-first cells
//     preload C' and rhs' from previous cell.
//
//     assumed send happens outside this routine, but that
//     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     outer most do loops - sweeping in i direction
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     multiply c(i,j,0) by b_inverse and copy back to c
//     multiply rhs(0) by b_inverse(0) and copy to rhs
//---------------------------------------------------------------------
        binvcrhs(&lhs(1,1,bb,0), &lhs(1,1,cc,0), &rtmp(1,0));


//---------------------------------------------------------------------
//     begin inner most do loop
//     do all the elements of the cell unless last
//---------------------------------------------------------------------
        // do k=1,ksize-1
        for (SizeT k = 1; k <= ksize-1; ++k) {

//---------------------------------------------------------------------
//     subtract A*lhs_vector(k-1) from lhs_vector(k)
//
//     rhs(k) = rhs(k) - A*rhs(k-1)
//---------------------------------------------------------------------
          matvec_sub(&lhs(1,1,aa,k), &rtmp(1,k-1), &rtmp(1,k));

//---------------------------------------------------------------------
//     B(k) = B(k) - C(k-1)*A(k)
//     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k)
//---------------------------------------------------------------------
          matmul_sub(&lhs(1,1,aa,k), &lhs(1,1,cc,k-1), &lhs(1,1,bb,k));

//---------------------------------------------------------------------
//     multiply c(i,j,k) by b_inverse and copy back to c
//     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
//---------------------------------------------------------------------
          binvcrhs(&lhs(1,1,bb,k), &lhs(1,1,cc,k), &rtmp(1,k));

        } // enddo

//---------------------------------------------------------------------
//     Now finish up special cases for last cell
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
//---------------------------------------------------------------------
        matvec_sub(&lhs(1,1,aa,ksize), &rtmp(1,ksize-1), &rtmp(1,ksize));

//---------------------------------------------------------------------
//     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
//     call matmul_sub(aa,i,j,ksize,c,
//     $              cc,i,j,ksize-1,c,bb,i,j,ksize)
//---------------------------------------------------------------------
        matmul_sub(&lhs(1,1,aa,ksize), &lhs(1,1,cc,ksize-1), &lhs(1,1,bb,ksize));

//---------------------------------------------------------------------
//     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
//---------------------------------------------------------------------
        binvrhs(&lhs(1,1,bb,ksize), &rtmp(1,ksize));


//---------------------------------------------------------------------
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     back solve: if last cell, then generate U(ksize)=rhs(ksize)
//     else assume U(ksize) is loaded in un pack backsub_info
//     so just use it
//     after call u(kstart) will be sent to next cell
//---------------------------------------------------------------------

        // do k=ksize-1,0,-1
        for (SizeT k = ksize-1; k >= 0; --k) {
          // do m=1,BLOCK_SIZE
          for (int m = 1; m <= BLOCK_SIZE; ++m) {
            // do n=1,BLOCK_SIZE
            for (int n = 1; n <= BLOCK_SIZE; ++n) {
              rtmp(m,k) = rtmp(m,k)
                   - lhs(m,n,cc,k)*rtmp(n,k+1);
            } // enddo
            rhs(m,i,j,k) = rtmp(m,k);
          } // enddo
        } // enddo

    } // enddo
  } // enddo
  });
//!$OMP END DO nowait
} // omp parallel
//!$OMP END PARALLEL
  if (timeron) timer_stop(t_zsolve);

  return;

}
/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_solve(rho_i, qs, square, u, rhs, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     Performs line solves in Z direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix,
c     and then performing back substitution to solve for the unknow
c     vectors of each line.
c
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c---------------------------------------------------------------------

      include 'header.h'
      include 'work_lhs.h'

      integer nx, nxmax, ny, nz
      double precision rho_i(  0:nxmax-1,0:ny-1,0:nz-1),
     $                 qs    ( 0:nxmax-1,0:ny-1,0:nz-1),
     $                 square( 0:nxmax-1,0:ny-1,0:nz-1),
     $                 u    (5,0:nxmax-1,0:ny-1,0:nz-1),
     $                 rhs  (5,0:nxmax-1,0:ny-1,0:nz-1)

      integer i, j, k, m, n, ksize

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      if (timeron) call timer_start(t_zsolve)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     This function computes the left hand side for the three z-factors
c---------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,m,k,i,j,ksize)
!$OMP&  SHARED(dz5,dz4,dz3,dz2,dz1,tz2,tz1,dt,c1345,c4,c3,con43,c3c4,c1,
!$OMP&         c2,nx,ny,nz)
      ksize = nz-1

c---------------------------------------------------------------------
c     Compute the indices for storing the block-diagonal matrix;
c     determine c (labeled f) and s jacobians
c---------------------------------------------------------------------
!$OMP DO
      do j = 1, ny-2
         do i = 1, nx-2
            do k = 0, ksize

               tmp1 = 1.d0 / u(1,i,j,k)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               fjac(1,1,k) = 0.d0
               fjac(1,2,k) = 0.d0
               fjac(1,3,k) = 0.d0
               fjac(1,4,k) = 1.d0
               fjac(1,5,k) = 0.d0

               fjac(2,1,k) = - ( u(2,i,j,k)*u(4,i,j,k) )
     >              * tmp2
               fjac(2,2,k) = u(4,i,j,k) * tmp1
               fjac(2,3,k) = 0.d0
               fjac(2,4,k) = u(2,i,j,k) * tmp1
               fjac(2,5,k) = 0.d0

               fjac(3,1,k) = - ( u(3,i,j,k)*u(4,i,j,k) )
     >              * tmp2
               fjac(3,2,k) = 0.d0
               fjac(3,3,k) = u(4,i,j,k) * tmp1
               fjac(3,4,k) = u(3,i,j,k) * tmp1
               fjac(3,5,k) = 0.d0

               fjac(4,1,k) = - (u(4,i,j,k)*u(4,i,j,k) * tmp2 )
     >              + c2 * qs(i,j,k)
               fjac(4,2,k) = - c2 *  u(2,i,j,k) * tmp1
               fjac(4,3,k) = - c2 *  u(3,i,j,k) * tmp1
               fjac(4,4,k) = ( 2.d0 - c2 )
     >              *  u(4,i,j,k) * tmp1
               fjac(4,5,k) = c2

               fjac(5,1,k) = ( c2 * 2.0d0 * square(i,j,k)
     >              - c1 * u(5,i,j,k) )
     >              * u(4,i,j,k) * tmp2
               fjac(5,2,k) = - c2 * ( u(2,i,j,k)*u(4,i,j,k) )
     >              * tmp2
               fjac(5,3,k) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
     >              * tmp2
               fjac(5,4,k) = c1 * ( u(5,i,j,k) * tmp1 )
     >              - c2
     >              * ( qs(i,j,k)
     >              + u(4,i,j,k)*u(4,i,j,k) * tmp2 )
               fjac(5,5,k) = c1 * u(4,i,j,k) * tmp1

               njac(1,1,k) = 0.d0
               njac(1,2,k) = 0.d0
               njac(1,3,k) = 0.d0
               njac(1,4,k) = 0.d0
               njac(1,5,k) = 0.d0

               njac(2,1,k) = - c3c4 * tmp2 * u(2,i,j,k)
               njac(2,2,k) =   c3c4 * tmp1
               njac(2,3,k) =   0.d0
               njac(2,4,k) =   0.d0
               njac(2,5,k) =   0.d0

               njac(3,1,k) = - c3c4 * tmp2 * u(3,i,j,k)
               njac(3,2,k) =   0.d0
               njac(3,3,k) =   c3c4 * tmp1
               njac(3,4,k) =   0.d0
               njac(3,5,k) =   0.d0

               njac(4,1,k) = - con43 * c3c4 * tmp2 * u(4,i,j,k)
               njac(4,2,k) =   0.d0
               njac(4,3,k) =   0.d0
               njac(4,4,k) =   con43 * c3 * c4 * tmp1
               njac(4,5,k) =   0.d0

               njac(5,1,k) = - (  c3c4
     >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
     >              - ( c3c4 - c1345 ) * tmp3 * (u(3,i,j,k)**2)
     >              - ( con43 * c3c4
     >              - c1345 ) * tmp3 * (u(4,i,j,k)**2)
     >              - c1345 * tmp2 * u(5,i,j,k)

               njac(5,2,k) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k)
               njac(5,3,k) = (  c3c4 - c1345 ) * tmp2 * u(3,i,j,k)
               njac(5,4,k) = ( con43 * c3c4
     >              - c1345 ) * tmp2 * u(4,i,j,k)
               njac(5,5,k) = ( c1345 )* tmp1


            enddo

c---------------------------------------------------------------------
c     now jacobians set, so form left hand side in z direction
c---------------------------------------------------------------------
            call lhsinit(lhs, ksize)
            do k = 1, ksize-1

               tmp1 = dt * tz1
               tmp2 = dt * tz2

               lhs(1,1,aa,k) = - tmp2 * fjac(1,1,k-1)
     >              - tmp1 * njac(1,1,k-1)
     >              - tmp1 * dz1
               lhs(1,2,aa,k) = - tmp2 * fjac(1,2,k-1)
     >              - tmp1 * njac(1,2,k-1)
               lhs(1,3,aa,k) = - tmp2 * fjac(1,3,k-1)
     >              - tmp1 * njac(1,3,k-1)
               lhs(1,4,aa,k) = - tmp2 * fjac(1,4,k-1)
     >              - tmp1 * njac(1,4,k-1)
               lhs(1,5,aa,k) = - tmp2 * fjac(1,5,k-1)
     >              - tmp1 * njac(1,5,k-1)

               lhs(2,1,aa,k) = - tmp2 * fjac(2,1,k-1)
     >              - tmp1 * njac(2,1,k-1)
               lhs(2,2,aa,k) = - tmp2 * fjac(2,2,k-1)
     >              - tmp1 * njac(2,2,k-1)
     >              - tmp1 * dz2
               lhs(2,3,aa,k) = - tmp2 * fjac(2,3,k-1)
     >              - tmp1 * njac(2,3,k-1)
               lhs(2,4,aa,k) = - tmp2 * fjac(2,4,k-1)
     >              - tmp1 * njac(2,4,k-1)
               lhs(2,5,aa,k) = - tmp2 * fjac(2,5,k-1)
     >              - tmp1 * njac(2,5,k-1)

               lhs(3,1,aa,k) = - tmp2 * fjac(3,1,k-1)
     >              - tmp1 * njac(3,1,k-1)
               lhs(3,2,aa,k) = - tmp2 * fjac(3,2,k-1)
     >              - tmp1 * njac(3,2,k-1)
               lhs(3,3,aa,k) = - tmp2 * fjac(3,3,k-1)
     >              - tmp1 * njac(3,3,k-1)
     >              - tmp1 * dz3
               lhs(3,4,aa,k) = - tmp2 * fjac(3,4,k-1)
     >              - tmp1 * njac(3,4,k-1)
               lhs(3,5,aa,k) = - tmp2 * fjac(3,5,k-1)
     >              - tmp1 * njac(3,5,k-1)

               lhs(4,1,aa,k) = - tmp2 * fjac(4,1,k-1)
     >              - tmp1 * njac(4,1,k-1)
               lhs(4,2,aa,k) = - tmp2 * fjac(4,2,k-1)
     >              - tmp1 * njac(4,2,k-1)
               lhs(4,3,aa,k) = - tmp2 * fjac(4,3,k-1)
     >              - tmp1 * njac(4,3,k-1)
               lhs(4,4,aa,k) = - tmp2 * fjac(4,4,k-1)
     >              - tmp1 * njac(4,4,k-1)
     >              - tmp1 * dz4
               lhs(4,5,aa,k) = - tmp2 * fjac(4,5,k-1)
     >              - tmp1 * njac(4,5,k-1)

               lhs(5,1,aa,k) = - tmp2 * fjac(5,1,k-1)
     >              - tmp1 * njac(5,1,k-1)
               lhs(5,2,aa,k) = - tmp2 * fjac(5,2,k-1)
     >              - tmp1 * njac(5,2,k-1)
               lhs(5,3,aa,k) = - tmp2 * fjac(5,3,k-1)
     >              - tmp1 * njac(5,3,k-1)
               lhs(5,4,aa,k) = - tmp2 * fjac(5,4,k-1)
     >              - tmp1 * njac(5,4,k-1)
               lhs(5,5,aa,k) = - tmp2 * fjac(5,5,k-1)
     >              - tmp1 * njac(5,5,k-1)
     >              - tmp1 * dz5

               lhs(1,1,bb,k) = 1.d0
     >              + tmp1 * 2.d0 * njac(1,1,k)
     >              + tmp1 * 2.d0 * dz1
               lhs(1,2,bb,k) = tmp1 * 2.d0 * njac(1,2,k)
               lhs(1,3,bb,k) = tmp1 * 2.d0 * njac(1,3,k)
               lhs(1,4,bb,k) = tmp1 * 2.d0 * njac(1,4,k)
               lhs(1,5,bb,k) = tmp1 * 2.d0 * njac(1,5,k)

               lhs(2,1,bb,k) = tmp1 * 2.d0 * njac(2,1,k)
               lhs(2,2,bb,k) = 1.d0
     >              + tmp1 * 2.d0 * njac(2,2,k)
     >              + tmp1 * 2.d0 * dz2
               lhs(2,3,bb,k) = tmp1 * 2.d0 * njac(2,3,k)
               lhs(2,4,bb,k) = tmp1 * 2.d0 * njac(2,4,k)
               lhs(2,5,bb,k) = tmp1 * 2.d0 * njac(2,5,k)

               lhs(3,1,bb,k) = tmp1 * 2.d0 * njac(3,1,k)
               lhs(3,2,bb,k) = tmp1 * 2.d0 * njac(3,2,k)
               lhs(3,3,bb,k) = 1.d0
     >              + tmp1 * 2.d0 * njac(3,3,k)
     >              + tmp1 * 2.d0 * dz3
               lhs(3,4,bb,k) = tmp1 * 2.d0 * njac(3,4,k)
               lhs(3,5,bb,k) = tmp1 * 2.d0 * njac(3,5,k)

               lhs(4,1,bb,k) = tmp1 * 2.d0 * njac(4,1,k)
               lhs(4,2,bb,k) = tmp1 * 2.d0 * njac(4,2,k)
               lhs(4,3,bb,k) = tmp1 * 2.d0 * njac(4,3,k)
               lhs(4,4,bb,k) = 1.d0
     >              + tmp1 * 2.d0 * njac(4,4,k)
     >              + tmp1 * 2.d0 * dz4
               lhs(4,5,bb,k) = tmp1 * 2.d0 * njac(4,5,k)

               lhs(5,1,bb,k) = tmp1 * 2.d0 * njac(5,1,k)
               lhs(5,2,bb,k) = tmp1 * 2.d0 * njac(5,2,k)
               lhs(5,3,bb,k) = tmp1 * 2.d0 * njac(5,3,k)
               lhs(5,4,bb,k) = tmp1 * 2.d0 * njac(5,4,k)
               lhs(5,5,bb,k) = 1.d0
     >              + tmp1 * 2.d0 * njac(5,5,k)
     >              + tmp1 * 2.d0 * dz5

               lhs(1,1,cc,k) =  tmp2 * fjac(1,1,k+1)
     >              - tmp1 * njac(1,1,k+1)
     >              - tmp1 * dz1
               lhs(1,2,cc,k) =  tmp2 * fjac(1,2,k+1)
     >              - tmp1 * njac(1,2,k+1)
               lhs(1,3,cc,k) =  tmp2 * fjac(1,3,k+1)
     >              - tmp1 * njac(1,3,k+1)
               lhs(1,4,cc,k) =  tmp2 * fjac(1,4,k+1)
     >              - tmp1 * njac(1,4,k+1)
               lhs(1,5,cc,k) =  tmp2 * fjac(1,5,k+1)
     >              - tmp1 * njac(1,5,k+1)

               lhs(2,1,cc,k) =  tmp2 * fjac(2,1,k+1)
     >              - tmp1 * njac(2,1,k+1)
               lhs(2,2,cc,k) =  tmp2 * fjac(2,2,k+1)
     >              - tmp1 * njac(2,2,k+1)
     >              - tmp1 * dz2
               lhs(2,3,cc,k) =  tmp2 * fjac(2,3,k+1)
     >              - tmp1 * njac(2,3,k+1)
               lhs(2,4,cc,k) =  tmp2 * fjac(2,4,k+1)
     >              - tmp1 * njac(2,4,k+1)
               lhs(2,5,cc,k) =  tmp2 * fjac(2,5,k+1)
     >              - tmp1 * njac(2,5,k+1)

               lhs(3,1,cc,k) =  tmp2 * fjac(3,1,k+1)
     >              - tmp1 * njac(3,1,k+1)
               lhs(3,2,cc,k) =  tmp2 * fjac(3,2,k+1)
     >              - tmp1 * njac(3,2,k+1)
               lhs(3,3,cc,k) =  tmp2 * fjac(3,3,k+1)
     >              - tmp1 * njac(3,3,k+1)
     >              - tmp1 * dz3
               lhs(3,4,cc,k) =  tmp2 * fjac(3,4,k+1)
     >              - tmp1 * njac(3,4,k+1)
               lhs(3,5,cc,k) =  tmp2 * fjac(3,5,k+1)
     >              - tmp1 * njac(3,5,k+1)

               lhs(4,1,cc,k) =  tmp2 * fjac(4,1,k+1)
     >              - tmp1 * njac(4,1,k+1)
               lhs(4,2,cc,k) =  tmp2 * fjac(4,2,k+1)
     >              - tmp1 * njac(4,2,k+1)
               lhs(4,3,cc,k) =  tmp2 * fjac(4,3,k+1)
     >              - tmp1 * njac(4,3,k+1)
               lhs(4,4,cc,k) =  tmp2 * fjac(4,4,k+1)
     >              - tmp1 * njac(4,4,k+1)
     >              - tmp1 * dz4
               lhs(4,5,cc,k) =  tmp2 * fjac(4,5,k+1)
     >              - tmp1 * njac(4,5,k+1)

               lhs(5,1,cc,k) =  tmp2 * fjac(5,1,k+1)
     >              - tmp1 * njac(5,1,k+1)
               lhs(5,2,cc,k) =  tmp2 * fjac(5,2,k+1)
     >              - tmp1 * njac(5,2,k+1)
               lhs(5,3,cc,k) =  tmp2 * fjac(5,3,k+1)
     >              - tmp1 * njac(5,3,k+1)
               lhs(5,4,cc,k) =  tmp2 * fjac(5,4,k+1)
     >              - tmp1 * njac(5,4,k+1)
               lhs(5,5,cc,k) =  tmp2 * fjac(5,5,k+1)
     >              - tmp1 * njac(5,5,k+1)
     >              - tmp1 * dz5

            enddo

            do k = 0, ksize
               rtmp(1,k) = rhs(1,i,j,k)
               rtmp(2,k) = rhs(2,i,j,k)
               rtmp(3,k) = rhs(3,i,j,k)
               rtmp(4,k) = rhs(4,i,j,k)
               rtmp(5,k) = rhs(5,i,j,k)
            enddo

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     performs gaussian elimination on this cell.
c
c     assumes that unpacking routines for non-first cells
c     preload C' and rhs' from previous cell.
c
c     assumed send happens outside this routine, but that
c     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     multiply c(i,j,0) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c---------------------------------------------------------------------
            call binvcrhs( lhs(1,1,bb,0),
     >                        lhs(1,1,cc,0),
     >                        rtmp(1,0) )


c---------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last
c---------------------------------------------------------------------
            do k=1,ksize-1

c---------------------------------------------------------------------
c     subtract A*lhs_vector(k-1) from lhs_vector(k)
c
c     rhs(k) = rhs(k) - A*rhs(k-1)
c---------------------------------------------------------------------
               call matvec_sub(lhs(1,1,aa,k),
     >                         rtmp(1,k-1),rtmp(1,k))

c---------------------------------------------------------------------
c     B(k) = B(k) - C(k-1)*A(k)
c     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k)
c---------------------------------------------------------------------
               call matmul_sub(lhs(1,1,aa,k),
     >                         lhs(1,1,cc,k-1),
     >                         lhs(1,1,bb,k))

c---------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
c---------------------------------------------------------------------
               call binvcrhs( lhs(1,1,bb,k),
     >                        lhs(1,1,cc,k),
     >                        rtmp(1,k) )

            enddo

c---------------------------------------------------------------------
c     Now finish up special cases for last cell
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
c---------------------------------------------------------------------
            call matvec_sub(lhs(1,1,aa,ksize),
     >                         rtmp(1,ksize-1),rtmp(1,ksize))

c---------------------------------------------------------------------
c     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
c     call matmul_sub(aa,i,j,ksize,c,
c     $              cc,i,j,ksize-1,c,bb,i,j,ksize)
c---------------------------------------------------------------------
            call matmul_sub(lhs(1,1,aa,ksize),
     >                         lhs(1,1,cc,ksize-1),
     >                         lhs(1,1,bb,ksize))

c---------------------------------------------------------------------
c     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
c---------------------------------------------------------------------
            call binvrhs( lhs(1,1,bb,ksize),
     >                       rtmp(1,ksize) )


c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     back solve: if last cell, then generate U(ksize)=rhs(ksize)
c     else assume U(ksize) is loaded in un pack backsub_info
c     so just use it
c     after call u(kstart) will be sent to next cell
c---------------------------------------------------------------------

            do k=ksize-1,0,-1
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rtmp(m,k) = rtmp(m,k)
     >                    - lhs(m,n,cc,k)*rtmp(n,k+1)
                  enddo
                  rhs(m,i,j,k) = rtmp(m,k)
               enddo
            enddo

         enddo
      enddo
!$OMP END DO nowait
!$OMP END PARALLEL
      if (timeron) call timer_stop(t_zsolve)

      return
      end

*/
