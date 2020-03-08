
#include "NArray.h"
#include "initialize.h"

template<typename SizeT, typename ValueT>
void y_solve(
  ValueT* rho_i_ptr,
  ValueT* qs_ptr,
  ValueT* square_ptr,
  ValueT* u_ptr,
  ValueT* rhs_ptr,
  SizeT nx,
  SizeT nxmax,
  SizeT ny,
  SizeT nz)
{

  ArrayViewT3<ValueT, SizeT> rho_i(rho_i_ptr, nxmax, ny, nz);
  ArrayViewT3<ValueT, SizeT> qs(qs_ptr, nxmax, ny, nz);
  ArrayViewT3<ValueT, SizeT> square(square_ptr, nxmax, ny, nz);
  ArrayViewT4<ValueT, SizeT> u(u_ptr, 5, nxmax, ny, nz);
  ArrayViewT4<ValueT, SizeT> rhs(rhs_ptr, 5, nxmax, ny, nz);

/*
c---------------------------------------------------------------------
c     Performs line solves in Y direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix,
c     and then performing back substitution to solve for the unknow
c     vectors of each line.
c
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c---------------------------------------------------------------------
*/

#pragma omp parallel
{
  // thread-local temporary variables, allocated on the stack
  StaticArrayT3<double, size_type, 5, 5,    problem_size+1, 1, 1, 0> fjac/*(5, 5, problem_size+1)*/;
  StaticArrayT3<double, size_type, 5, 5,    problem_size+1, 1, 1, 0> njac/*(5, 5, problem_size+1)*/;
  StaticArrayT4<double, size_type, 5, 5, 3, problem_size+1, 1, 1, 1, 0> lhs /*(5, 5, 3, problem_size+1)*/;

/*
c---------------------------------------------------------------------
c     This function computes the left hand side for the three y-factors
c---------------------------------------------------------------------
*/

  SizeT jsize = ny-1;

  if (timeron) timer_start(t_ysolve);

/*
c---------------------------------------------------------------------
c     Compute the indices for storing the tri-diagonal matrix;
c     determine a (labeled f) and n jacobians for cell c
c---------------------------------------------------------------------
*/
  //do k = 1, nz-2
#pragma omp for
  for (SizeT k = 1; k <= nz-2; ++k) {
    //do i = 1, nx-2
    for (SizeT i = 1; i <= nx-2; ++i) {
      // do j = 0, jsize
      for (SizeT j = 0; j <= jsize; ++j) {
        ValueT tmp1, tmp2, tmp3;

        tmp1 = rho_i(i,j,k);
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;

        fjac(1,1,j) = 0.0;
        fjac(1,2,j) = 0.0;
        fjac(1,3,j) = 1.0;
        fjac(1,4,j) = 0.0;
        fjac(1,5,j) = 0.0;

        fjac(2,1,j) = - ( u(2,i,j,k)*u(3,i,j,k) )
            * tmp2;
        fjac(2,2,j) = u(3,i,j,k) * tmp1;
        fjac(2,3,j) = u(2,i,j,k) * tmp1;
        fjac(2,4,j) = 0.0;
        fjac(2,5,j) = 0.0;

        fjac(3,1,j) = - ( u(3,i,j,k)*u(3,i,j,k)*tmp2)
            + c2 * qs(i,j,k);
        fjac(3,2,j) = - c2 *  u(2,i,j,k) * tmp1;
        fjac(3,3,j) = ( 2.0 - c2 )
            *  u(3,i,j,k) * tmp1;
        fjac(3,4,j) = - c2 * u(4,i,j,k) * tmp1;
        fjac(3,5,j) = c2;

        fjac(4,1,j) = - ( u(3,i,j,k)*u(4,i,j,k) )
            * tmp2;
        fjac(4,2,j) = 0.0;
        fjac(4,3,j) = u(4,i,j,k) * tmp1;
        fjac(4,4,j) = u(3,i,j,k) * tmp1;
        fjac(4,5,j) = 0.0;

        fjac(5,1,j) = ( c2 * 2.0 * square(i,j,k)
            - c1 * u(5,i,j,k) )
            * u(3,i,j,k) * tmp2;
        fjac(5,2,j) = - c2 * u(2,i,j,k)*u(3,i,j,k)
            * tmp2;
        fjac(5,3,j) = c1 * u(5,i,j,k) * tmp1
            - c2
            * ( qs(i,j,k)
            + u(3,i,j,k)*u(3,i,j,k) * tmp2 );
        fjac(5,4,j) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
            * tmp2;
        fjac(5,5,j) = c1 * u(3,i,j,k) * tmp1;

        njac(1,1,j) = 0.0;
        njac(1,2,j) = 0.0;
        njac(1,3,j) = 0.0;
        njac(1,4,j) = 0.0;
        njac(1,5,j) = 0.0;

        njac(2,1,j) = - c3c4 * tmp2 * u(2,i,j,k);
        njac(2,2,j) =   c3c4 * tmp1;
        njac(2,3,j) =   0.0;
        njac(2,4,j) =   0.0;
        njac(2,5,j) =   0.0;

        njac(3,1,j) = - con43 * c3c4 * tmp2 * u(3,i,j,k);
        njac(3,2,j) =   0.0;
        njac(3,3,j) =   con43 * c3c4 * tmp1;
        njac(3,4,j) =   0.0;
        njac(3,5,j) =   0.0;

        njac(4,1,j) = - c3c4 * tmp2 * u(4,i,j,k);
        njac(4,2,j) =   0.0;
        njac(4,3,j) =   0.0;
        njac(4,4,j) =   c3c4 * tmp1;
        njac(4,5,j) =   0.0;

        njac(5,1,j) = - (  c3c4
            - c1345 ) * tmp3 * (u(2,i,j,k)*u(2,i,j,k))
            - ( con43 * c3c4
            - c1345 ) * tmp3 * (u(3,i,j,k)*u(3,i,j,k))
            - ( c3c4 - c1345 ) * tmp3 * (u(4,i,j,k)*u(4,i,j,k))
            - c1345 * tmp2 * u(5,i,j,k);

        njac(5,2,j) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k);
        njac(5,3,j) = ( con43 * c3c4
            - c1345 ) * tmp2 * u(3,i,j,k);
        njac(5,4,j) = ( c3c4 - c1345 ) * tmp2 * u(4,i,j,k);
        njac(5,5,j) = ( c1345 ) * tmp1;
        //enddo
      }

/*
c---------------------------------------------------------------------
c     now joacobians set, so form left hand side in y direction
c---------------------------------------------------------------------
*/
      lhsinit(lhs.begin(), jsize);
      //do j = 1, jsize-1
      for (int j = 1; j <= jsize-1; ++j) {

        ValueT tmp1 = dt * ty1;
        ValueT tmp2 = dt * ty2;

        lhs(1,1,aa,j) = - tmp2 * fjac(1,1,j-1)
               - tmp1 * njac(1,1,j-1)
               - tmp1 * dy1;
        lhs(1,2,aa,j) = - tmp2 * fjac(1,2,j-1)
               - tmp1 * njac(1,2,j-1);
        lhs(1,3,aa,j) = - tmp2 * fjac(1,3,j-1)
               - tmp1 * njac(1,3,j-1);
        lhs(1,4,aa,j) = - tmp2 * fjac(1,4,j-1)
               - tmp1 * njac(1,4,j-1);
        lhs(1,5,aa,j) = - tmp2 * fjac(1,5,j-1)
               - tmp1 * njac(1,5,j-1);

        lhs(2,1,aa,j) = - tmp2 * fjac(2,1,j-1)
               - tmp1 * njac(2,1,j-1);
        lhs(2,2,aa,j) = - tmp2 * fjac(2,2,j-1)
               - tmp1 * njac(2,2,j-1)
               - tmp1 * dy2;
        lhs(2,3,aa,j) = - tmp2 * fjac(2,3,j-1)
               - tmp1 * njac(2,3,j-1);
        lhs(2,4,aa,j) = - tmp2 * fjac(2,4,j-1)
               - tmp1 * njac(2,4,j-1);
        lhs(2,5,aa,j) = - tmp2 * fjac(2,5,j-1)
               - tmp1 * njac(2,5,j-1);

        lhs(3,1,aa,j) = - tmp2 * fjac(3,1,j-1)
               - tmp1 * njac(3,1,j-1);
        lhs(3,2,aa,j) = - tmp2 * fjac(3,2,j-1)
               - tmp1 * njac(3,2,j-1);
        lhs(3,3,aa,j) = - tmp2 * fjac(3,3,j-1)
               - tmp1 * njac(3,3,j-1)
               - tmp1 * dy3;
        lhs(3,4,aa,j) = - tmp2 * fjac(3,4,j-1)
               - tmp1 * njac(3,4,j-1);
        lhs(3,5,aa,j) = - tmp2 * fjac(3,5,j-1)
               - tmp1 * njac(3,5,j-1);

        lhs(4,1,aa,j) = - tmp2 * fjac(4,1,j-1)
               - tmp1 * njac(4,1,j-1);
        lhs(4,2,aa,j) = - tmp2 * fjac(4,2,j-1)
               - tmp1 * njac(4,2,j-1);
        lhs(4,3,aa,j) = - tmp2 * fjac(4,3,j-1)
               - tmp1 * njac(4,3,j-1);
        lhs(4,4,aa,j) = - tmp2 * fjac(4,4,j-1)
               - tmp1 * njac(4,4,j-1)
               - tmp1 * dy4;
        lhs(4,5,aa,j) = - tmp2 * fjac(4,5,j-1)
               - tmp1 * njac(4,5,j-1);

        lhs(5,1,aa,j) = - tmp2 * fjac(5,1,j-1)
               - tmp1 * njac(5,1,j-1);
        lhs(5,2,aa,j) = - tmp2 * fjac(5,2,j-1)
               - tmp1 * njac(5,2,j-1);
        lhs(5,3,aa,j) = - tmp2 * fjac(5,3,j-1)
               - tmp1 * njac(5,3,j-1);
        lhs(5,4,aa,j) = - tmp2 * fjac(5,4,j-1)
               - tmp1 * njac(5,4,j-1);
        lhs(5,5,aa,j) = - tmp2 * fjac(5,5,j-1)
               - tmp1 * njac(5,5,j-1)
               - tmp1 * dy5;

        lhs(1,1,bb,j) = 1.0
               + tmp1 * 2.0 * njac(1,1,j)
               + tmp1 * 2.0 * dy1;
        lhs(1,2,bb,j) = tmp1 * 2.0 * njac(1,2,j);
        lhs(1,3,bb,j) = tmp1 * 2.0 * njac(1,3,j);
        lhs(1,4,bb,j) = tmp1 * 2.0 * njac(1,4,j);
        lhs(1,5,bb,j) = tmp1 * 2.0 * njac(1,5,j);

        lhs(2,1,bb,j) = tmp1 * 2.0 * njac(2,1,j);
        lhs(2,2,bb,j) = 1.0
               + tmp1 * 2.0 * njac(2,2,j)
               + tmp1 * 2.0 * dy2;
        lhs(2,3,bb,j) = tmp1 * 2.0 * njac(2,3,j);
        lhs(2,4,bb,j) = tmp1 * 2.0 * njac(2,4,j);
        lhs(2,5,bb,j) = tmp1 * 2.0 * njac(2,5,j);

        lhs(3,1,bb,j) = tmp1 * 2.0 * njac(3,1,j);
        lhs(3,2,bb,j) = tmp1 * 2.0 * njac(3,2,j);
        lhs(3,3,bb,j) = 1.0
               + tmp1 * 2.0 * njac(3,3,j)
               + tmp1 * 2.0 * dy3;
        lhs(3,4,bb,j) = tmp1 * 2.0 * njac(3,4,j);
        lhs(3,5,bb,j) = tmp1 * 2.0 * njac(3,5,j);

        lhs(4,1,bb,j) = tmp1 * 2.0 * njac(4,1,j);
        lhs(4,2,bb,j) = tmp1 * 2.0 * njac(4,2,j);
        lhs(4,3,bb,j) = tmp1 * 2.0 * njac(4,3,j);
        lhs(4,4,bb,j) = 1.0
               + tmp1 * 2.0 * njac(4,4,j)
               + tmp1 * 2.0 * dy4;
        lhs(4,5,bb,j) = tmp1 * 2.0 * njac(4,5,j);

        lhs(5,1,bb,j) = tmp1 * 2.0 * njac(5,1,j);
        lhs(5,2,bb,j) = tmp1 * 2.0 * njac(5,2,j);
        lhs(5,3,bb,j) = tmp1 * 2.0 * njac(5,3,j);
        lhs(5,4,bb,j) = tmp1 * 2.0 * njac(5,4,j);
        lhs(5,5,bb,j) = 1.0
               + tmp1 * 2.0 * njac(5,5,j)
               + tmp1 * 2.0 * dy5;

        lhs(1,1,cc,j) =  tmp2 * fjac(1,1,j+1)
               - tmp1 * njac(1,1,j+1)
               - tmp1 * dy1;
        lhs(1,2,cc,j) =  tmp2 * fjac(1,2,j+1)
               - tmp1 * njac(1,2,j+1);
        lhs(1,3,cc,j) =  tmp2 * fjac(1,3,j+1)
               - tmp1 * njac(1,3,j+1);
        lhs(1,4,cc,j) =  tmp2 * fjac(1,4,j+1)
               - tmp1 * njac(1,4,j+1);
        lhs(1,5,cc,j) =  tmp2 * fjac(1,5,j+1)
               - tmp1 * njac(1,5,j+1);

        lhs(2,1,cc,j) =  tmp2 * fjac(2,1,j+1)
               - tmp1 * njac(2,1,j+1);
        lhs(2,2,cc,j) =  tmp2 * fjac(2,2,j+1)
               - tmp1 * njac(2,2,j+1)
               - tmp1 * dy2;
        lhs(2,3,cc,j) =  tmp2 * fjac(2,3,j+1)
               - tmp1 * njac(2,3,j+1);
        lhs(2,4,cc,j) =  tmp2 * fjac(2,4,j+1)
               - tmp1 * njac(2,4,j+1);
        lhs(2,5,cc,j) =  tmp2 * fjac(2,5,j+1)
               - tmp1 * njac(2,5,j+1);

        lhs(3,1,cc,j) =  tmp2 * fjac(3,1,j+1)
               - tmp1 * njac(3,1,j+1);
        lhs(3,2,cc,j) =  tmp2 * fjac(3,2,j+1)
               - tmp1 * njac(3,2,j+1);
        lhs(3,3,cc,j) =  tmp2 * fjac(3,3,j+1)
               - tmp1 * njac(3,3,j+1)
               - tmp1 * dy3;
        lhs(3,4,cc,j) =  tmp2 * fjac(3,4,j+1)
               - tmp1 * njac(3,4,j+1);
        lhs(3,5,cc,j) =  tmp2 * fjac(3,5,j+1)
               - tmp1 * njac(3,5,j+1);

        lhs(4,1,cc,j) =  tmp2 * fjac(4,1,j+1)
               - tmp1 * njac(4,1,j+1);
        lhs(4,2,cc,j) =  tmp2 * fjac(4,2,j+1)
               - tmp1 * njac(4,2,j+1);
        lhs(4,3,cc,j) =  tmp2 * fjac(4,3,j+1)
               - tmp1 * njac(4,3,j+1);
        lhs(4,4,cc,j) =  tmp2 * fjac(4,4,j+1)
               - tmp1 * njac(4,4,j+1)
               - tmp1 * dy4;
        lhs(4,5,cc,j) =  tmp2 * fjac(4,5,j+1)
               - tmp1 * njac(4,5,j+1);

        lhs(5,1,cc,j) =  tmp2 * fjac(5,1,j+1)
               - tmp1 * njac(5,1,j+1);
        lhs(5,2,cc,j) =  tmp2 * fjac(5,2,j+1)
               - tmp1 * njac(5,2,j+1);
        lhs(5,3,cc,j) =  tmp2 * fjac(5,3,j+1)
               - tmp1 * njac(5,3,j+1);
        lhs(5,4,cc,j) =  tmp2 * fjac(5,4,j+1)
               - tmp1 * njac(5,4,j+1);
        lhs(5,5,cc,j) =  tmp2 * fjac(5,5,j+1)
               - tmp1 * njac(5,5,j+1)
               - tmp1 * dy5;

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
//     c'(JMAX) and rhs'(JMAX) will be sent to next cell
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//     multiply c(i,0,k) by b_inverse and copy back to c
//     multiply rhs(0) by b_inverse(0) and copy to rhs
//---------------------------------------------------------------------
      binvcrhs( &lhs(1,1,bb,0), &lhs(1,1,cc,0), &rhs(1,i,0,k) );



//---------------------------------------------------------------------
//     begin inner most do loop
//     do all the elements of the cell unless last
//---------------------------------------------------------------------
      //do j=1,jsize-1
      for (SizeT j = 1; j <= jsize-1; ++j) {

//---------------------------------------------------------------------
//     subtract A*lhs_vector(j-1) from lhs_vector(j)
//
//     rhs(j) = rhs(j) - A*rhs(j-1)
//---------------------------------------------------------------------
        matvec_sub(&lhs(1,1,aa,j), &rhs(1,i,j-1,k), &rhs(1,i,j,k));

//---------------------------------------------------------------------
//     B(j) = B(j) - C(j-1)*A(j)
//---------------------------------------------------------------------
        matmul_sub(&lhs(1,1,aa,j), &lhs(1,1,cc,j-1), &lhs(1,1,bb,j));

//---------------------------------------------------------------------
//     multiply c(i,j,k) by b_inverse and copy back to c
//     multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
//---------------------------------------------------------------------
        binvcrhs(&lhs(1,1,bb,j), &lhs(1,1,cc,j), &rhs(1,i,j,k));

      } // enddo


//---------------------------------------------------------------------
//     rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
//---------------------------------------------------------------------
      matvec_sub(&lhs(1,1,aa,jsize), &rhs(1,i,jsize-1,k), &rhs(1,i,jsize,k));

//---------------------------------------------------------------------
//     B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
//     call matmul_sub(aa,i,jsize,k,c,
//     $              cc,i,jsize-1,k,c,bb,i,jsize,k)
//---------------------------------------------------------------------
      matmul_sub(&lhs(1,1,aa,jsize), &lhs(1,1,cc,jsize-1), &lhs(1,1,bb,jsize));

//---------------------------------------------------------------------
//     multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
//---------------------------------------------------------------------
      binvrhs(&lhs(1,1,bb,jsize), &rhs(1,i,jsize,k));


//---------------------------------------------------------------------
//     back solve: if last cell, then generate U(jsize)=rhs(jsize)
//     else assume U(jsize) is loaded in un pack backsub_info
//     so just use it
//     after call u(jstart) will be sent to next cell
//---------------------------------------------------------------------

      //do j=jsize-1,0,-1
      for (SizeT j = jsize-1; j >= 0; --j) {
        //do m=1,BLOCK_SIZE
        for (int m = 1; m <= BLOCK_SIZE; ++m) {
          // do n=1,BLOCK_SIZE
          for (int n = 1; n <= BLOCK_SIZE; ++n) {
            rhs(m,i,j,k) = rhs(m,i,j,k)
                    - lhs(m,n,cc,j)*rhs(n,i,j+1,k);
          } // enddo
        } // enddo
      } // enddo

    } // enddo
  } // enddo

} // omp parallel

  if (timeron) timer_stop(t_ysolve);


  return;

}
