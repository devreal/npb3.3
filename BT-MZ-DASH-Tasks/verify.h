#ifndef HAVE_VERIFY_H
#define HAVE_VERIFY_H

#include <mpi.h>

#include "NArray.h"
#include "mpi_stuff.h"
#include "header.h"
#include "error.h"

/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

        subroutine verify(no_time_steps, verified, num_zones,
     $                    rho_i, us, vs, ws, qs, square,
     $                    rhs, forcing, u, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c  verification routine
c---------------------------------------------------------------------
*/

template<typename ValueT, typename SizeT>
void verify(int no_time_steps, bool& verified, int num_zones,
            ArrayT1<ValueT, size_type>& rho_i,
            ArrayT1<ValueT, size_type>& us,
            ArrayT1<ValueT, size_type>& vs,
            ArrayT1<ValueT, size_type>& ws,
            ArrayT1<ValueT, size_type>& qs,
            ArrayT1<ValueT, size_type>& square,
            ArrayT1<ValueT, size_type>& rhs,
            ArrayT1<ValueT, size_type>& forcing,
            ArrayT1<ValueT, size_type>& u,
            ArrayT1<SizeT , size_type>& nx,
            ArrayT1<SizeT , size_type>& nxmax,
            ArrayT1<SizeT , size_type>& ny,
            ArrayT1<SizeT , size_type>& nz)
{
/*
        include 'header.h'
        include 'mpi_stuff.h'

        integer zone, num_zones
        double precision rho_i(*), us(*), vs(*), ws(*), qs(*),
     $                   square(*), rhs(*), forcing(*), u(*)

        double precision xcrref(5),xceref(5),xcrdif(5),xcedif(5),
     >                   epsilon, xce(5), xcr(5), dtref,
     $                   xce_sub(5), xcr_sub(5)
        integer m, no_time_steps, niterref, iz, ip
        integer nx(*), nxmax(*), ny(*), nz(*)
        logical verified
*/
  StaticArrayT1<ValueT, SizeT, 5, 1>
    xcrref,xceref,xcrdif,xcedif, xce, xcr, xce_sub, xcr_sub;

  ValueT dtref, epsilon;

  int niterref;

//---------------------------------------------------------------------
//   tolerance level
//---------------------------------------------------------------------
  epsilon = 1.0e-08;

//---------------------------------------------------------------------
//   compute the error norm and the residual norm, and exit if not printing
//---------------------------------------------------------------------

  //do m = 1, 5
  for (int m = 1; m <= 5; ++m) {
    xcr(m) = 0.0;
    xce(m) = 0.0;
  } // end do

  //do iz = 1, proc_num_zones
  for (int iz = 1; iz <= proc_num_zones; ++iz) {
    int zone = proc_zone_id(iz);
    error_norm (xce_sub.begin(), &u(start5(iz)), nx(zone), nxmax(zone), ny(zone), nz(zone));
    compute_rhs(&rho_i(start1(iz)), &us(start1(iz)),
                &vs(start1(iz)), &ws(start1(iz)),
                &qs(start1(iz)), &square(start1(iz)),
                &rhs(start5(iz)), &forcing(start5(iz)),
                &u(start5(iz)),
                nx(zone), nxmax(zone), ny(zone), nz(zone));
    dash::tasks::complete(true);
    rhs_norm(xcr_sub.begin(), &rhs(start5(iz)), nx(zone), nxmax(zone), ny(zone), nz(zone));

    // do m = 1, 5
    for (int m = 1; m <= 5; ++m) {
      xcr(m) = xcr(m) + xcr_sub(m) / dt;
      xce(m) = xce(m) + xce_sub(m);
    } // end do
  } //end do

  //do m = 1, 5
  for (int m = 1; m <= 5; ++m) {
    xcr_sub(m) = xcr(m);
    xce_sub(m) = xce(m);
  } // end do

  MPI_Reduce(xcr_sub.begin(), xcr.begin(), 5,
             MPI_DOUBLE, MPI_SUM, root, comm_setup);
  MPI_Reduce(xce_sub.begin(), xce.begin(), 5,
             MPI_DOUBLE, MPI_SUM, root, comm_setup);

  if (myid != root) return;

  verified = true;

  //do m = 1,5
  for (int m = 1; m <= 5; ++m) {
    xcrref(m) = 1.0;
    xceref(m) = 1.0;
  } // end do

//---------------------------------------------------------------------
//    reference data for class S
//---------------------------------------------------------------------
  if ( CLASS == 'S' ) {
    dtref = 1.0e-2;
    niterref = 60;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//---------------------------------------------------------------------
    xcrref(1) = 0.1047687395830e+04;
    xcrref(2) = 0.9419911314792e+02;
    xcrref(3) = 0.2124737403068e+03;
    xcrref(4) = 0.1422173591794e+03;
    xcrref(5) = 0.1135441572375e+04;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//---------------------------------------------------------------------
    xceref(1) = 0.1775416062982e+03;
    xceref(2) = 0.1875540250835e+02;
    xceref(3) = 0.3863334844506e+02;
    xceref(4) = 0.2634713890362e+02;
    xceref(5) = 0.1965566269675e+03;

//---------------------------------------------------------------------
//    reference data for class W
//---------------------------------------------------------------------
  } else if ( CLASS == 'W' ) {
    dtref = 0.8e-3;
    niterref = 200;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//---------------------------------------------------------------------
    xcrref(1) = 0.5562611195402e+05;
    xcrref(2) = 0.5151404119932e+04;
    xcrref(3) = 0.1080453907954e+05;
    xcrref(4) = 0.6576058591929e+04;
    xcrref(5) = 0.4528609293561e+05;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//---------------------------------------------------------------------
    xceref(1) = 0.7185154786403e+04;
    xceref(2) = 0.7040472738068e+03;
    xceref(3) = 0.1437035074443e+04;
    xceref(4) = 0.8570666307849e+03;
    xceref(5) = 0.5991235147368e+04;

//---------------------------------------------------------------------
//    reference data for class A
//---------------------------------------------------------------------
  } else if ( CLASS == 'A' ) {
    dtref = 0.8e-3;
    niterref = 200;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//---------------------------------------------------------------------
    xcrref(1) = 0.5536703889522e+05;
    xcrref(2) = 0.5077835038405e+04;
    xcrref(3) = 0.1067391361067e+05;
    xcrref(4) = 0.6441179694972e+04;
    xcrref(5) = 0.4371926324069e+05;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//---------------------------------------------------------------------
    xceref(1) = 0.6716797714343e+04;
    xceref(2) = 0.6512687902160e+03;
    xceref(3) = 0.1332930740128e+04;
    xceref(4) = 0.7848302089180e+03;
    xceref(5) = 0.5429053878818e+04;

//---------------------------------------------------------------------
//    reference data for class B
//---------------------------------------------------------------------
  } else if ( CLASS == 'B' ) {
    dtref = 3.0e-4;
    niterref = 200;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//---------------------------------------------------------------------
    xcrref(1) = 0.4461388343844e+06;
    xcrref(2) = 0.3799759138035e+05;
    xcrref(3) = 0.8383296623970e+05;
    xcrref(4) = 0.5301970201273e+05;
    xcrref(5) = 0.3618106851311e+06;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//---------------------------------------------------------------------
    xceref(1) = 0.4496733567600e+05;
    xceref(2) = 0.3892068540524e+04;
    xceref(3) = 0.8763825844217e+04;
    xceref(4) = 0.5599040091792e+04;
    xceref(5) = 0.4082652045598e+05;

//---------------------------------------------------------------------
//    reference data class C
//---------------------------------------------------------------------
  } else if ( CLASS == 'C' ) {
    dtref = 1.0e-4;
    niterref = 200;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//---------------------------------------------------------------------
    xcrref(1) = 0.3457703287806e+07;
    xcrref(2) = 0.3213621375929e+06;
    xcrref(3) = 0.7002579656870e+06;
    xcrref(4) = 0.4517459627471e+06;
    xcrref(5) = 0.2818715870791e+07;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//---------------------------------------------------------------------
    xceref(1) = 0.2059106993570e+06;
    xceref(2) = 0.1680761129461e+05;
    xceref(3) = 0.4080731640795e+05;
    xceref(4) = 0.2836541076778e+05;
    xceref(5) = 0.2136807610771e+06;

//---------------------------------------------------------------------
//    reference data class D
//---------------------------------------------------------------------
  } else if ( CLASS == 'D' ) {
    dtref = 2.0e-5;
    niterref = 250;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//---------------------------------------------------------------------
    xcrref(1) = 0.4250417034981e+08;
    xcrref(2) = 0.4293882192175e+07;
    xcrref(3) = 0.9121841878270e+07;
    xcrref(4) = 0.6201357771439e+07;
    xcrref(5) = 0.3474801891304e+08;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//---------------------------------------------------------------------
    xceref(1) = 0.9462418484583e+06;
    xceref(2) = 0.7884728947105e+05;
    xceref(3) = 0.1902874461259e+06;
    xceref(4) = 0.1361858029909e+06;
    xceref(5) = 0.9816489456253e+06;

//---------------------------------------------------------------------
//    reference data class E
//---------------------------------------------------------------------
  } else if ( CLASS == 'E' ) {
    dtref = 4.0e-6;
    niterref = 250;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//---------------------------------------------------------------------
    xcrref(1) = 0.5744815962469e+09;
    xcrref(2) = 0.6088696479719e+08;
    xcrref(3) = 0.1276325224438e+09;
    xcrref(4) = 0.8947040105616e+08;
    xcrref(5) = 0.4726115284807e+09;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//---------------------------------------------------------------------
    xceref(1) = 0.4114447054461e+07;
    xceref(2) = 0.3570776728190e+06;
    xceref(3) = 0.8465106191458e+06;
    xceref(4) = 0.6147182273817e+06;
    xceref(5) = 0.4238908025163e+07;

//---------------------------------------------------------------------
//    reference data class F
//---------------------------------------------------------------------
  } else if ( CLASS == 'F' ) {
    dtref = 1.0e-6;
    niterref = 250;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//---------------------------------------------------------------------
    xcrref(1) = 0.6524078317845e+10;
    xcrref(2) = 0.7020439279514e+09;
    xcrref(3) = 0.1467588422194e+10;
    xcrref(4) = 0.1042973064137e+10;
    xcrref(5) = 0.5411102201141e+10;

//---------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//---------------------------------------------------------------------
    xceref(1) = 0.1708795375347e+08;
    xceref(2) = 0.1514359936802e+07;
    xceref(3) = 0.3552878359250e+07;
    xceref(4) = 0.2594549582184e+07;
    xceref(5) = 0.1749809607845e+08;

    if (no_time_steps == 25) {
      niterref = 25;
      xcrref(1) = 0.3565049484400e+11;
      xcrref(2) = 0.3752029586145e+10;
      xcrref(3) = 0.7805935552197e+10;
      xcrref(4) = 0.5685995438056e+10;
      xcrref(5) = 0.2908811276266e+11;

      xceref(1) = 0.1805995755490e+08;
      xceref(2) = 0.1632306899424e+07;
      xceref(3) = 0.3778610439036e+07;
      xceref(4) = 0.2749319818549e+07;
      xceref(5) = 0.1814401049296e+08;
    }
  }
  else {
    dtref = 0.0e0;
    niterref = 0;
    verified = false;
  }

//---------------------------------------------------------------------
//    Compute the difference of solution values and the known reference values.
//---------------------------------------------------------------------
  //do m = 1, 5
  for (int m = 1; m <= 5; ++m) {

    xcrdif(m) = dabs((xcr(m)-xcrref(m))/xcrref(m));
    xcedif(m) = dabs((xce(m)-xceref(m))/xceref(m));

  } // enddo

//---------------------------------------------------------------------
//    Output the comparison of computed results to known cases.
//---------------------------------------------------------------------

  std::cout << std::scientific;
  std::cout << " Verification being performed for class " << CLASS << std::endl;
  std::cout << " accuracy setting for epsilon = " << epsilon << std::endl;
  if (dabs(dt-dtref) > epsilon) {
    verified = false;
    std::cout << " DT does not match the reference value of " << dtref << std::endl;
  } else if (no_time_steps != niterref) {
    std::cout << " NITER does not match the reference value of " << niterref << std::endl;
  }

  std::cout << "\n Comparison of RMS-norms of residual" << std::endl;
  for (int m = 1; m <= 5; ++m) {

    if (xcrdif(m) <= epsilon) {
      std::cout << "          " << m << " " << xcr(m) << " " << xcrref(m) << " " << xcrdif(m) << std::endl;
    } else {
      verified = false;
      std::cout << " FAILURE: " << m << " " << xcr(m) << " " << xcrref(m) << " " << xcrdif(m) << std::endl;
    }
  }

  std::cout <<  "\n Comparison of RMS-norms of solution error" << std::endl;

  for (int m = 1; m <= 5; ++m) {
    if (xcedif(m) <= epsilon) {
      std::cout << "          " << m << " " << xce(m) << " " << xceref(m) << " " << xcedif(m) << std::endl;
    } else {
      verified = false;
      std::cout << " FAILURE: " << m << " " << xce(m) << " " << xceref(m) << " " << xcedif(m) << std::endl;
    }
  }

  std::cout << " Verification " << std::string(verified ? " Successful" : " failed") << std::endl;

  /*
        write(*, 1990) class
 1990   format(' Verification being performed for class ', a)
        write (*,2000) epsilon
 2000   format(' accuracy setting for epsilon = ', E20.13)
        if (dabs(dt-dtref) .gt. epsilon) then
           verified = .false.
           write (*,1000) dtref
 1000      format(' DT does not match the reference value of ',
     >              E15.8)
        else if (no_time_steps .ne. niterref) then
           verified = .false.
           write (*,1002) niterref
 1002      format(' NITER does not match the reference value of ',
     >              I5)
        endif

        write (*,2001)

 2001   format(' Comparison of RMS-norms of residual')
        do m = 1, 5
           if (xcrdif(m) .le. epsilon) then
              write (*,2011) m,xcr(m),xcrref(m),xcrdif(m)
           else
              verified = .false.
              write (*,2010) m,xcr(m),xcrref(m),xcrdif(m)
           endif
        enddo

        write (*,2002)

 2002   format(' Comparison of RMS-norms of solution error')

        do m = 1, 5
           if (xcedif(m) .le. epsilon) then
              write (*,2011) m,xce(m),xceref(m),xcedif(m)
           else
              verified = .false.
              write (*,2010) m,xce(m),xceref(m),xcedif(m)
           endif
        enddo

 2010   format(' FAILURE: ', i2, E20.13, E20.13, E20.13)
 2011   format('          ', i2, E20.13, E20.13, E20.13)

        if (verified) then
           write(*, 2020)
 2020      format(' Verification Successful')
        else
           write(*, 2021)
 2021      format(' Verification failed')
        endif

  */

  return;
}

#endif // HAVE_VERIFY_H
