/*
!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         !
!                                                                         !
!             M P I    M U L T I - Z O N E    V E R S I O N               !
!                                                                         !
!                           B T - M Z - M P I                             !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is an MPI+OpenMP version of the NPB BT code.          !
!    Refer to NAS Technical Reports 95-020 and 99-011 for details.        !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.  We           !
!    request, however, that all derived work reference the NAS            !
!    Parallel Benchmarks 3.3. This software is provided "as is"           !
!    without express or implied warranty.                                 !
!                                                                         !
!    Information on NPB 3.3, including the technical report, the          !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!           http://www.nas.nasa.gov/Software/NPB/                         !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          NAS Parallel Benchmarks Group                                  !
!          NASA Ames Research Center                                      !
!          Mail Stop: T27A-1                                              !
!          Moffett Field, CA   94035-1000                                 !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!-------------------------------------------------------------------------!

c---------------------------------------------------------------------
c
c Authors: R. Van der Wijngaart
c          T. Harris
c          M. Yarrow
c          H. Jin
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       program BT
c---------------------------------------------------------------------
*/

#include <iostream>

#include "NArray.h"
#include "mpi_stuff.h"
#include "header.h"

#include "zone_setup.h"

#include "adi.h"
#include "verify.h"

#include "../common/print_results.h"

int main(int argc, char *argv[])
{

  constexpr const int num_zones = x_zones*y_zones;

  ArrayT1<size_type, size_type> nx{num_zones},
                                nxmax{num_zones},
                                ny{num_zones},
                                nz{num_zones};

/*
c---------------------------------------------------------------------
c   Define all field arrays as one-dimenional arrays, to be reshaped
c---------------------------------------------------------------------
       double precision
     >   u       (proc_max_size5),
     >   us      (proc_max_size ),
     >   vs      (proc_max_size ),
     >   ws      (proc_max_size ),
     >   qs      (proc_max_size ),
     >   rho_i   (proc_max_size ),
     >   square  (proc_max_size ),
     >   rhs     (proc_max_size5),
     >   forcing (proc_max_size5),
     >   qbc_ou  (proc_max_bcsize),
     >   qbc_in  (proc_max_bcsize)

       common /fields/ u, us, vs, ws, qs, rho_i, square,
     >                 rhs, forcing, qbc_ou, qbc_in
*/
  ArrayT1<double, size_type>
        u       {proc_max_size5},
        us      {proc_max_size },
        vs      {proc_max_size },
        ws      {proc_max_size },
        qs      {proc_max_size },
        rho_i   {proc_max_size },
        square  {proc_max_size },
        rhs     {proc_max_size5},
        forcing {proc_max_size5},
        qbc_ou  {proc_max_bcsize},
        qbc_in  {proc_max_bcsize};

  int niter, step, fstatus, zone, iz, tot_threads, itimer;

  double tmax, t;
  StaticArrayT1<double, int, t_last, 1> trecs;
  //character        t_names(t_last)*8
  StaticArrayT1<const char*, int, t_last, 1> t_names;
  bool verified;


  mpi_setup();
  if (active) {

/*
c---------------------------------------------------------------------
c      Root node reads input file (if it exists) else takes
c      defaults from parameters
c---------------------------------------------------------------------
*/
    if (myid == root) {

#if 0
         write(*, 1000)
         open (unit=2,file='inputbt-mz.data',status='old',
     >         iostat=fstatus)

         timeron = .false.
         if (fstatus .eq. 0) then
           write(*,*) 'Reading from input file inputbt-mz.data'
           read (2,*) niter
           read (2,*) dt
           read (2,*) itimer
           close(2)

           if (niter .eq. 0)  niter = niter_default
           if (dt .eq. 0.d0)  dt    = dt_default
           if (itimer .gt. 0) timeron = .true.

         else
#endif // 0
           // TODO: reading from data file not yet supported!
           niter = niter_default;
           dt    = dt_default;
#if 0
         endif
#endif // 0

      std::cout << " NAS Parallel Benchmarks (NPB3.3-MZ-MPI)"
                << " - BT-MZ MPI+OpenMP Benchmark" << std::endl;
      std::cout << " Number of zones: " << x_zones << " x " << y_zones << std::endl;
      std::cout << " Iterations: " << niter << "    dt: " << dt << std::endl;
      std::cout << " Number of active processes: " << num_procs << std::endl;
    } // endif
    MPI_Bcast(&niter,   1, MPI_INT,    root, comm_setup);
    MPI_Bcast(&dt,      1, MPI_DOUBLE, root, comm_setup);
    MPI_Bcast(&timeron, 1, MPI_INT,    root, comm_setup);

    if (timeron) {
      t_names(t_total) = "total";
      t_names(t_rhsx) = "rhsx";
      t_names(t_rhsy) = "rhsy";
      t_names(t_rhsz) = "rhsz";
      t_names(t_rhs) = "rhs";
      t_names(t_xsolve) = "xsolve";
      t_names(t_ysolve) = "ysolve";
      t_names(t_zsolve) = "zsolve";
      t_names(t_rdis1) = "qbc_copy";
      t_names(t_rdis2) = "qbc_comm";
      t_names(t_add) = "add";
    } // endif

    env_setup(tot_threads);

    zone_setup(nx, nxmax, ny, nz);

    map_zones(num_zones, nx, ny, nz, tot_threads);

    zone_starts(num_zones, nx, nxmax, ny, nz);

    set_constants();

    // do iz = 1, proc_num_zones
    for (int iz = 1; iz <= proc_num_zones; ++iz) {
      int zone = proc_zone_id(iz);

      initialize(&u(start5(iz)), nx(zone), nxmax(zone), ny(zone), nz(zone));
      exact_rhs(&forcing(start5(iz)), nx(zone), nxmax(zone), ny(zone), nz(zone));

      
      print_zone(&u(start5(iz)), nx(zone), nxmax(zone), ny(zone), nz(zone), zone, -1);

    } // end do

    doloop(1, t_last, [](auto i){
      timer_clear(i);
    }); // end do

/*
c---------------------------------------------------------------------
c      do one time step to touch all code, and reinitialize
c---------------------------------------------------------------------
*/
#pragma omp parallel
{
#pragma omp master
{
    exch_qbc(u, qbc_ou, qbc_in, nx, nxmax, ny, nz, npb_verbose);

    // wait for all tasks to complete
#pragma omp taskwait

    doloop(1, proc_num_zones, [&](auto iz){
      int zone = proc_zone_id(iz);
      adi(&rho_i(start1(iz)), &us(start1(iz)),
          &vs(start1(iz)), &ws(start1(iz)),
          &qs(start1(iz)), &square(start1(iz)),
          &rhs(start5(iz)), &forcing(start5(iz)),
          &u(start5(iz)),
          nx(zone), nxmax(zone), ny(zone), nz(zone));
    });// end do
} // pragma omp master
} // pragma omp parallel

    doloop(1, proc_num_zones, [&](auto iz) {
      int zone = proc_zone_id(iz);
      initialize(&u(start5(iz)),
                  nx(zone), nxmax(zone), ny(zone), nz(zone));
    }); // end do


    doloop(1, t_last, [](auto i){
      timer_clear(i);
    }); // end do

    MPI_Barrier(comm_setup);

    timer_start(1);

/*
c---------------------------------------------------------------------
c      start the benchmark time step loop
c---------------------------------------------------------------------
*/

#pragma omp parallel default(shared)
{
#pragma omp master
{
    //doloop(1, niter, [&](auto step) {
    for (int step = 1; step <= niter; ++step) {
      if (step % 20 == 0 || step == 1) {
        if (myid == root) {
          std::cout << "Time step " << step << std::endl;
        }
      } // endif

      exch_qbc(u, qbc_ou, qbc_in, nx, nxmax, ny, nz, 0);

// TODO: is this needed?
//#pragma omp taskwait

      //doloop(1, proc_num_zones, [&](auto iz){
      for (int iz = 1; iz <= proc_num_zones; ++iz) {
        auto *u_ptr = &u(start5(iz));
#pragma omp task depend(out:u_ptr[0])
{
        int zone = proc_zone_id(iz);

        print_zone(&u(start5(iz)), nx(zone), nxmax(zone), ny(zone), nz(zone), zone, step);
        adi(&rho_i(start1(iz)), &us(start1(iz)),
            &vs(start1(iz)), &ws(start1(iz)),
            &qs(start1(iz)), &square(start1(iz)),
            &rhs(start5(iz)), &forcing(start5(iz)),
            &u(start5(iz)),
            nx(zone), nxmax(zone), ny(zone), nz(zone));
#pragma omp taskwait
}
      }//); // end do

    } //); // end do
} // pragma omp master
} // pragma omp parallel

    timer_stop(1);
    tmax = timer_read(1);

/*
c---------------------------------------------------------------------
c      perform verification and print results
c---------------------------------------------------------------------
*/
    verify(niter, verified, num_zones, rho_i, us, vs, ws,
          qs, square, rhs, forcing, u, nx, nxmax, ny, nz);

    t = tmax;
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, root, comm_setup);

    if (myid == root) {

      double mflops = 0.0;
      if( tmax != 0.0 ) {
        doloop(1, num_zones, [&](auto zone){
          size_type n3 = nx(zone)*ny(zone)*nz(zone);
          double navg = (nx(zone) + ny(zone) + nz(zone))/3.0;
          double nsur = (nx(zone)*ny(zone) + nx(zone)*nz(zone) +
                  ny(zone)*nz(zone))/3.0;
          mflops = mflops + 1.0e-6*niter *
                    (3478.8 * n3 - 17655.7 * nsur + 28023.7 * navg)
                    / tmax;
        }); // end do
      } // endif

      print_results("BT-MZ", CLASS, gx_size, gy_size, gz_size,
                    niter, tmax, mflops,
                    "          floating point",
                    verified, npbversion, compiletime, cs1, cs2,
                    cs3, cs4, cs5, cs6);

    }
/*
c---------------------------------------------------------------------
c      More timers
c---------------------------------------------------------------------
*/
    if (timeron) {

       doloop(1, t_last, [&](auto i){
          trecs(i) = timer_read(i);
       }); // end do

       if (myid > 0) {
#if 0
        // TODO: what are these send/recv pairs for??
        MPI_Recv(i, 1, MPI_INTEGER, 0, 1000,
                       comm_setup, statuses);
          call mpi_send(trecs, t_last, dp_type, 0, 1001,
                       comm_setup, ierror)
#endif // 0
       } else {

        int ip = 0;
        if (tmax == 0.0) tmax = 1.0;
        std::cout << " Myid =" << ip << "   num_threads =" << proc_num_threads(ip+1)
              << "  SECTION   Time (secs)" << std::endl;
        doloop(1, t_last, [&](auto i){
          std::cout << t_names(i) << ':' << trecs(i) << "  (" << trecs(i)*100./tmax << "%)" << std::endl;
            if (i==t_rhs) {
              t = trecs(t_rhsx) + trecs(t_rhsy) + trecs(t_rhsz);
              std::cout << "    --> total " << "sub-rhs" << ':' << t << "  (" << t*100./tmax << "%)" << std::endl;
              t = trecs(t_rhs) - t;
              std::cout << "    --> total " << "rest-rhs" << ':' << t << "  (" << t*100./tmax << "%)" << std::endl;
            } else if (i == t_rdis2) {
              t = trecs(t_rdis1) + trecs(t_rdis2);
              std::cout << "    --> total " << "exch_qbc" << ':' << t << "  (" << t*100./tmax << "%)" << std::endl;
            } // endif
        }); // end do
#if 0
        // TODO: what are these send/recv pairs for??
        ip = ip + 1;
        if (ip < num_procs) {
            call mpi_send(myid, 1, MPI_INTEGER, ip, 1000,
                        comm_setup, ierror);
            call mpi_recv(trecs, t_last, dp_type, ip, 1001,
                        comm_setup, statuses, ierror);
            write(*,*)
            goto 910
        }
#endif // 0
      }

    }
 //999   continue
  }
  //mpi_barrier(MPI_COMM_WORLD, ierror)
  MPI_Finalize();

  return 0;
}

