
#include <iostream>
#include <cstdlib>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "header.h"
#include "mpi_stuff.h"


int proc_num_zones;
ArrayT1<int, int> zone_proc_id;
ArrayT1<int, int> proc_zone_id;
ArrayT1<int, int> proc_zone_count;
ArrayT1<int, int> proc_num_threads;
ArrayT1<int, int> proc_group;

ArrayT1<double, int> proc_zone_size;


int myid, root;
MPI_Comm comm_setup;
int num_threads, mz_bload, max_threads;
bool active;

int npb_verbose, timeron;

ArrayT1<int, int> pcomm_group;
ArrayT1<int, int> qcomm_size, qstart2_west , qstart2_east,
                  qstart2_south, qstart2_north;

/*
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>c
c
      subroutine mpi_setup
c
c  Set up MPI stuff, including
c     - define the active set of processes
c     - set up new communicator
c
      include 'header.h'
c
      include 'mpi_stuff.h'
c
*/

void mpi_setup()
{

  int no_nodes, color;

 // initialize MPI parameters
  int provided;
  MPI_Init_thread(0, NULL, MPI_THREAD_MULTIPLE, &provided);
  assert(provided == MPI_THREAD_MULTIPLE);

  MPI_Comm_size(MPI_COMM_WORLD, &no_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  MPI_Continue_init(&request_completion_cb, &cont_req);


  zone_proc_id.allocate(max_zones);
  proc_zone_id.allocate(max_zones);
  proc_zone_count.allocate(num_procs);
  proc_num_threads.allocate(num_procs);
  proc_group.allocate(num_procs);

  proc_zone_size.allocate(num_procs);


  pcomm_group.allocate(num_procs2);
  qcomm_size.allocate(num_procs);
  qstart2_west.allocate(max_zones);
  qstart2_east.allocate(max_zones);
  qstart2_south.allocate(max_zones);
  qstart2_north.allocate(max_zones);

//---------------------------------------------------------------------
//     let node 0 be the root for the group (there is only one)
//---------------------------------------------------------------------
  root = 0;

  if (no_nodes < num_procs) {
      if (myid == root) {
        std::cout << " Requested MPI processes " << no_nodes
                  << " less than the compiled value " << num_procs << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
  } // endif

  if (myid >= num_procs) {
      active = false;
      color = 1;
  } else {
      active = true;
      color = 0;
  } // end if

  MPI_Comm_split(MPI_COMM_WORLD,color,myid,&comm_setup);
  if (!active) return;

  MPI_Comm_rank(comm_setup, &myid);
  if (no_nodes != num_procs) {
      if (myid == root) {
        std::cout << "Warning: Requested " << no_nodes << "MPI processes, "
                  << "but the compiled value is " << num_procs << ". "
                  <<  "The compiled value is used for benchmarking" << std::endl;
      }
  } // endif

  return;
}

/*
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>c
c
      subroutine env_setup(tot_threads)
c
c  Set up from environment variables
c
c ... common variables
      include 'header.h'
c
      include 'mpi_stuff.h'
c
c
*/

void env_setup(int tot_threads)
{
  int ios, curr_threads, ip, mp, group, ip1, ip2;
  ArrayT1<int, int> entry_counts{num_procs};

  std::string envstr;
  const char *tmp, *line;

  if (myid == root) {
    mp = 1;
#ifdef _OPENMP
    mp = omp_get_max_threads();
#endif

#if 0
    //c ... master sets up parameters
      envstr = getenv('OMP_NUM_THREADS');
      if (envstr .ne. ' ' .and. mp .gt. 0) then
         read(envstr,*,iostat=ios) num_threads
         if (ios.ne.0 .or. num_threads.lt.1) num_threads = 1
         if (mp .ne. num_threads) then
            write(*, 10) num_threads, mp
   10       format(' Warning: Requested ',i4,' threads per process,',
     &             ' but the active value is ',i4)
            num_threads = mp
         endif
      else
         num_threads = 1
      endif
#endif // 0
    num_threads = mp;

      tmp = getenv("NPB_MZ_BLOAD");
      envstr = tmp ? tmp : "";
      if (!envstr.empty()) {
        if (envstr == "on" || envstr == "ON") {
          mz_bload = 1;
        } else if (envstr[0] == 't' || envstr[0] == 'T') {
          mz_bload = 1;
        } else {
          mz_bload = atoi(envstr.c_str());
        } // endif
      } else {
        mz_bload = 1;
      } // endif

      tmp = getenv("NPB_MAX_THREADS");
      envstr = tmp ? tmp : "";
      max_threads = 0;
      if (mz_bload > 0 && !envstr.empty()) {
        max_threads = atoi(envstr.c_str());
        if (max_threads < 0) max_threads = 0;
        if (max_threads > 0 && max_threads < num_threads) {
        std::cout << "Error: max_threads " << max_threads
                  <<  " is less than num_threads " << num_threads
                  << "Please redefine the value for NPB_MAX_THREADS"
                  << "or OMP_NUM_THREADS";
        MPI_Abort(MPI_COMM_WORLD, 1);
        } // endif
      } // endif

      tmp = getenv("NPB_VERBOSE");
      envstr = tmp ? tmp : "";
      npb_verbose = 0;
      if (!envstr.empty()) {
         npb_verbose = atoi(envstr.c_str());
      }

      DO(ip, 1, num_procs, {
         proc_num_threads(ip) = num_threads;
         proc_group(ip) = 0;
      }); // end do

#if 0
      // TODO: load balancing from file is not supported atm
      open(2, file='loadbt-mz.data', status='old', iostat=ios)
      if (ios.eq.0) then
         write(*,*) 'Reading load factors from loadbt-mz.data'

         if (mz_bload .ge. 1) then
            mz_bload = -mz_bload
         endif

         do ip = 1, num_procs
            entry_counts(ip) = 0
         end do

         do while (.true.)
   25       read(2,'(a)',end=40,err=40) line
            if (line.eq.' ' .or. line(1:1).eq.'#') goto 25

            call decode_line(line, ip1, ip2, curr_threads, group, ios)
            if (ios .ne. 0) goto 40

            if (mz_bload .lt. 0 .and. group .gt. 0) then
               mz_bload = -mz_bload
            endif

            if (curr_threads .lt. 1) curr_threads = 1
            if (mp .le. 0) curr_threads = 1
            if (ip1.lt.0) ip1 = 0
            if (ip2.ge.num_procs) ip2 = num_procs - 1

            do ip = ip1+1, ip2+1
               proc_num_threads(ip) = curr_threads
               proc_group(ip) = group
               entry_counts(ip) = entry_counts(ip) + 1
            end do
         end do
   40    close(2)

         do ip = 1, num_procs
            if (entry_counts(ip) .eq. 0) then
               write(*,*) '*** Error: Missing entry for proc ',ip-1
               call mpi_abort(MPI_COMM_WORLD, 1, ierror)
               stop
            else if (entry_counts(ip) .gt. 1) then
               write(*,*) '*** Warning: Multiple entries for proc ',
     &                    ip-1, ', only the last one used'
            endif
         end do

         ip1 = 1
      else
#endif // 0
      std::cout << "Use the default load factors with threads" << std::endl;
      ip1 = 0;
#if 0
      endif
#endif // 0

      if (ip1 > 0 || npb_verbose > 0) {
        ip1 = 0;
        DO(ip, 1, num_procs, {
          if (ip == 1 ||
              proc_num_threads(ip) != curr_threads ||
              proc_group(ip) != group) {

            ip2 = ip-2;
            if (ip2 > ip1+1) std::cout << " ..." << std::endl; // write(*,*) '    ...'
            if (ip2 > ip1) {
              std::cout << "  proc " << ip2 << "  num_threads = " << curr_threads
                        << "  group = " << group << std::endl;
            // write(*,30) ip2, curr_threads, group
            }

            curr_threads = proc_num_threads(ip);
            group = proc_group(ip);

            ip1 = ip - 1;
            std::cout << "  proc " << ip2 << "  num_threads = " << curr_threads
                      << "  group = " << group << std::endl;
              //write(*,30) ip1, curr_threads, group

          } else if (ip == num_procs) {
            ip2 = ip-1;
            // if (ip2 .gt. ip1+1) write(*,*) '    ...'
            if (ip2 > ip1+1) std::cout << " ..." << std::endl;
            std::cout << "  proc " << ip2 << "  num_threads = " << curr_threads
                      << "  group = " << group << std::endl;
          } // endif
        }); // end do
      } // endif

// c ... broadcast parameters to all processes
  }
#if 0
   80 call mpi_bcast(num_threads, 1, mpi_integer, root,
     &               comm_setup, ierror)
      call mpi_bcast(mz_bload, 1, mpi_integer, root,
     &               comm_setup, ierror)
      call mpi_bcast(max_threads, 1, mpi_integer, root,
     &               comm_setup, ierror)
      call mpi_bcast(proc_num_threads, num_procs, mpi_integer, root,
     &               comm_setup, ierror)
      call mpi_bcast(proc_group, num_procs, mpi_integer, root,
     &               comm_setup, ierror)
#endif // 0
  MPI_Bcast(&num_threads, 1, MPI_INT, root, comm_setup);
  MPI_Bcast(&mz_bload,    1, MPI_INT, root, comm_setup);
  MPI_Bcast(&max_threads, 1, MPI_INT, root, comm_setup);
  MPI_Bcast(proc_num_threads.begin(), num_procs, MPI_INT, root, comm_setup);
  MPI_Bcast(proc_group.begin(),  num_procs, MPI_INT, root, comm_setup);
  MPI_Bcast(&npb_verbose, 1, MPI_INT, root, comm_setup);

#if 0
      tot_threads = 0
      do ip = 1, num_procs
         tot_threads = tot_threads + proc_num_threads(ip)
      end do
      if (myid .eq. root) then
         if (mp .gt. 0)
     &      write(*, 1004) tot_threads, dble(tot_threads)/num_procs
      endif
 1004 format(' Total number of threads: ', i6,
     &       '  (', f5.1, ' threads/process)')
c
#endif // 0
  return;
}

#if 0
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>c
c
      subroutine decode_line(line, ip1, ip2, curr_threads, group, ios)
      implicit none
c
c  decode a line from the load data file
c  format:  ip1[:ip2] curr_threads group
c
      character line*(*)
      integer ip1, ip2, curr_threads, group, ios
c
      integer is, n
c
      ios = -1
c
      n  = len(line)
      is = 1
      do while (is.le.n .and. line(is:is).ne.':')
         if (line(is:is).eq.'!') n = is
         is = is + 1
      end do
c
      if (is .gt. n) then
         read(line,*,err=90,end=90) ip1, curr_threads, group
         ip2 = ip1
      else if (is.eq.1 .or. is.eq.n) then
         go to 90
      else
         read(line(:is-1),*,err=90,end=90) ip1
         read(line(is+1:),*,err=90,end=90) ip2, curr_threads, group
      endif
c
      if (ip2 .lt. ip1) then
         is  = ip2
         ip2 = ip1
         ip1 = is
      endif
      ios = 0
c
   90 return
      end

#endif


/*
c
c  Calculate the communication index of a zone within a processor group
c
*/
static int get_comm_index(int zone, int iproc)
{
  int comm_index;
  int izone, jzone;
  jzone  = (zone - 1)/x_zones + 1;
  izone  = (zone - 1)%x_zones + 1;

  comm_index = 0;
  if (zone_proc_id(iz_west(zone)) == iproc) {
    comm_index = comm_index + y_size(jzone);
  }
  if (zone_proc_id(iz_east(zone)) == iproc) {
    comm_index = comm_index + y_size(jzone);
  }
  if (zone_proc_id(iz_south(zone)) == iproc) {
    comm_index = comm_index + x_size(izone);
  }
  if (zone_proc_id(iz_north(zone)) == iproc) {
    comm_index = comm_index + x_size(izone);
  }

  return comm_index;
}

/*
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>c
c
      subroutine get_comm_index(zone, iproc, comm_index)
c
      include 'header.h'
c
      include 'mpi_stuff.h'
c
c  Calculate the communication index of a zone within a processor group
c
      integer zone, iproc, comm_index
c
c     local variables
      integer izone, jzone
c
      jzone  = (zone - 1)/x_zones + 1
      izone  = mod(zone - 1, x_zones) + 1
c
      comm_index = 0
      if (zone_proc_id(iz_west(zone)) .eq. iproc)
     $   comm_index = comm_index + y_size(jzone)
      if (zone_proc_id(iz_east(zone)) .eq. iproc)
     $   comm_index = comm_index + y_size(jzone)
      if (zone_proc_id(iz_south(zone)) .eq. iproc)
     $   comm_index = comm_index + x_size(izone)
      if (zone_proc_id(iz_north(zone)) .eq. iproc)
     $   comm_index = comm_index + x_size(izone)
c
      return
      end
*/


void map_zones(
  int num_zones,
  ArrayT1<size_type, size_type>& nx,
  ArrayT1<size_type, size_type>& ny,
  ArrayT1<size_type, size_type>& nz,
  int tot_threads)
{

/*
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>c
c
      subroutine map_zones(num_zones, nx, ny, nz, tot_threads)
c
c  Perform zone-process mapping for load balance
c
      include 'header.h'
c
      integer num_zones, nx(*), ny(*), nz(*), tot_threads
c
      include 'mpi_stuff.h'
c
c     local variables
*/
  ArrayT1<int, int> z_order(max_zones);
  int zone, iz, z2, mz, np, ip, zone_comm, comm_index;
  int imx, imn, inc, icur_size;
  double tot_size, cur_size, max_size, ave_size;
  ArrayT1<double, int> zone_size(max_zones);
  double diff_ratio, tot_group_size;

  int group, ipg, tot_group_threads;
  ArrayT1<int, int> proc_group_flag(num_procs);

  // sort the zones in decending order
  tot_size = 0.0;
  //do iz = 1, num_zones
  for (int iz = 1; iz <= num_zones; ++iz) {
    zone_size(iz) = 1.00*nx(iz)*ny(iz)*nz(iz);
    z_order(iz) = iz;
    tot_size = tot_size + zone_size(iz);
  } // end do

  //do iz = 1, num_zones-1
  for (int iz = 1; iz < num_zones; ++iz) {
      cur_size = zone_size(z_order(iz));
      mz = iz;
      //do z2 = iz+1, num_zones
      for (int z2 = iz+1; z2 <= num_zones; ++z2) {
        if (cur_size < zone_size(z_order(z2))) {
          cur_size = zone_size(z_order(z2));
          mz = z2;
        } // endif
      } // end do
      if (mz != iz) {
        z2 = z_order(iz);
        z_order(iz) = z_order(mz);
        z_order(mz) = z2;
      } // endif
  } // end do

  if (npb_verbose > 1 && myid == root) {
      std::cout << "Sorted zones: " << "seq. zone    nx    ny    nz    size" << std::endl;
      //do iz = 1, num_zones
      for (int iz = 1; iz <= num_zones; ++iz) {
        z2 = z_order(iz);
        std::cout << iz << ": " << z2 << " " << nx(z2) << " " << ny(z2)
                  << " " << nz(z2) << " " << zone_size(z2) << std::endl;
      } // end do
  } // endif

  // use a bin-packing scheme to balance the load among processes
  //do ip = 1, num_procs
  for (int ip = 1; ip <= num_procs; ++ip) {
    proc_zone_count(ip) = 0;
    proc_zone_size(ip) = 0.0;
  } // end do

  //do iz = 1, num_zones
  for (int iz = 1; iz <= num_zones; ++iz) {
    zone_proc_id(iz) = -1;
  } // end do

  iz = 1;
  //do while (iz .le. num_zones)
  while (iz <= num_zones) {
    // the current most empty processor
    np = 1;
    cur_size = proc_zone_size(1);
    // do ip = 2, num_procs
    for (int ip = 2; ip <= num_procs; ++ip) {
      if (cur_size > proc_zone_size(ip)) {
        np = ip;
        cur_size = proc_zone_size(ip);
      } // endif
    } // end do
    ip = np - 1;

    // get a zone that has the largest communication index with
    // the current group and does not worsen the computation balance
    mz = z_order(iz);
    if (iz < num_zones) {
      zone_comm = get_comm_index(mz, ip);
      // do z2 = iz+1, num_zones
      for (int z2 = iz+1; z2 <= num_zones; ++z2) {
        zone = z_order(z2);

        diff_ratio = (zone_size(z_order(iz)) -
                    zone_size(zone)) / zone_size(z_order(iz));
        if (diff_ratio > 0.05) break;

        if (zone_proc_id(zone) < 0) {
          comm_index = get_comm_index(zone, ip);
          if (comm_index > zone_comm) {
            mz = zone;
            zone_comm = comm_index;
          } // endif
        } // endif
      } // end do
    } // endif

    // assign the zone to the current processor group
    zone_proc_id(mz) = ip;
    proc_zone_size(np) = proc_zone_size(np) + zone_size(mz);
    proc_zone_count(np) = proc_zone_count(np) + 1;
    // skip the previously assigned zones
    while (iz <= num_zones) {
      if (zone_proc_id(z_order(iz)) < 0) break;
      iz = iz + 1;
    } // end do
  } // end do


  // move threads around if needed
  mz = 1;
  if (tot_threads == num_procs || mz_bload < 1) mz = 0;

      if (mz != 0) {

        DO(ipg, 1, num_procs, {
          proc_group_flag(ipg) = 0;
        }); // end do

        ipg = 1;

        // balance load within a processor group
        do {
  // 200    do while (ipg .le. num_procs)
          while (ipg <= num_procs) {
            if (proc_group_flag(ipg) == 0) break;
            ipg = ipg + 1;
          } // end do

          if (ipg > num_procs) break;


          group = proc_group(ipg);
          tot_group_size = 0.0;
          tot_group_threads = 0;
          DO(ip, ipg, num_procs, {
            if (proc_group(ip) == group) {
              proc_group_flag(ip) = 1;
              tot_group_size = tot_group_size + proc_zone_size(ip);
              tot_group_threads = tot_group_threads +
                                  proc_num_threads(ip);
            } // endif
          }); //end do

          ave_size = tot_group_size / tot_group_threads;

          // distribute size evenly among threads
          icur_size = 0;
          for(ip = 1; ip <= num_procs; ip++) {
            if (proc_group(ip) != group) continue;
            proc_num_threads(ip) = proc_zone_size(ip) / ave_size;
            if (proc_num_threads(ip) < 1) {
              proc_num_threads(ip) = 1;
            }
            if (max_threads > 0 &&
                proc_num_threads(ip) > max_threads) {
              proc_num_threads(ip) = max_threads;
            }
            icur_size = icur_size + proc_num_threads(ip);
          } // end do
          mz = tot_group_threads - icur_size;

          // take care of any remainers
          inc = 1;
          if (mz < 0) inc = -1;
          //do while (mz .ne. 0)
          while (mz != 0) {
            max_size = 0.0;
            imx = 0;
            DO(ip, 1, num_procs, {
              if (proc_group(ip) == group) {
                if (mz > 0) {
                  cur_size = proc_zone_size(ip) / proc_num_threads(ip);
                  if (cur_size > max_size && (max_threads <= 0 || proc_num_threads(ip) < max_threads)) {
                    max_size = cur_size;
                    imx = ip;
                  } // endif
                } else if (proc_num_threads(ip) > 1) {
                  cur_size = proc_zone_size(ip) /
                              (proc_num_threads(ip)-1);
                  if (max_size == 0 || cur_size < max_size) {
                    max_size = cur_size;
                    imx = ip;
                  } // endif
                } // endif
  //230          continue
              }
            }); // end do
            proc_num_threads(imx) = proc_num_threads(imx) + inc;
            mz = mz - inc;
          } // end do

         //goto 200
        } while (true);
      } // endif

  // print the mapping
  // 300 if (npb_verbose .gt. 0 .and. myid .eq. root) then
  if (npb_verbose > 0 && myid == root) {
    std::cout << "Zone-process mapping: proc nzones  zone_size nthreads size_per_thread" << std::endl;
    // do ip = 1, num_procs
    for (int ip = 1; ip <= num_procs; ++ip) {
      std::cout << ip-1 << " " << proc_zone_count(ip) << " "
                << proc_zone_size(ip) << " " << proc_num_threads(ip) << " "
                << proc_zone_size(ip)/proc_num_threads(ip) << std::endl;
      //do iz = 1, num_zones
      for (int iz = 1; iz <= num_zones; ++iz) {
        if (zone_proc_id(iz) == ip-1) {
          std::cout << iz << " " << zone_size(iz) << std::endl;
        } // endif
      } // end do
    } // end do
  } // endif

  if (myid == root) {
    imx = 1;
    max_size = proc_zone_size(1)/proc_num_threads(1);
    imn = imx;
    ave_size = max_size;
    DO(ip, 2, num_procs, {
      cur_size = proc_zone_size(ip)/proc_num_threads(ip);
      if (cur_size > max_size) {
        imx = ip;
        max_size = cur_size;
      } // endif
      if (cur_size < ave_size) {
        imn = ip;
        ave_size = cur_size;
      } // endif
    }); // end do

    if (npb_verbose > 0) {
      std::cout << "Max: proc=" << imx-1 << " nzones=" << proc_zone_count(imx)
                << " size=" << proc_zone_size(imx) << " nthreads=" << proc_num_threads(imx)
                << std::endl;
      std::cout << "Min: proc=" << imn-1 << " nzones=" << proc_zone_count(imn)
                << " size=" << proc_zone_size(imn) << " nthreads=" << proc_num_threads(imn)
                << std::endl;
    } // endif

    std::cout << "Calculated speedup = " << tot_size / max_size << std::endl;
  } // endif

  // reorganize list of zones for this process
  zone = 0;
  DO(iz, 1, num_zones, {
    if (zone_proc_id(iz) == myid) {
      zone = zone + 1;
      proc_zone_id(zone) = iz;
    } // endif
  }); // end do
  proc_num_zones = zone;
  if (zone != proc_zone_count(myid+1)) {
    std::cout << "Warning: " << myid << ": mis-matched zone counts - "
              << zone << " " << proc_zone_count(myid+1);
  } // endif

  // set number of threads for this process
  group = proc_group(myid+1);
  np = 0;
  DO(ip, 1, num_procs, {
    if (proc_group(ip) == group) {
      proc_group(ip) = np;
      np = np + 1;
      proc_num_threads(np) = proc_num_threads(ip);
    } // endif
  }); // end do
  ipg = proc_group(myid+1);
  if (npb_verbose > 1) {
    std::cout << "myid " << myid << " group "<< group << " group_size " << np <<
        " group_pid "<< ipg << " threads "<< proc_num_threads(ipg+1);
  } // endif
#if 0
c$    call omp_set_num_threads(proc_num_threads(ipg+1))
c
c ... pin-to-node within one process group
c      call smp_pinit_thread(np, ipg, proc_num_threads)
c
#endif // 0
  return;
}
