
#include<iostream>
#include<assert.h>

#include "NArray.h"
#include "mpi_stuff.h"
#include "header.h"


/* 3D column-major Array, 0-index-based */
using QBCArrayT = NArrayView<double, 3, NArrayShape<3, size_type, NARRAY_MEMORDER_COL, 1, 1, 1>>;
/* 4D column-major Array, 0-index-based except on the fastest dimension */
using UArrayT = NArrayView<double, 4, NArrayShape<4, size_type, NARRAY_MEMORDER_COL, 1, 0, 0, 0>>;

enum {
  DIR_IN  = 0,
  DIR_OUT = 1
};

enum{
  ZONE_WEST = 0,
  ZONE_EAST,
  ZONE_SOUTH,
  ZONE_NORTH
};

static void
copy_y_face(
  double* u_ptr,
  double* qbc_ptr,
  int nx,
  int nxmax,
  int ny,
  int nz,
  int jloc,
  int dir);

static void
copy_x_face(
  double* u_ptr,
  double* qbc_ptr,
  int nx,
  int nxmax,
  int ny,
  int nz,
  int iloc,
  int dir);

void exch_qbc(
  ArrayT1<double, size_type>& u,
  ArrayT1<double, size_type>& qbc_ou,
  ArrayT1<double, size_type>& qbc_in,
  ArrayT1<size_type, size_type>& nx,
  ArrayT1<size_type, size_type>& nxmax,
  ArrayT1<size_type, size_type>& ny,
  ArrayT1<size_type, size_type>& nz,
  int iprn_msg)
{

/*
       subroutine exch_qbc(u, qbc_ou, qbc_in, nx, nxmax, ny, nz,
     $                     iprn_msg)

       include 'header.h'
       include 'mpi_stuff.h'

       integer   nx(*), nxmax(*), ny(*), nz(*), iprn_msg
       double precision u(*), qbc_ou(*), qbc_in(*)
*/
/*
  int c_size, m_size, num_msgs;
  constexpr const int MSG_TAG = 10000;
  StaticArrayT1<MPI_Request, int, MAX_REQS, 1> requests;
*/

  // wait for all processes to complete previous reads
  dash::barrier();

  // copy data to qbc buffer
  if (timeron) timer_start(t_rdis1);
  // do iz = 1, proc_num_zones
  for (int iz = 1; iz <= proc_num_zones; ++iz) {
      int zone_no = proc_zone_id(iz);
      int nnx    = nx(zone_no);
      int nnxmax = nxmax(zone_no);
      int nny    = ny(zone_no);
      int nnz    = nz(zone_no);
      auto lzone = iz-1;

      copy_x_face(&u(start5(iz)),
                  &halo_array.local(0, lzone, ZONE_WEST, 0),
                  nnx, nnxmax, nny, nnz, 1, DIR_OUT);

      copy_x_face(&u(start5(iz)),
                  &halo_array.local(0, lzone, ZONE_EAST, 0),
                  nnx, nnxmax, nny, nnz, nnx-2, DIR_OUT);

      copy_y_face(&u(start5(iz)),
                  &halo_array.local(0, lzone, ZONE_SOUTH, 0),
                  nnx, nnxmax, nny, nnz, 1, DIR_OUT);

      copy_y_face(&u(start5(iz)),
                  &halo_array.local(0, lzone, ZONE_NORTH, 0),
                  nnx, nnxmax, nny, nnz, nny-2, DIR_OUT);

  } // end do
  if (timeron) timer_stop(t_rdis1);


  // exchange qbc buffers
  if (timeron) timer_start(t_rdis2);


  // wait for all processes to complete filling their buffers
  dash::barrier();

  auto& pattern = halo_array.pattern();

  for (int iz = 1; iz <= proc_num_zones; ++iz) {
      int zone_no = proc_zone_id(iz);

      int ip_west   = zone_proc_id(iz_west(zone_no));
      int ip_east   = zone_proc_id(iz_east(zone_no));
      int ip_south  = zone_proc_id(iz_south(zone_no));
      int ip_north  = zone_proc_id(iz_north(zone_no));

      size_type x_face_size = (ny(zone_no)-2)*(nz(zone_no)-2)*5;
      size_type y_face_size = (nx(zone_no)-2)*(nz(zone_no)-2)*5;

      if (ip_west != myid) {
        auto gzone = iz_west(zone_no);
        auto lzone = static_cast<long>(proc_gzone_id(gzone))-1;
        auto offset = pattern.global_at({ip_west, lzone, ZONE_EAST, 0});
        auto gptr = halo_array.begin() + offset;
        //std::cout << "Zone " << zone_no << " WEST global " << gzone << " local "
        //          << lzone << " at unit " << ip_west << " offset " << offset << " gptr " << gptr << std::endl;
        dash::copy_async(gptr, gptr + x_face_size, &qbc_in(qstart2_west(zone_no)));
      }

      if (ip_east != myid) {
        auto gzone = iz_east(zone_no);
        auto lzone = static_cast<long>(proc_gzone_id(gzone))-1;
        auto offset = pattern.global_at({ip_east, lzone, ZONE_WEST, 0});
        auto gptr = halo_array.begin() + offset;
        //std::cout << "Zone " << zone_no << " EAST global " << gzone << " local "
        //          << lzone << " offset " << offset << " gptr " << gptr << std::endl;
        dash::copy_async(gptr, gptr + x_face_size, &qbc_in(qstart2_east(zone_no)));
      }

      if (ip_north != myid) {
        auto gzone = iz_north(zone_no);
        auto lzone = static_cast<long>(proc_gzone_id(gzone))-1;
        auto offset = pattern.global_at({ip_north, lzone, ZONE_SOUTH, 0});
        auto gptr = halo_array.begin() + offset;
        //std::cout << "Zone " << zone_no << " NORTH global " << gzone << " local "
        //          << lzone << " offset " << offset << " gptr " << gptr << std::endl;
        dash::copy_async(gptr, gptr + y_face_size, &qbc_in(qstart2_north(zone_no)));
      }

      if (ip_south != myid) {
        auto gzone = iz_south(zone_no);
        auto lzone = static_cast<long>(proc_gzone_id(gzone))-1;
        auto offset = pattern.global_at({ip_south, lzone, ZONE_NORTH, 0});
        auto gptr = halo_array.begin() + offset;
        //std::cout << "Zone " << zone_no << " SOUTH global " << gzone << " local "
        //          << lzone << " offset " << offset << " gptr " << gptr << std::endl;
        dash::copy_async(gptr, gptr + y_face_size, &qbc_in(qstart2_south(zone_no)));
      }
  }

  // wait for all transfers to complete
  halo_array.flush();

#if 0
  for (int ig = 1; ig <= num_procs; ++ig) {
    int ip = pcomm_group(ig);

    if (ip == 0) {
        c_size = qcomm_size(ip+1);
    } else {
        c_size = qcomm_size(ip+1) - qcomm_size(ip);
    } // endif

    int nr = 0;
    if (c_size > 0) {
        num_msgs = c_size / MSG_SIZE;
        if (num_msgs == 0) num_msgs = 1;
        m_size = (c_size + num_msgs - 1)/ num_msgs;

        if (iprn_msg > 0) {
          std::cout << "myid,msgs" << myid << " " << ip << " "
                    << num_msgs << " " << m_size << std::endl;
        }
        int qoffset = qcomm_size(ip+1) - c_size + 1;
        int tag = MSG_TAG;
        for (int n = 1; n <= num_msgs; ++n) {

          if (nr >= MAX_REQS) {
            MPI_Waitall(nr, requests.begin(), MPI_STATUSES_IGNORE);
            nr = 0;
            tag = MSG_TAG;
          } // endif

          if (qoffset+m_size-1 > qcomm_size(ip+1)) {
            m_size = qcomm_size(ip+1) - qoffset + 1;
          } // endif

          MPI_Isend(&qbc_ou(qoffset), m_size,
                    MPI_DOUBLE, ip, tag+myid,
                    comm_setup, &requests(nr+1));

          MPI_Irecv(&qbc_in(qoffset), m_size,
                    MPI_DOUBLE, ip, tag+ip,
                    comm_setup, &requests(nr+2));

          nr = nr + 2;
          qoffset = qoffset + m_size;
          tag = tag + num_procs;
        } // end do
    } else if (c_size < 0) {
        std::cout << "error: integer overflow" << myid << " " << ip << " " << c_size << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    } // endif

    if (nr > 0) {
      MPI_Waitall(nr, requests.begin(), MPI_STATUSES_IGNORE);
    } // endif

  } // enddo
#endif // 0

  if (timeron) timer_stop(t_rdis2);

  // copy data from qbc buffer
  if (timeron) timer_start(t_rdis1);
  for (int iz = 1; iz <= proc_num_zones; ++iz) {
      int zone_no = proc_zone_id(iz);
      int nnx    = nx(zone_no);
      int nnxmax = nxmax(zone_no);
      int nny    = ny(zone_no);
      int nnz    = nz(zone_no);

      int ip_west   = zone_proc_id(iz_west(zone_no));
      int ip_east   = zone_proc_id(iz_east(zone_no));
      int ip_south  = zone_proc_id(iz_south(zone_no));
      int ip_north  = zone_proc_id(iz_north(zone_no));


      if (ip_west != myid) {
        copy_x_face(&u(start5(iz)),
                    &qbc_in(qstart2_west(zone_no)),
                    nnx, nnxmax, nny, nnz, 0, DIR_IN);
      } else {
        int izone_west = iz_west(zone_no);
        copy_x_face(&u(start5(iz)),
                    //&qbc_ou(qstart_east(izone_west)),
                    &halo_array.local(0, proc_gzone_id(izone_west)-1, ZONE_EAST, 0),
                    nnx, nnxmax, nny, nnz, 0, DIR_IN);
      } // endif


      if (ip_east != myid) {
        copy_x_face(&u(start5(iz)),
                    &qbc_in(qstart2_east(zone_no)),
                    nnx, nnxmax, nny, nnz, nnx-1, DIR_IN);
      } else {
        int izone_east = iz_east(zone_no);
        copy_x_face(&u(start5(iz)),
                    //&qbc_ou(qstart_west(izone_east)),
                    &halo_array.local(0, proc_gzone_id(izone_east)-1, ZONE_WEST, 0),
                    nnx, nnxmax, nny, nnz, nnx-1, DIR_IN);
      } // endif


      if (ip_south != myid) {
          copy_y_face(&u(start5(iz)),
                      &qbc_in(qstart2_south(zone_no)),
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
      } else {
          int jzone_south = iz_south(zone_no);
          copy_y_face(&u(start5(iz)),
                      //&qbc_ou(qstart_north(jzone_south)),
                      &halo_array.local(0, proc_gzone_id(jzone_south)-1, ZONE_NORTH, 0),
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
      } // endif


      if (ip_north != myid) {
        copy_y_face(&u(start5(iz)),
                    &qbc_in(qstart2_north(zone_no)),
                    nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
      } else {
        int jzone_north = iz_north(zone_no);
        copy_y_face(&u(start5(iz)),
                    //&qbc_ou(qstart_south(jzone_north)),
                    &halo_array.local(0, proc_gzone_id(jzone_north)-1, ZONE_SOUTH, 0),
                    nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
      } // endif

  } // end do
  if (timeron) timer_stop(t_rdis1);

  return;
}

static void
copy_y_face(
  double* u_ptr,
  double* qbc_ptr,
  int nx,
  int nxmax,
  int ny,
  int nz,
  int jloc,
  int dir)
{
  UArrayT u(u_ptr, 5, nxmax, ny, nz);
  QBCArrayT qbc(qbc_ptr, 5, nx-2, nz-2);
/*
       subroutine copy_y_face(u, qbc, nx, nxmax, ny, nz, jloc, dir)

       implicit         none

       integer          nx, nxmax, ny, nz, i, j, k, jloc, m
       double precision u(5,0:nxmax-1,0:ny-1,0:nz-1), qbc(5,nx-2,nz-2)
       character        dir*(*)
*/
  int j = jloc;
  if (dir == DIR_IN) {
//!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m,i,k)
//!$OMP&  SHARED(j,nx,nz)
#pragma omp parallel for
    for (int k = 1; k <= nz-2; k++) {
      for (int i = 1; i <= nx-2; i++) {
        for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = qbc(m,i,k);
        } // end do
      } // end do
    } // end do
//!$OMP END PARALLEL DO
  } else if (dir == DIR_OUT) {
//!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m,i,k)
//!$OMP&  SHARED(j,nx,nz)
#pragma omp parallel for
    for (int k = 1; k <= nz-2; ++k) {
      for (int i = 1; i <= nx-2; i++) {
        for (int m = 1; m <= 5; ++m) {
          qbc(m,i,k) = u(m,i,j,k);
        } // end do
      } // end do
    } // end do
//!$OMP END PARALLEL DO
  } else {
    std::cout << "Erroneous data designation: " << dir << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 6);
  } // endif

  return;
}


static void
copy_x_face(
  double* u_ptr,
  double* qbc_ptr,
  int nx,
  int nxmax,
  int ny,
  int nz,
  int iloc,
  int dir)
{
  UArrayT u(u_ptr, 5, nxmax, ny, nz);
  QBCArrayT qbc(qbc_ptr, 5, ny-2, nz-2);
/*
       subroutine copy_x_face(u, qbc, nx, nxmax, ny, nz, iloc, dir)

       implicit         none

       integer          nx, nxmax, ny, nz, i, j, k, iloc, m
       double precision u(5,0:nxmax-1,0:ny-1,0:nz-1), qbc(5,ny-2,nz-2)
       character        dir*(*)
*/
  int i = iloc;
  if (dir == DIR_IN) {
//!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m,j,k)
//!$OMP&  SHARED(i,ny,nz)
#pragma omp parallel for
    for (int k = 1; k <= nz-2; ++k) {
      for (int j = 1; j <= ny-2; ++j) {
        for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = qbc(m,j,k);
        } // end do
      } // end do
    } // end do
//!$OMP END PARALLEL DO
  } else if (dir == DIR_OUT) {
//!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m,j,k)
// !$OMP&  SHARED(i,ny,nz)
#pragma omp parallel for
    for (int k = 1; k <= nz-2; ++k) {
      for (int j = 1; j <= ny-2; ++j) {
        for (int m = 1; m <= 5; ++m) {
          qbc(m,j,k) = u(m,i,j,k);
        } // end do
      } // end do
    } // end do
//!$OMP END PARALLEL DO
  } else {
    std::cout << "Erroneous data designation: " << dir << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 6);
  } // endif

  return;

}

