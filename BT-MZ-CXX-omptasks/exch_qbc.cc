
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

  int c_size, m_size, num_msgs;
  constexpr const int MSG_TAG = 10000;

  StaticArrayT1<MPI_Request, int, MAX_REQS, 1> requests;

  // copy data to qbc buffer
  if (timeron) timer_start(t_rdis1);

  // do iz = 1, proc_num_zones
  for (int iz = 1; iz <= proc_num_zones; ++iz) {
      auto *u_ptr = &u(start5(iz));
      ArrayT1<double, size_type>* qbc_ou_ = &qbc_ou;
      ArrayT1<double, size_type>* qbc_in_ = &qbc_in;
      int zone_no = proc_zone_id(iz);
      int nnx    = nx(zone_no);
      int nnxmax = nxmax(zone_no);
      int nny    = ny(zone_no);
      int nnz    = nz(zone_no);

#pragma omp task depend(in:u_ptr[0])
{
      auto &qbc_ou = *qbc_ou_;
      auto &qbc_in = *qbc_in_;
  
      copy_x_face(u_ptr,
                  &qbc_ou(qstart_west(zone_no)),
                  nnx, nnxmax, nny, nnz, 1, DIR_OUT);

      copy_x_face(u_ptr,
                  &qbc_ou(qstart_east(zone_no)),
                  nnx, nnxmax, nny, nnz, nnx-2, DIR_OUT);


      copy_y_face(u_ptr,
                  &qbc_ou(qstart_south(zone_no)),
                  nnx, nnxmax, nny, nnz, 1, DIR_OUT);

      copy_y_face(u_ptr,
                  &qbc_ou(qstart_north(zone_no)),
                  nnx, nnxmax, nny, nnz, nny-2, DIR_OUT);
#pragma omp taskwait
      print_zone(u_ptr, nnx, nnxmax, nny, nnz, zone_no, -10);
} // pragma omp task

  } // end do

#pragma omp taskwait

  if (timeron) timer_stop(t_rdis1);


  // exchange qbc buffers
  if (timeron) timer_start(t_rdis2);

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

  if (timeron) timer_stop(t_rdis2);


  // copy data from qbc buffer
  if (timeron) timer_start(t_rdis1);

  for (int iz = 1; iz <= proc_num_zones; ++iz) {
      auto *u_ptr = &u(start5(iz));
      ArrayT1<double, size_type>* qbc_ou_ = &qbc_ou;
      ArrayT1<double, size_type>* qbc_in_ = &qbc_in;
      int zone_no = proc_zone_id(iz);
      int nnx    = nx(zone_no);
      int nnxmax = nxmax(zone_no);
      int nny    = ny(zone_no);
      int nnz    = nz(zone_no);
      int ip_west   = zone_proc_id(iz_west(zone_no));
      int ip_east   = zone_proc_id(iz_east(zone_no));
      int ip_south  = zone_proc_id(iz_south(zone_no));
      int ip_north  = zone_proc_id(iz_north(zone_no));

#pragma omp task depend(out:u_ptr[0])
{
      auto &qbc_ou = *qbc_ou_;
      auto &qbc_in = *qbc_in_;
      if (ip_west != myid) {
        copy_x_face(u_ptr,
                    &qbc_in(qstart2_west(zone_no)),
                    nnx, nnxmax, nny, nnz, 0, DIR_IN);
      } else {
        int izone_west = iz_west(zone_no);
        copy_x_face(u_ptr,
                    &qbc_ou(qstart_east(izone_west)),
                    nnx, nnxmax, nny, nnz, 0, DIR_IN);
      } // endif

      if (ip_east != myid) {
        copy_x_face(u_ptr,
                    &qbc_in(qstart2_east(zone_no)),
                    nnx, nnxmax, nny, nnz, nnx-1, DIR_IN);
      } else {
        int izone_east = iz_east(zone_no);
        copy_x_face(u_ptr,
                    &qbc_ou(qstart_west(izone_east)),
                    nnx, nnxmax, nny, nnz, nnx-1, DIR_IN);
      } // endif

      if (ip_south != myid) {
          copy_y_face(u_ptr,
                      &qbc_in(qstart2_south(zone_no)),
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
      } else {
          int jzone_south = iz_south(zone_no);
          copy_y_face(u_ptr,
                      &qbc_ou(qstart_north(jzone_south)),
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
      } // endif

      if (ip_north != myid) {
        copy_y_face(u_ptr,
                    &qbc_in(qstart2_north(zone_no)),
                    nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
      } else {
        int jzone_north = iz_north(zone_no);
        copy_y_face(u_ptr,
                    &qbc_ou(qstart_south(jzone_north)),
                    nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
      } // endif

#pragma omp taskwait

      print_zone(u_ptr, nnx, nnxmax, nny, nnz, zone_no, -11);
} // pragma omp task
  } // end do
  if (timeron) timer_stop(t_rdis1);
#pragma omp taskwait
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
  int j = jloc;
  if (dir == DIR_IN) {
#pragma omp taskloop
    for (int k = 1; k <= nz-2; k++) {
      for (int i = 1; i <= nx-2; i++) {
        for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = qbc(m,i,k);
        } // end do
      } // end do
    } // end do
  } else if (dir == DIR_OUT) {
#pragma omp taskloop
    for (int k = 1; k <= nz-2; ++k) {
      for (int i = 1; i <= nx-2; i++) {
        for (int m = 1; m <= 5; ++m) {
          qbc(m,i,k) = u(m,i,j,k);
        } // end do
      } // end do
    } // end do
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
  int i = iloc;
  if (dir == DIR_IN) {
#pragma omp taskloop
    for (int k = 1; k <= nz-2; ++k) {
      for (int j = 1; j <= ny-2; ++j) {
        for (int m = 1; m <= 5; ++m) {
          u(m,i,j,k) = qbc(m,j,k);
        } // end do
      } // end do
    } // end do
  } else if (dir == DIR_OUT) {
#pragma omp taskloop
    for (int k = 1; k <= nz-2; ++k) {
      for (int j = 1; j <= ny-2; ++j) {
        for (int m = 1; m <= 5; ++m) {
          qbc(m,j,k) = u(m,i,j,k);
        } // end do
      } // end do
    } // end do
  } else {
    std::cout << "Erroneous data designation: " << dir << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 6);
  } // endif

  return;

}

