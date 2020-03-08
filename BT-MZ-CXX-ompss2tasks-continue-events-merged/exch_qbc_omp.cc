
#include<iostream>
#include<assert.h>
#include <mpi.h>
#include <atomic>
#include <nanos6.h>

#include "npbparams_cc.h"
#include "exch_qbc.h"
#include "print_zone.h"

void
copy_y_face_impl(
  double* u_ptr,
  double* qbc_ptr,
  int k,
  int nx,
  int nxmax,
  int ny,
  int nz,
  int jloc,
  int dir);

void
copy_x_face_impl(
  double* u_ptr,
  double* qbc_ptr,
  int k,
  int nx,
  int nxmax,
  int ny,
  int nz,
  int iloc,
  int dir);


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


static std::atomic<int> active_reqs{0};

static thread_local int complete_reqs = 0;

uint64_t num_sends = 0;

int max_completed = 0;

MPI_Request cont_req = MPI_REQUEST_NULL;

int mpi_poll_service(void* data)
{
  int completed = 0;
  if (active_reqs > 0)
  {
    do {
      int flag; // ignored
      complete_reqs = 0;
      MPI_Test(&cont_req, &flag, MPI_STATUS_IGNORE);
      completed += complete_reqs;
    } while (complete_reqs > 0 && active_reqs > 0);
    max_completed = std::max(max_completed, completed);
  }
  return false; /* signal that we need to continue to be called */
}

void request_completion_cb(void *task)
{
  //--active_reqs;
  auto ar = std::atomic_fetch_add_explicit(&active_reqs, -1, std::memory_order_relaxed);
  ++complete_reqs;
  nanos6_decrease_task_event_counter(task, 1);
}

static void wait_for_reqs(MPI_Request *reqs)
{
  MPI_Request tmp[2] = {reqs[0], reqs[1]};
  int flag = 0;
  void *task = nanos6_get_current_event_counter();
  
nanos6_increase_current_task_event_counter(task, 1);
  MPI_Continueall(2, reqs, &flag, &request_completion_cb, task, MPI_STATUSES_IGNORE, cont_req);
  if (!flag) {
    auto ar = std::atomic_fetch_add_explicit(&active_reqs, 1, std::memory_order_relaxed) +1;
  } else {
    nanos6_decrease_task_event_counter(task, 1);
  }
}



void exch_qbc(
  double* u,
  double* qbc_ou,
  double* qbc_in,
  size_type* nx,
  size_type* nxmax,
  size_type* ny,
  size_type* nz,
  size_type* start5,
  int      * proc_zone_id,
  size_type* qstart_west,
  size_type* qstart_east,
  size_type* qstart_south,
  size_type* qstart_north,
  size_type* qstart2_west,
  size_type* qstart2_east,
  size_type* qstart2_south,
  size_type* qstart2_north,
  size_type* iz_west,
  size_type* iz_east,
  size_type* iz_south,
  size_type* iz_north,
  MPI_Comm comm_setup,
  int myid,
  int* pcomm_group,
  size_type* qcomm_size,
  int* zone_proc_id,
  int num_procs,
  int proc_num_zones,
 int iprn_msg)
{

  // do iz = 1, proc_num_zones
  for (int iz = 0; iz < proc_num_zones; ++iz) {
      auto *u_ptr = &u[start5[iz]-1];
      int zone_no = proc_zone_id[iz]-1;
      int nnx    = nx[zone_no];
      int nnxmax = nxmax[zone_no];
      int nny    = ny[zone_no];
      int nnz    = nz[zone_no];
      size_type x_face_size = (nny-2)*(nnz-2)*5;
      size_type y_face_size = (nnx-2)*(nnz-2)*5;

      /* send west: use only an input dependency, which will sync with the
       * output dependency at the corresponding adi task */
#pragma oss task depend(in:u_ptr[0])  \
                 depend(out: qbc_ou[qstart_west[zone_no]-1]) no_copy_deps wait
{
      auto *out_buf = &qbc_ou[qstart_west[zone_no]-1];
      copy_x_face(&u[start5[iz]-1],
                  out_buf,
                  nnx, nnxmax, nny, nnz, 1, DIR_OUT);
}
      int ip_west   = zone_proc_id[iz_west[zone_no]-1];
      if (ip_west != myid) {
        auto *out_buf = &qbc_ou[qstart_west[zone_no]-1];
#pragma oss task depend(in: qbc_ou[qstart_west[zone_no]-1]) depend(out: qbc_in[qstart2_west[zone_no]-1]) no_copy_deps
{
        MPI_Request reqs[2];
        int stag = 10000 + zone_no;
        MPI_Isend(out_buf, x_face_size, MPI_DOUBLE,
                  ip_west, stag, comm_setup, &reqs[0]);
        int rtag = 20000 + iz_west[zone_no]-1;
        size_t roffset = qstart2_west[zone_no]-1;
        MPI_Irecv(&qbc_in[roffset], x_face_size, MPI_DOUBLE,
                  ip_west, rtag, comm_setup, &reqs[1]);

        wait_for_reqs(reqs);
} // pragma oss task
      }


      /* send east */
#pragma oss task depend(in:u_ptr[0])                       \
                 depend(out: qbc_ou[qstart_east[zone_no]-1]) no_copy_deps wait
{
      auto *out_buf = &qbc_ou[qstart_east[zone_no]-1];
      copy_x_face(&u[start5[iz]-1],
                  out_buf,
                  nnx, nnxmax, nny, nnz, nnx-2, DIR_OUT);
}
      int ip_east   = zone_proc_id[iz_east[zone_no]-1];

      if (ip_east != myid) {
        auto *out_buf = &qbc_ou[qstart_east[zone_no]-1];
#pragma oss task depend(in: qbc_ou[qstart_east[zone_no]-1]) depend(out: qbc_in[qstart2_east[zone_no]-1]) no_copy_deps
{
        MPI_Request reqs[2];
        int stag = 20000 + zone_no;
        MPI_Isend(out_buf, x_face_size, MPI_DOUBLE,
                  ip_east, stag, comm_setup, &reqs[0]);
        int rtag = 10000 + iz_east[zone_no]-1;
        size_t roffset = qstart2_east[zone_no]-1;
        MPI_Irecv(&qbc_in[roffset], x_face_size, MPI_DOUBLE,
                  ip_east, rtag, comm_setup, &reqs[1]);
        wait_for_reqs(reqs);
} // pragma oss task
      }


      /* send south */
#pragma oss task depend(in:u_ptr[0])                  \
                 depend(out: qbc_ou[qstart_south[zone_no]-1]) no_copy_deps wait
{
      auto *out_buf = &qbc_ou[qstart_south[zone_no]-1];
      copy_y_face(&u[start5[iz]-1],
                  out_buf,
                  nnx, nnxmax, nny, nnz, 1, DIR_OUT);
}
      int ip_south   = zone_proc_id[iz_south[zone_no]-1];

      if (ip_south != myid) {
        auto *out_buf = &qbc_ou[qstart_south[zone_no]-1];
#pragma oss task depend(in: qbc_ou[qstart_south[zone_no]-1]) depend(out: qbc_in[qstart2_south[zone_no]-1]) no_copy_deps
{
        MPI_Request reqs[2];
        int stag = 30000 + zone_no;
        MPI_Isend(out_buf, y_face_size, MPI_DOUBLE,
                  ip_south, stag, comm_setup, &reqs[0]);
        int rtag = 40000 + iz_south[zone_no]-1;
        size_t roffset = qstart2_south[zone_no]-1;
        MPI_Irecv(&qbc_in[roffset], y_face_size, MPI_DOUBLE,
                  ip_south, rtag, comm_setup, &reqs[1]);
        wait_for_reqs(reqs);
} // pragma oss task
      }


      /* send north */
#pragma oss task depend(in:u_ptr[0]) \
                 depend(out: qbc_ou[qstart_north[zone_no]-1]) no_copy_deps wait
{
      auto *out_buf = &qbc_ou[qstart_north[zone_no]-1];
      copy_y_face(&u[start5[iz]-1],
                  out_buf,
                  nnx, nnxmax, nny, nnz, nny-2, DIR_OUT);
}
      int ip_north   = zone_proc_id[iz_north[zone_no]-1];

      if (ip_north != myid) {
        auto *out_buf = &qbc_ou[qstart_north[zone_no]-1];
#pragma oss task depend(in: qbc_ou[qstart_north[zone_no]-1]) depend(out: qbc_in[qstart2_north[zone_no]-1]) no_copy_deps
{
        MPI_Request reqs[2];
        int stag = 40000 + zone_no;
        MPI_Isend(out_buf, y_face_size, MPI_DOUBLE,
                  ip_north, stag, comm_setup, &reqs[0]);
        int rtag = 30000 + iz_north[zone_no]-1;
        size_t roffset = qstart2_north[zone_no]-1;
        MPI_Irecv(&qbc_in[roffset], y_face_size, MPI_DOUBLE,
                  ip_north, rtag, comm_setup, &reqs[1]);
        wait_for_reqs(reqs);
} // pragma oss task
      }

  } // end do


  for (int iz = 0; iz < proc_num_zones; ++iz) {
      auto *u_ptr = &u[start5[iz]-1];
      int zone_no = proc_zone_id[iz]-1;
      int nnx    = nx[zone_no];
      int nnxmax = nxmax[zone_no];
      int nny    = ny[zone_no];
      int nnz    = nz[zone_no];

      int ip_west   = zone_proc_id[iz_west[zone_no]-1];
      int ip_east   = zone_proc_id[iz_east[zone_no]-1];
      int ip_south  = zone_proc_id[iz_south[zone_no]-1];
      int ip_north  = zone_proc_id[iz_north[zone_no]-1];

      size_type x_face_size = (nny-2)*(nnz-2)*5;
      size_type y_face_size = (nnx-2)*(nnz-2)*5;

      if (ip_west != myid) {
        #pragma oss task depend(in:u_ptr[0]) \
                        depend(in: qbc_in[qstart2_west[zone_no]-1]) no_copy_deps wait
        {
          auto *in_buf = &qbc_in[qstart2_west[zone_no]-1];
          copy_x_face(&u[start5[iz]-1],
                      in_buf,
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
        } // oss task
      } else {
        #pragma oss task depend(in:u_ptr[0]) \
                        depend(in: qbc_ou[qstart_east[iz_west[zone_no]-1]-1]) no_copy_deps wait
        {
          int izone_west = iz_west[zone_no]-1;
          copy_x_face(&u[start5[iz]-1],
                      &qbc_ou[qstart_east[izone_west]-1],
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
        } // oss task
      } // endif

      if (ip_east != myid) {
        #pragma oss task depend(in:u_ptr[0]) \
                        depend(in: qbc_in[qstart2_east[zone_no]-1]) no_copy_deps wait
        {
          auto *in_buf = &qbc_in[qstart2_east[zone_no]-1];
          copy_x_face(&u[start5[iz]-1],
                      in_buf,
                      nnx, nnxmax, nny, nnz, nnx-1, DIR_IN);
        } // oss task
      } else {
        #pragma oss task depend(in:u_ptr[0]) \
                         depend(in: qbc_ou[qstart_west[iz_east[zone_no]-1]-1]) no_copy_deps wait
        {
          //std::cout << "LCOPYIN zone " << zone_no+1 << " WEST" << std::endl;
          int izone_east = iz_east[zone_no]-1;
          copy_x_face(&u[start5[iz]-1],
                      &qbc_ou[qstart_west[izone_east]-1],
                      nnx, nnxmax, nny, nnz, nnx-1, DIR_IN);
        } // oss task
      } // endif

      if (ip_south != myid) {
        #pragma oss task depend(in:u_ptr[0]) \
                        depend(in: qbc_in[qstart2_south[zone_no]-1]) no_copy_deps wait
        {
          auto *in_buf = &qbc_in[qstart2_south[zone_no]-1];
          copy_y_face(&u[start5[iz]-1],
                      in_buf,
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
        } // oss task
      } else {
        #pragma oss task depend(in:u_ptr[0]) \
                        depend(in: qbc_ou[qstart_north[iz_south[zone_no]-1]-1]) no_copy_deps wait
        {
          //std::cout << "LCOPYIN zone " << zone_no+1 << " WEST" << std::endl;
          int jzone_south = iz_south[zone_no]-1;
          copy_y_face(&u[start5[iz]-1],
                      &qbc_ou[qstart_north[jzone_south]-1],
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
        } // oss task
      } // endif

      if (ip_north != myid) {
        #pragma oss task depend(in:u_ptr[0]) \
                        depend(in: qbc_in[qstart2_north[zone_no]-1]) no_copy_deps wait
        {
          auto *in_buf =  &qbc_in[qstart2_north[zone_no]-1];
          copy_y_face(&u[start5[iz]-1],
                      in_buf,
                      nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
        } // oss task
      } else {
        #pragma oss task depend(in:u_ptr[0]) \
                        depend(in: qbc_ou[qstart_south[iz_north[zone_no]-1]-1]) no_copy_deps wait
        {
          int jzone_north = iz_north[zone_no]-1;
          copy_y_face(&u[start5[iz]-1],
                      &qbc_ou[qstart_south[jzone_north]-1],
                      nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
        } // oss task
      } // endif
  } // end do


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
#pragma oss task for
  for (int k = 1; k <= nz-2; k++) {
    copy_y_face_impl(u_ptr, qbc_ptr, k, nx, nxmax, ny, nz, jloc, dir);
  } // end do

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
#pragma oss task for
  for (int k = 1; k <= nz-2; k++) {
    copy_x_face_impl(u_ptr, qbc_ptr, k, nx, nxmax, ny, nz, iloc, dir);
  } // end do

  return;
}

