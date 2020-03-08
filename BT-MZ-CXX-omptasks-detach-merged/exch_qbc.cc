
#include<iostream>
#include<assert.h>
#include <mutex>

#include <atomic>
#include <omp.h>

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

static std::atomic<int> active_reqs{0};

static thread_local int complete_reqs = 0;

uint64_t num_sends = 0;

int max_completed = 0;

MPI_Request cont_req = MPI_REQUEST_NULL;

bool mpi_poll_service(void* data)
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

int request_completion_cb(void *data)
{
  int flag;
  omp_event_handle_t *event_ptr = static_cast<omp_event_handle_t*>(data);
  //--active_reqs;
  std::atomic_fetch_add_explicit(&active_reqs, -1, std::memory_order_relaxed);
  ++complete_reqs;
  omp_fulfill_event(*event_ptr);

  delete event_ptr;
  
  return MPI_SUCCESS;
}

static void wait_for_reqs(omp_event_handle_t event, MPI_Request *reqs)
{
  int flag = 0;
  omp_event_handle_t * event_ptr = new omp_event_handle_t{event};
  MPI_Continueall(2, reqs, &flag, event_ptr, MPI_STATUSES_IGNORE, cont_req);
  if (!flag) {
    //++active_reqs;
    auto ar = std::atomic_fetch_add_explicit(&active_reqs, 1, std::memory_order_relaxed) +1;
  } else {
    omp_fulfill_event(event);
    delete event_ptr;
  }
}



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

  for (int iz = 1; iz <= proc_num_zones; ++iz) {
      auto *u_ptr = &u(start5(iz));
      int zone_no = proc_zone_id(iz);
      int nnx    = nx(zone_no);
      int nnxmax = nxmax(zone_no);
      int nny    = ny(zone_no);
      int nnz    = nz(zone_no);
      omp_event_handle_t event;
      
      size_type x_face_size = (nny-2)*(nnz-2)*5;
      size_type y_face_size = (nnx-2)*(nnz-2)*5;

      /* send west: use only an input dependency, which will sync with the
       * output dependency at the corresponding adi task */
      {
      auto *out_buf = &qbc_ou(qstart_west(zone_no));
#pragma omp task depend(in:u_ptr[0])  \
                 depend(out: out_buf[0])
{
      copy_x_face(u_ptr,
                  out_buf,
                  nnx, nnxmax, nny, nnz, 1, DIR_OUT);
      #pragma omp taskwait
}
      int ip_west   = zone_proc_id(iz_west(zone_no));
      if (ip_west != myid) {
        auto *in_buf  = &qbc_in(qstart2_west(zone_no));
#pragma omp task depend(in: out_buf[0]) depend(out: in_buf[0]) detach(event)
{
        MPI_Request reqs[2];
        int stag = 10000 + zone_no;
        MPI_Isend(out_buf, x_face_size, MPI_DOUBLE,
                 ip_west, stag, comm_setup, &reqs[0]);
        int rtag = 20000 + iz_west(zone_no);
        MPI_Irecv(in_buf, x_face_size, MPI_DOUBLE,
                 ip_west, rtag, comm_setup, &reqs[1]);
        wait_for_reqs(event, reqs);

} // pragma omp task
      }
      }

      /* send east */
      {
      auto *out_buf = &qbc_ou(qstart_east(zone_no));
#pragma omp task depend(in:u_ptr[0])                       \
                 depend(out: out_buf[0])
{
      copy_x_face(u_ptr,
                  out_buf,
                  nnx, nnxmax, nny, nnz, nnx-2, DIR_OUT);
      #pragma omp taskwait
}
      int ip_east   = zone_proc_id(iz_east(zone_no));

      if (ip_east != myid) {
        auto *in_buf  = &qbc_in(qstart2_east(zone_no));
#pragma omp task depend(in: out_buf[0]) depend(out: in_buf[0]) detach(event)
{
        MPI_Request reqs[2];
        int stag = 20000 + zone_no;
        MPI_Isend(out_buf, x_face_size, MPI_DOUBLE,
                 ip_east, stag, comm_setup, &reqs[0]);
        int rtag = 10000 + iz_east(zone_no);
        MPI_Irecv(in_buf, x_face_size, MPI_DOUBLE,
                 ip_east, rtag, comm_setup, &reqs[1]);
        wait_for_reqs(event, reqs);
} // pragma omp task
      }
      }


      /* send south */
      {
      auto *out_buf = &qbc_ou(qstart_south(zone_no));
#pragma omp task depend(in:u_ptr[0])                  \
                 depend(out:out_buf[0])
{
      copy_y_face(u_ptr,
                  out_buf,
                  nnx, nnxmax, nny, nnz, 1, DIR_OUT);
      #pragma omp taskwait
}
      int ip_south   = zone_proc_id(iz_south(zone_no));

      if (ip_south != myid) {
        auto *in_buf  = &qbc_in(qstart2_south(zone_no));
#pragma omp task depend(in: out_buf[0]) depend(out:in_buf[0]) detach(event)
{
        MPI_Request reqs[2];
        int stag = 30000 + zone_no;
        MPI_Isend(out_buf, y_face_size, MPI_DOUBLE,
                  ip_south, stag, comm_setup, &reqs[0]);
        int rtag = 40000 + iz_south(zone_no);
        MPI_Irecv(in_buf, y_face_size, MPI_DOUBLE,
                 ip_south, rtag, comm_setup, &reqs[1]);
        wait_for_reqs(event, reqs);
} // pragma omp task
      }
      }


      /* send north */
      {
      auto *out_buf = &qbc_ou(qstart_north(zone_no));
#pragma omp task depend(in:u_ptr[0]) \
                 depend(out:out_buf[0])
{
      copy_y_face(u_ptr,
                  out_buf,
                  nnx, nnxmax, nny, nnz, nny-2, DIR_OUT);
      #pragma omp taskwait
}
      int ip_north   = zone_proc_id(iz_north(zone_no));

      if (ip_north != myid) {
        auto *in_buf  = &qbc_in(qstart2_north(zone_no));
#pragma omp task depend(in:out_buf[0]) depend(out:in_buf[0]) detach(event)
{
        MPI_Request reqs[2];
        int stag = 40000 + zone_no;
        MPI_Isend(out_buf, y_face_size, MPI_DOUBLE,
                 ip_north, stag, comm_setup, &reqs[0]);
        int rtag = 30000 + iz_north(zone_no);
        MPI_Irecv(in_buf, y_face_size, MPI_DOUBLE,
                 ip_north, rtag, comm_setup, &reqs[1]);
        wait_for_reqs(event, reqs);
} // pragma omp task
      }
      }

  } // end do



  for (int iz = 1; iz <= proc_num_zones; ++iz) {
      auto *u_ptr = &u(start5(iz));
      int zone_no = proc_zone_id(iz);
      int nnx    = nx(zone_no);
      int nnxmax = nxmax(zone_no);
      int nny    = ny(zone_no);
      int nnz    = nz(zone_no);

      int ip_west   = zone_proc_id(iz_west(zone_no));
      int ip_east   = zone_proc_id(iz_east(zone_no));
      int ip_south  = zone_proc_id(iz_south(zone_no));
      int ip_north  = zone_proc_id(iz_north(zone_no));

      size_type x_face_size = (nny-2)*(nnz-2)*5;
      size_type y_face_size = (nnx-2)*(nnz-2)*5;

      if (ip_west != myid) {
        auto *in_buf = &qbc_in(qstart2_west(zone_no));
        #pragma omp task depend(in:u_ptr[0]) \
                         depend(in:in_buf[0])
        {
          copy_x_face(u_ptr,
                      in_buf,
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
          #pragma omp taskwait
        } // omp task
      } else {
        int izone_west = iz_west(zone_no);
        auto *ou_buf = &qbc_ou(qstart_east(izone_west));
        #pragma omp task depend(in:u_ptr[0]) \
                         depend(in:ou_buf[0])
        {
          copy_x_face(u_ptr,
                      ou_buf,
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
          #pragma omp taskwait
        } // omp task
      } // endif

      if (ip_east != myid) {
        auto *in_buf = &qbc_in(qstart2_east(zone_no));
        #pragma omp task depend(in:u_ptr[0]) \
                         depend(in:in_buf[0])
        {
          copy_x_face(u_ptr,
                      in_buf,
                      nnx, nnxmax, nny, nnz, nnx-1, DIR_IN);
          #pragma omp taskwait
        } // omp task
      } else {
        int izone_east = iz_east(zone_no);
        auto *ou_buf = &qbc_ou(qstart_west(izone_east));
        #pragma omp task depend(in:u_ptr[0]) \
                         depend(in:ou_buf[0])
        {
          copy_x_face(u_ptr,
                      ou_buf,
                      nnx, nnxmax, nny, nnz, nnx-1, DIR_IN);
          #pragma omp taskwait
        } // omp task
      } // endif

      if (ip_south != myid) {
        auto *in_buf = &qbc_in(qstart2_south(zone_no));
        #pragma omp task depend(in:u_ptr[0]) \
                         depend(in:in_buf[0])
        {
          copy_y_face(u_ptr,
                      in_buf,
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
          #pragma omp taskwait
        } // omp task
      } else {
        int jzone_south = iz_south(zone_no);
        auto *ou_buf = &qbc_ou(qstart_north(jzone_south));
        #pragma omp task depend(in:u_ptr[0]) \
                         depend(in:ou_buf[0])
        {
          copy_y_face(u_ptr,
                      ou_buf,
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
          #pragma omp taskwait
        } // omp task
      } // endif

      if (ip_north != myid) {
        auto *in_buf =  &qbc_in(qstart2_north(zone_no));
        #pragma omp task depend(in:u_ptr[0]) \
                         depend(in:in_buf[0])
        {
          copy_y_face(u_ptr,
                      in_buf,
                      nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
          #pragma omp taskwait
        } // omp task
      } else {
        int jzone_north = iz_north(zone_no);
        auto *ou_buf = &qbc_ou(qstart_south(jzone_north));
        #pragma omp task depend(in:u_ptr[0]) \
                         depend(in:ou_buf[0])
        {
          copy_y_face(u_ptr,
                      ou_buf,
                      nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
          #pragma omp taskwait
        } // omp task
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

