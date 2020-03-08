
#include<iostream>
#include<assert.h>
#include <mpi.h>
#include <TAMPI.h>

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

/*
       subroutine exch_qbc(u, qbc_ou, qbc_in, nx, nxmax, ny, nz,
     $                     iprn_msg)

       include 'header.h'
       include 'mpi_stuff.h'

       integer   nx(*), nxmax(*), ny(*), nz(*), iprn_msg
       double precision u(*), qbc_ou(*), qbc_in(*)
*/
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

//#pragma oss task depend(out:u_ptr[0]) no_copy_deps
//{
//      print_zone(u_ptr, nnx, nnxmax, nny, nnz, zone_no+1, -10);
//}

//      std::cout << myid << ": Local zone " << iz << " is global zone "
//                << zone_no << " with size " << x_face_size << "x" 
//                << y_face_size << ", " << nnx << " x" << nny << " x" << nnz << ", neighbors " << iz_west[zone_no]-1 << " " << iz_east[zone_no]-1 << " " << iz_south[zone_no]-1 << " " << iz_north[zone_no]-1 << std::endl; 

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
        //std::cout << myid << ": Send " << x_face_size << "B, dest " << ip_west << ", tag " << stag << std::endl;
        MPI_Isend(out_buf, x_face_size, MPI_DOUBLE,
                  ip_west, stag, comm_setup, &reqs[0]);
        int rtag = 20000 + iz_west[zone_no]-1;
        //std::cout << myid << ": Recv " << x_face_size << "B, src " << ip_west << ", tag " << rtag << std::endl;
        size_t roffset = qstart2_west[zone_no]-1;
        MPI_Irecv(&qbc_in[roffset], x_face_size, MPI_DOUBLE,
                  ip_west, rtag, comm_setup, &reqs[1]);

        TAMPI_Iwaitall(2, reqs, MPI_STATUSES_IGNORE);
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
#if 0
      std::cout << "COPYOUT zone " << zone_no+1 << " EAST, sending to " << ip_east << ", receiving from zone " << iz_east[zone_no] << std::endl;
      for (int i = 0; i < x_face_size; ++i) {
        std::cout << out_buf[i]  << " ";
      }
      std::cout << std::endl;
#endif

      if (ip_east != myid) {
        auto *out_buf = &qbc_ou[qstart_east[zone_no]-1];
#pragma oss task depend(in: qbc_ou[qstart_east[zone_no]-1]) depend(out: qbc_in[qstart2_east[zone_no]-1]) no_copy_deps
{
        MPI_Request reqs[2];
        int stag = 20000 + zone_no;
        //std::cout << myid << ": Send " << x_face_size << "B, dest " << ip_east << ", tag " << stag << std::endl;
        MPI_Isend(out_buf, x_face_size, MPI_DOUBLE,
                  ip_east, stag, comm_setup, &reqs[0]);
        int rtag = 10000 + iz_east[zone_no]-1;
        //std::cout << myid << ": Recv " << x_face_size << "B, src " << ip_east << ", tag " << rtag << std::endl;
        size_t roffset = qstart2_east[zone_no]-1;
        MPI_Irecv(&qbc_in[roffset], x_face_size, MPI_DOUBLE,
                  ip_east, rtag, comm_setup, &reqs[1]);
        TAMPI_Iwaitall(2, reqs, MPI_STATUSES_IGNORE);
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
#if 0
      std::cout << "COPYOUT zone " << zone_no+1 << " SOUTH, sending to " << ip_south << ", receiving from zone " << iz_south[zone_no] << std::endl;
      for (int i = 0; i < y_face_size; ++i) {
        std::cout << out_buf[i]  << " ";
      }
      std::cout << std::endl;
#endif
      if (ip_south != myid) {
        auto *out_buf = &qbc_ou[qstart_south[zone_no]-1];
#pragma oss task depend(in: qbc_ou[qstart_south[zone_no]-1]) depend(out: qbc_in[qstart2_south[zone_no]-1]) no_copy_deps
{
        MPI_Request reqs[2];
        int stag = 30000 + zone_no;
        //std::cout << myid << ": Send " << y_face_size << "B, dest " << ip_south << ", tag " << stag << std::endl;
        MPI_Isend(out_buf, y_face_size, MPI_DOUBLE,
                  ip_south, stag, comm_setup, &reqs[0]);
        int rtag = 40000 + iz_south[zone_no]-1;
        size_t roffset = qstart2_south[zone_no]-1;
        //std::cout << myid << ": Recv " << y_face_size << "B, src " << ip_south << ", tag " << rtag << std::endl;
        MPI_Irecv(&qbc_in[roffset], y_face_size, MPI_DOUBLE,
                  ip_south, rtag, comm_setup, &reqs[1]);
        TAMPI_Iwaitall(2, reqs, MPI_STATUSES_IGNORE);
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
#if 0
      std::cout << "COPYOUT zone " << zone_no+1 << " NORTH, sending to " << ip_north << ", receiving from zone " << iz_north[zone_no] << std::endl;
      for (int i = 0; i < y_face_size; ++i) {
        std::cout << out_buf[i]  << " ";
      }
      std::cout << std::endl;
#endif

      if (ip_north != myid) {
        auto *out_buf = &qbc_ou[qstart_north[zone_no]-1];
#pragma oss task depend(in: qbc_ou[qstart_north[zone_no]-1]) depend(out: qbc_in[qstart2_north[zone_no]-1]) no_copy_deps
{
        MPI_Request reqs[2];
        int stag = 40000 + zone_no;
        //std::cout << myid << ": Send " << y_face_size << "B, dest " << ip_north << ", tag " << stag << std::endl;
        MPI_Isend(out_buf, y_face_size, MPI_DOUBLE,
                  ip_north, stag, comm_setup, &reqs[0]);
        int rtag = 30000 + iz_north[zone_no]-1;
        //std::cout << myid << ": Recv " << y_face_size << "B, src " << ip_north << ", tag " << rtag << std::endl;
        size_t roffset = qstart2_north[zone_no]-1;
        MPI_Irecv(&qbc_in[roffset], y_face_size, MPI_DOUBLE,
                  ip_north, rtag, comm_setup, &reqs[1]);
        // let tampi do the task magic
        TAMPI_Iwaitall(2, reqs, MPI_STATUSES_IGNORE);
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
#if 0
          std::cout << "COPYIN zone " << zone_no+1 << " WEST" << std::endl;
          for (int i = 0; i < x_face_size; ++i) {
            std::cout << in_buf[i]  << " ";
          }
          std::cout << std::endl;
#endif
          copy_x_face(&u[start5[iz]-1],
                      in_buf,
                      nnx, nnxmax, nny, nnz, 0, DIR_IN);
        } // oss task
      } else {
        #pragma oss task depend(in:u_ptr[0]) \
                        depend(in: qbc_ou[qstart_east[iz_west[zone_no]-1]-1]) no_copy_deps wait
        {
          //std::cout << "LCOPYIN zone " << zone_no+1 << " WEST" << std::endl;
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
#if 0
          std::cout << "COPYIN zone " << zone_no+1 << " EAST" << std::endl;
          for (int i = 0; i < x_face_size; ++i) {
            std::cout << in_buf[i]  << " ";
          }
          std::cout << std::endl;
#endif
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
#if 0
          std::cout << "COPYIN zone " << zone_no+1 << " SOUTH" << std::endl;
          for (int i = 0; i < y_face_size; ++i) {
            std::cout << in_buf[i]  << " ";
          }
          std::cout << std::endl;
#endif

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
#if 0
          std::cout << "COPYIN zone " << zone_no+1 << " NORTH" << std::endl;
          for (int i = 0; i < y_face_size; ++i) {
            std::cout << in_buf[i]  << " ";
          }
          std::cout << std::endl;
#endif
          copy_y_face(&u[start5[iz]-1],
                      in_buf,
                      nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
        } // oss task
      } else {
        #pragma oss task depend(in:u_ptr[0]) \
                        depend(in: qbc_ou[qstart_south[iz_north[zone_no]-1]-1]) no_copy_deps wait
        {
          //std::cout << "LCOPYIN zone " << zone_no+1 << " WEST" << std::endl;
          int jzone_north = iz_north[zone_no]-1;
          copy_y_face(&u[start5[iz]-1],
                      &qbc_ou[qstart_south[jzone_north]-1],
                      nnx, nnxmax, nny, nnz, nny-1, DIR_IN);
        } // oss task
      } // endif
//#pragma oss task depend(inout:u_ptr[0]) no_copy_deps
//{
//      print_zone(&u[start5[iz]-1], nnx, nnxmax, nny, nnz, zone_no+1, -11);
//}


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
/*
       subroutine copy_y_face(u, qbc, nx, nxmax, ny, nz, jloc, dir)

       implicit         none

       integer          nx, nxmax, ny, nz, i, j, k, jloc, m
       double precision u(5,0:nxmax-1,0:ny-1,0:nz-1), qbc(5,nx-2,nz-2)
       character        dir*(*)
*/

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
/*
       subroutine copy_x_face(u, qbc, nx, nxmax, ny, nz, iloc, dir)

       implicit         none

       integer          nx, nxmax, ny, nz, i, j, k, iloc, m
       double precision u(5,0:nxmax-1,0:ny-1,0:nz-1), qbc(5,ny-2,nz-2)
       character        dir*(*)
*/
#pragma oss task for
  for (int k = 1; k <= nz-2; k++) {
    copy_x_face_impl(u_ptr, qbc_ptr, k, nx, nxmax, ny, nz, iloc, dir);
  } // end do

  return;
}

