#include <iostream>
#include "npbparams_cc.h"
#include "adi.h"
#include "exch_qbc.h"
#include "print_zone.h"

constexpr static const int root = 0;

void do_timesteps(
  int niter,
  double *u,
  double *us,
  double *vs,
  double *ws,
  double *qs,
  double *rho_i,
  double *square,
  double *rhs,
  double *forcing,
  double *qbc_ou,
  double *qbc_in,
  size_type* nx,
  size_type* nxmax,
  size_type* ny,
  size_type* nz,
  size_type* start1,
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

    //doloop(1, niter, [&](auto step) {
    for (int step = 1; step <= niter; ++step) {
      exch_qbc(u, qbc_ou, qbc_in, nx, nxmax, ny, nz, start5,
               proc_zone_id, qstart_west, qstart_east, qstart_south, qstart_north,
               qstart2_west, qstart2_east, qstart2_south, qstart2_north,
               iz_west, iz_east, iz_south, iz_north,
               comm_setup, myid, pcomm_group, qcomm_size, zone_proc_id,
               num_procs, proc_num_zones, iprn_msg);

      //doloop(1, proc_num_zones, [&](auto iz){
      for (int iz = 0; iz < proc_num_zones; ++iz) {
        auto *u_ptr = &u[start5[iz]-1];
#pragma oss task depend(out:u_ptr[0]) wait
{
      if (iz == 0 && (step % 20 == 0 || step == 1)) {
        if (myid == root) {
          std::cout << "Time step " << step << std::endl;
        }
      } // endif


        int zone = proc_zone_id[iz]-1;
        print_zone(&u[start5[iz]-1], nx[zone], nxmax[zone], ny[zone], nz[zone], zone+1, step);


        adi(&rho_i[start1[iz]-1], &us[start1[iz]-1],
            &vs[start1[iz]-1], &ws[start1[iz]-1],
            &qs[start1[iz]-1], &square[start1[iz]-1],
            &rhs[start5[iz]-1], &forcing[start5[iz]-1],
            &u[start5[iz]-1],
            nx[zone], nxmax[zone], ny[zone], nz[zone]);
}
      }//); // end do

    } //); // end do

#pragma oss taskwait

}
