#ifndef HAVE_EXCH_QBC_H
#define HAVE_EXCH_QBC_H

#include <mpi.h>
#include "npbparams_cc.h"

constexpr const int MAX_REQS = 32;
constexpr const int MSG_SIZE = 400000;

enum {
  DIR_IN  = 0,
  DIR_OUT = 1
};


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
  int iprn_msg);

#endif // HAVE_EXCH_QBC_H
