#ifndef HAVE_MPI_STUFF_H
#define HAVE_MPI_STUFF_H

#include <mpi.h>

#include "header.h"
#include "NArray.h"

/*
c     zone_proc_id(MZ)     - process id each zone assigned to
c     proc_zone_id(MZ)     - list of zones assigned to this process
c     proc_num_zones       - number of zones assigned to this process
c     proc_zone_count(NP)  - number of zones assigned to each process
c     proc_num_threads(NP) - number of threads assigned to each process
c     proc_group(NP)       - group id each process assigned to
c
*/


/*
      integer   zone_proc_id(max_zones), proc_zone_id(max_zones),
     &          proc_num_zones, proc_zone_count(num_procs),
     &          proc_num_threads(num_procs), proc_group(num_procs)
 */
extern int proc_num_zones;
extern ArrayT1<int, int> zone_proc_id;
extern ArrayT1<int, int> proc_zone_id;
extern ArrayT1<int, int> proc_zone_count;
extern ArrayT1<int, int> proc_num_threads;
extern ArrayT1<int, int> proc_group;

/*
 * double precision proc_zone_size(num_procs)
 */
extern ArrayT1<double, int> proc_zone_size;

/*
      common /mpi_cmn1/ proc_zone_size, proc_zone_id, zone_proc_id,
     &          proc_zone_count, proc_num_threads, proc_num_zones,
     &          proc_group
c
      integer   myid, root, comm_setup, ierror, dp_type
      integer   num_threads, mz_bload, max_threads
      logical   active
      common /mpi_cmn2/ myid, root, comm_setup, ierror, active,
     &          dp_type, num_threads, mz_bload, max_threads
c
*/

extern int myid, root;
extern MPI_Comm comm_setup;
extern int num_threads, mz_bload, max_threads;
extern bool active;

/*
c ... Two adjustable parameters for MPI communication
c     max_reqs  -- max. number of async message requests
c     MSG_SIZE  -- optimal message size (in words) for communication
      integer   max_reqs, MSG_SIZE
      parameter (max_reqs=32, MSG_SIZE=400000)
*/
constexpr const int MAX_REQS = 32;
constexpr const int MSG_SIZE = 400000;

/*
c
      integer   requests(max_reqs), statuses(MPI_STATUS_SIZE,max_reqs)
      common /mpi_cmn3/ requests, statuses
c
*/

/*
      integer   pcomm_group(num_procs2)
      dimension qcomm_size(num_procs),
     &          qstart2_west (max_zones), qstart2_east (max_zones),
     &          qstart2_south(max_zones), qstart2_north(max_zones)
      common /mpi_cmn4/ qcomm_size, qstart2_west, qstart2_east,
     &          qstart2_south, qstart2_north, pcomm_group
*/

extern ArrayT1<int, int> pcomm_group;
extern ArrayT1<int, int>
                      qcomm_size, qstart2_west , qstart2_east,
                      qstart2_south, qstart2_north;

void mpi_setup();

extern MPI_Request cont_req;
int request_completion_cb(void *data);


void map_zones(
  int num_zones,
  ArrayT1<size_type, size_type>& nx,
  ArrayT1<size_type, size_type>& ny,
  ArrayT1<size_type, size_type>& nz,
  int tot_threads);

void env_setup(int tot_threads);


void exch_qbc(
  ArrayT1<double, size_type>& u,
  ArrayT1<double, size_type>& qbc_ou,
  ArrayT1<double, size_type>& qbc_in,
  ArrayT1<size_type, size_type>& nx,
  ArrayT1<size_type, size_type>& nxmax,
  ArrayT1<size_type, size_type>& ny,
  ArrayT1<size_type, size_type>& nz,
  int iprn_msg);

#endif // HAVE_MPI_STUFF_H
