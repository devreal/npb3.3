
#include <iostream>

#include "NArray.h"

#include "mpi_stuff.h"

#include "header.h"

#include "zone_setup.h"

#include "zone_setup.h"

StaticArrayT1<size_type, size_type, x_zones, 1>
                x_start, x_end, x_size;
StaticArrayT1<size_type, size_type, y_zones, 1>
                y_start, y_end, y_size;
StaticArrayT1<size_type, size_type, max_zones, 1>
                iz_west , iz_east ,
                iz_south, iz_north;

StaticArrayT1<size_type, size_type, max_zones, 1>
                start1, start5,
                qstart_west , qstart_east,
                qstart_south, qstart_north;


void
zone_setup(
  ArrayT1<size_type, size_type>& nx,
  ArrayT1<size_type, size_type>& nxmax,
  ArrayT1<size_type, size_type>& ny,
  ArrayT1<size_type, size_type>& nz)
{
  if (dabs(ratio-1.0) > 1.e-10) {

    // compute zone stretching only if the prescribed zone size ratio
    // is substantially larger than unity

    double x_r   = dexp(dlog(ratio)/(x_zones-1));
    double y_r   = dexp(dlog(ratio)/(y_zones-1));
    double x_smallest = dble(gx_size)*(x_r-1.0)/(std::pow(x_r, x_zones)-1.0);
    double y_smallest = dble(gy_size)*(y_r-1.0)/(std::pow(y_r, y_zones)-1.0);

    // compute tops of intervals, using a slightly tricked rounding
    // to make sure that the intervals are increasing monotonically
    // in size

    //do i = 1, x_zones
    for (int i = 1; i <= x_zones; ++i) {
      x_end(i) = x_smallest*(std::pow(x_r, i)-1.0)/(x_r-1.0)+0.45;
    } // end do

    //do j = 1, y_zones
    for (int j = 1; j <= y_zones; ++j) {
      y_end(j) = y_smallest*(std::pow(y_r, j)-1.0)/(y_r-1.0)+0.45;
    } // end do

  } else {

    // compute essentially equal sized zone dimensions

    //do i = 1, x_zones
    for (int i = 1; i <= x_zones; ++i) {
      x_end(i)   = (i*gx_size)/x_zones;
    } // end do

    // do j = 1, y_zones
    for (int j = 1; j <= y_zones; ++j) {
      y_end(j)   = (j*gy_size)/y_zones;
    } // end do
  }

  x_start(1) = 1;
  //do i = 1, x_zones
  for (int i = 1; i <= x_zones; ++i) {
    if (i != x_zones) x_start(i+1) = x_end(i) + 1;
    x_size(i)  = x_end(i) - x_start(i) + 1;
  } // end do

  y_start(1) = 1;
  //do j = 1, y_zones
  for (int j = 1; j <= y_zones; ++j) {
    if (j != y_zones) y_start(j+1) = y_end(j) + 1;
    y_size(j) = y_end(j) - y_start(j) + 1;
  } // end do

  //if (npb_verbose .gt. 1 .and. myid .eq. root) then
  if (npb_verbose == 1 && myid == root) {
    std::cout << " Zone sizes: " << std::endl;
  }

  // do j = 1, y_zones
  for (int j = 1; j <= y_zones; ++j) {
    //do i = 1, x_zones
    for (int i = 1; i <= x_zones; ++i) {
      int zone_no;
      zone_no = (i-1)+(j-1)*x_zones+1;
      nx(zone_no) = x_size(i);
      nxmax(zone_no) = nx(zone_no) + 1 - mod(nx(zone_no),2);
      ny(zone_no) = y_size(j);
      nz(zone_no) = gz_size;


      int id_west, id_east, jd_south, jd_north;
      id_west  = (i-2+x_zones) % x_zones;
      id_east  = i             % x_zones;
      jd_south = (j-2+y_zones) % y_zones;
      jd_north = j             % y_zones;
      iz_west (zone_no) = id_west +  (j-1)*x_zones + 1;
      iz_east (zone_no) = id_east +  (j-1)*x_zones + 1;
      iz_south(zone_no) = (i-1) + jd_south*x_zones + 1;
      iz_north(zone_no) = (i-1) + jd_north*x_zones + 1;

      if (npb_verbose == 1 && myid == root) {
        std::cout << zone_no << ":  " << nx(zone_no)
                             << " x"  << ny(zone_no)
                             << " x"  << nz(zone_no) << std::endl;
      }
    } // end do
  } // end do

  return;
}


/*
       subroutine zone_starts(num_zones, nx, nxmax, ny, nz)

       include 'header.h'
       include 'mpi_stuff.h'

       integer   num_zones
       integer   nx(*), nxmax(*), ny(*), nz(*)

       integer   zone, zone_size, iz, ip, zone2, ig, id, itmp
       integer   x_face_size, y_face_size
       integer   ip_west, ip_east, ip_south, ip_north
*/

void
zone_starts(
  int num_zones,
  ArrayT1<size_type, size_type>& nx,
  ArrayT1<size_type, size_type>& nxmax,
  ArrayT1<size_type, size_type>& ny,
  ArrayT1<size_type, size_type>& nz)
{

  int ip_east, ip_west, ip_south, ip_north;
  size_type x_face_size;
  size_type y_face_size;
  //do iz = 1, proc_num_zones
  DO(iz, 1, proc_num_zones, {
    int zone = proc_zone_id(iz);
    size_type zone_size = nxmax(zone)*ny(zone)*nz(zone);
    if (iz == 1) {
      start1(iz) = 1;
      start5(iz) = 1;
    }
    if (iz != proc_num_zones) {
      start1(iz+1) = start1(iz) + zone_size;
      start5(iz+1) = start5(iz) + zone_size*5;
    } else {
      if (start1(iz)+zone_size-1 > proc_max_size) {
        std::cout << " Error in size: zone " << zone
                  << " proc_max_size " << proc_max_size
                  << " access size " <<  start1(iz)+zone_size-1 << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
  }); // enddo

  size_type qoffset = 1;
  //do 10 ip = 0, num_procs-1
  for (int ip = 0; ip <num_procs; ++ip) {

    if (ip == myid) continue;

    //do 15 zone = 1, num_zones
    for (int zone = 1; zone <= num_zones; ++zone) {
      if (zone_proc_id(zone) != ip) continue;

      x_face_size = (ny(zone)-2)*(nz(zone)-2)*5;
      y_face_size = (nx(zone)-2)*(nz(zone)-2)*5;

      int zone2 = iz_west(zone);
      ip_east  = zone_proc_id(zone2);
      if (ip_east == myid) {
        qstart2_east(zone2) = qoffset;
        qoffset = qoffset + x_face_size;
      }

      zone2 = iz_east(zone);
      ip_west  = zone_proc_id(zone2);
      if (ip_west == myid) {
        qstart2_west(zone2) = qoffset;
        qoffset = qoffset + x_face_size;
      }

      zone2 = iz_south(zone);
      ip_north = zone_proc_id(zone2);
      if (ip_north == myid) {
        qstart2_north(zone2) = qoffset;
        qoffset = qoffset + y_face_size;
      }

      zone2 = iz_north(zone);
      ip_south = zone_proc_id(zone2);
      if (ip_south == myid) {
        qstart2_south(zone2) = qoffset;
        qoffset = qoffset + y_face_size;
      }
    }
  }

  qoffset = 1;
  //do 20 ip = 0, num_procs-1
  for (int ip = 0; ip < num_procs; ++ip) {

    if (ip != myid) {

      //do 30 zone = 1, num_zones
      for (int zone = 1; zone <= num_zones; ++zone) {
        if (zone_proc_id(zone) != myid) continue;

        ip_west  = zone_proc_id(iz_west(zone));
        ip_east  = zone_proc_id(iz_east(zone));
        ip_south = zone_proc_id(iz_south(zone));
        ip_north = zone_proc_id(iz_north(zone));

        x_face_size = (ny(zone)-2)*(nz(zone)-2)*5;
        y_face_size = (nx(zone)-2)*(nz(zone)-2)*5;

        if (ip_west == ip) {
          qstart_west(zone) = qoffset;
          qoffset = qoffset + x_face_size;
        }

        if (ip_east == ip) {
          qstart_east(zone) = qoffset;
          qoffset = qoffset + x_face_size;
        }

        if (ip_south == ip) {
          qstart_south(zone) = qoffset;
          qoffset = qoffset + y_face_size;
        }

        if (ip_north == ip) {
          qstart_north(zone) = qoffset;
          qoffset = qoffset + y_face_size;
        }
      }
    }
    qcomm_size(ip+1) = qoffset - 1;
//20  continue
  }

//c ...  for intra-process zone copy
  //do 40 zone = 1, num_zones
  for (int zone = 1; zone <= num_zones; ++zone) {
    if (zone_proc_id(zone) != myid) continue;

    ip_west  = zone_proc_id(iz_west(zone));
    ip_east  = zone_proc_id(iz_east(zone));
    ip_south = zone_proc_id(iz_south(zone));
    ip_north = zone_proc_id(iz_north(zone));

    x_face_size = (ny(zone)-2)*(nz(zone)-2)*5;
    y_face_size = (nx(zone)-2)*(nz(zone)-2)*5;

    if (ip_west == myid) {
      qstart_west(zone)  = qoffset;
      qoffset = qoffset + x_face_size;
    }
    if (ip_east == myid) {
      qstart_east(zone)  = qoffset;
      qoffset = qoffset + x_face_size;
    }
    if (ip_south == myid) {
      qstart_south(zone) = qoffset;
      qoffset = qoffset + y_face_size;
    }
    if (ip_north == myid) {
      qstart_north(zone) = qoffset;
      qoffset = qoffset + y_face_size;
    }
  }

//c ...  set up cyclic communication group
  int iz = 1;
  //do while (iz .lt. myid)
  while (iz < myid) {
    iz = iz * 2;
  }
  if (iz > myid) iz = iz / 2;

  //do ig = 1, num_procs
  for (int ig = 1; ig <= num_procs; ++ig) {
    pcomm_group(ig) = ig - 1;
  } // enddo
  int ig = num_procs + 1;
  //do while (ig .le. num_procs2)
  while (ig <= num_procs2) {
    pcomm_group(ig) = -1;
    ig = ig + 1;
  }

  int id = 0;
  //do while (iz .ge. 1)
  while (iz >= 1) {
    if (id+iz <= myid) {
      //do ig = 1, num_procs, iz*2
      for (int ig = 1; ig <= num_procs; ig += iz*2) {
        // do ip = ig, ig+iz-1
        for (int ip = ig; ip <= ig+iz-1; ++ip) {
          int itmp = pcomm_group(ip);
          pcomm_group(ip) = pcomm_group(ip+iz);
          pcomm_group(ip+iz) = itmp;
        }
      } // enddo
      id = id + iz;
    }
    iz = iz/2;
  }

  int ip = 1;
  //do ig = 1, num_procs
  for (ig = 1; ig <= num_procs; ++ig) {
    //do while (pcomm_group(ip) .lt. 0)
    while (pcomm_group(ip) < 0) {
      ip = ip + 1;
    } // enddo
    pcomm_group(ig) = pcomm_group(ip);
    ip = ip + 1;
  } // enddo

  if (npb_verbose > 1) {
    //do iz = 1, proc_num_zones
    for (int iz = 1; iz <= proc_num_zones; ++iz) {
      int zone = proc_zone_id(iz);
      std::cout << " myid " << myid << " iz="  << iz
                << " zone=" << zone << " start1=" << start1(iz)
                << " start5=" << start5(iz) << std::endl;
    } // enddo
    std::cout << " myid " << myid << " qcomm_size=" << qoffset-1 << std::endl;

    //do ig = 1, num_procs
    for (int ig = 1; ig <= num_procs; ++ig) {
      int ip = pcomm_group(ig) + 1;
      if (ip == 1) {
        qoffset = qcomm_size(ip);
      } else {
        qoffset = qcomm_size(ip) - qcomm_size(ip-1);
      }
      std::cout << " myid " << myid << " proc " << ip-1 << " qcomm_size " << qcomm_size(ip) << std::endl;
    } // enddo
  }

  return;
}
