#ifndef HAVE_ZONE_SETUP_H
#define HAVE_ZONE_SETUP_H

#include <cmath>

#include "NArray.h"

#include "mpi_stuff.h"

void
zone_setup(
  ArrayT1<size_type, size_type>& nx,
  ArrayT1<size_type, size_type>& nxmax,
  ArrayT1<size_type, size_type>& ny,
  ArrayT1<size_type, size_type>& nz);


void
zone_starts(
  int num_zones,
  ArrayT1<size_type, size_type>& nx,
  ArrayT1<size_type, size_type>& nxmax,
  ArrayT1<size_type, size_type>& ny,
  ArrayT1<size_type, size_type>& nz);

#endif // HAVE_ZONE_SETUP_H
