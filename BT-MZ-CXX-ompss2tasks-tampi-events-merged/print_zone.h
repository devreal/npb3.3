#ifndef PRINT_ZONE_H
#define PRINT_ZONE_H

// uncomment to enable dumping of buffers/zones
//#define DO_PRINT

void print_zone(
  double* u_ptr,
  size_type nx, size_type nxmax, size_type ny, size_type nz, int zone, int step);

#endif // PRINT_ZONE_H
