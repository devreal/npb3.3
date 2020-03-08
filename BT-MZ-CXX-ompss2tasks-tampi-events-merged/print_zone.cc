
#include <iostream>
#include <iomanip>
#include "NArray.h"
#include "header.h"
#include "print_zone.h"

void print_zone(
  double* u_ptr,
  size_type nx, size_type nxmax, size_type ny, size_type nz, int zone, int step)
{
#ifdef DO_PRINT
  ArrayViewT4<double, size_type> u(u_ptr, 5, nxmax, ny, nz);

  std::cout << "------ Printing U zone " << zone << " at step " << step
            << " -------" << std::fixed << std::endl;
  for (int k = 0; k <= nz-1; ++k) {
    for (int j = 0; j <= ny-1; ++j) {
      for (int i = 0; i <= nx-1; ++i) {
        for (int m = 1; m <= 5; ++m) {
          std::cout << std::setw( 20 ) << std::setprecision( 13 ) << u(m, i, j, k);
        }
        std::cout << std::endl;
      }
    }
  }
  std::cout << "------ End of U ------- \n\n" << std::endl;
#endif
}

