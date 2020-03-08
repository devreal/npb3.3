
#include "header.h"

/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine exact_solution(xi,eta,zeta,dtemp)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     this function returns the exact solution at point xi, eta, zeta
c---------------------------------------------------------------------
*/

template<typename ValueT>
void exact_solution(
  ValueT xi, ValueT eta, ValueT zeta,
  ValueT* dtemp_ptr)
{
  StaticArrayViewT1<ValueT, size_type, 5, 1> dtemp(dtemp_ptr);
  // do m = 1, 5
  for (int m = 1; m <= 5; ++m) {
      dtemp(m) =  ce(m,1) +
        xi*(ce(m,2) + xi*(ce(m,5) + xi*(ce(m,8) + xi*ce(m,11)))) +
        eta*(ce(m,3) + eta*(ce(m,6) + eta*(ce(m,9) + eta*ce(m,12))))+
        zeta*(ce(m,4) + zeta*(ce(m,7) + zeta*(ce(m,10) +
        zeta*ce(m,13))));
  } // enddo
}


