
/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  work_lhs.h
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------
*/

#include "NArray.h"

StaticArrayT3<double, size_type, 5, 5,    problem_size+1, 1, 1, 0> fjac/*(5, 5, problem_size+1)*/;
StaticArrayT3<double, size_type, 5, 5,    problem_size+1, 1, 1, 0> njac/*(5, 5, problem_size+1)*/;
StaticArrayT4<double, size_type, 5, 5, 3, problem_size+1, 1, 1, 1, 0> lhs /*(5, 5, 3, problem_size+1)*/;
StaticArrayT2<double, size_type, 5, problem_size+1, 1, 0> rtmp/*(5, problem_size+1)*/;
//thread_local NArray<double, 2, size_t, NARRAY_MEMORDER_COL, 1, 1> ce(5, 13);
// common /work_lhs/ fjac, njac, lhs, rtmp, tmp1, tmp2, tmp3

/*
      double precision fjac(5, 5, 0:problem_size),
     >                 njac(5, 5, 0:problem_size),
     >                 lhs (5, 5, 3, 0:problem_size),
     >                 rtmp(5, 0:problem_size),
     >                 tmp1, tmp2, tmp3
      common /work_lhs/ fjac, njac, lhs, rtmp, tmp1, tmp2, tmp3
!$OMP THREADPRIVATE(/work_lhs/)
*/
