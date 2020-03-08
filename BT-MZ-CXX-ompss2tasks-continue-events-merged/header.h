#ifndef HAVE_HEADER_H
#define HAVE_HEADER_H

/*
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  header.h
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none

c---------------------------------------------------------------------
c The following include file is generated automatically by the
c "setparams" utility. It defines
c      problem_size:  maximum overall grid size
c      dt_default:    default time step for this problem size if no
c                     config file
c      niter_default: default number of iterations for this problem size
c---------------------------------------------------------------------
*/

#include <cmath>
#include <mpi.h>

#include "npbparams_cc.h"
#include "NArray.h"
#include "../common/timers.h"

constexpr const int aa = 1, bb = 2, cc = 3, BLOCK_SIZE = 5;

extern int     npb_verbose;
extern double  elapsed_time;
extern int     timeron;

extern double  tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3,
               dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4,
               dy5, dz1, dz2, dz3, dz4, dz5, dssp, dt,
               dxmax, dymax, dzmax, xxcon1, xxcon2,
               xxcon3, xxcon4, xxcon5, dx1tx1, dx2tx1, dx3tx1,
               dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4,
               yycon5, dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
               zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1,
               dz2tz1, dz3tz1, dz4tz1, dz5tz1, dnxm1, dnym1,
               dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2,
               c3, c4, c5, c4dssp, c5dssp, dtdssp, dttx1,
               dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1,
               c2dtty1, c2dttz1, comz1, comz4, comz5, comz6,
               c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16;

extern
StaticArrayT2<double, size_type, 5, 13, 1, 1> ce/*(5,13)*/;

constexpr const int  max_zones = x_zones*y_zones;

extern
StaticArrayT1<size_type, size_type, x_zones, 1> x_start;
extern
StaticArrayT1<size_type, size_type, x_zones, 1> x_end;
extern
StaticArrayT1<size_type, size_type, x_zones, 1> x_size;

extern
StaticArrayT1<size_type, size_type, y_zones, 1> y_start;
extern
StaticArrayT1<size_type, size_type, y_zones, 1> y_end;
extern
StaticArrayT1<size_type, size_type, y_zones, 1> y_size;
extern
StaticArrayT1<size_type, size_type, max_zones, 1> iz_west;
extern
StaticArrayT1<size_type, size_type, max_zones, 1> iz_east;
extern
StaticArrayT1<size_type, size_type, max_zones, 1> iz_south;

extern
StaticArrayT1<size_type, size_type, max_zones, 1> iz_north;

extern
StaticArrayT1<size_type, size_type, max_zones, 1> start1;

extern
StaticArrayT1<size_type, size_type, max_zones, 1> start5;

extern
StaticArrayT1<size_type, size_type, max_zones, 1> qstart_west;

extern
StaticArrayT1<size_type, size_type, max_zones, 1> qstart_east;

extern
StaticArrayT1<size_type, size_type, max_zones, 1> qstart_south;

extern
StaticArrayT1<size_type, size_type, max_zones, 1> qstart_north;


/*
c-----------------------------------------------------------------------
c   Timer constants
c-----------------------------------------------------------------------
*/
constexpr const int t_total = 1;
constexpr const int t_rhsx  = 2;
constexpr const int t_rhsy = 3;
constexpr const int t_rhsz = 4;
constexpr const int t_rhs = 5;
constexpr const int t_xsolve = 6;
constexpr const int t_ysolve = 7;
constexpr const int t_zsolve = 8;
constexpr const int t_rdis1 = 9;
constexpr const int t_rdis2 = 10;
constexpr const int t_add = 11;
constexpr const int t_last = 11;

/****
 * Some wrappers for Fortran functions throughout the code
 */


/**
 * Absolute value of a double
 */

static constexpr inline
double dabs(double value){
  return std::abs(value);
}


/**
 * Cast to a double
 */
template<typename T>
static constexpr inline
double dble(T value){
  return static_cast<double>(value);
}

/**
 * Natural logarithm
 */
static inline
double dlog(double value){
  return std::log(value);
}

/**
 * exp()
 */
static inline
double dexp(double value){
  return std::exp(value);
}

/**
 * sqrt()
 */
static inline
double dsqrt(double value){
  return std::sqrt(value);
}

/**
 * Cast to a double
 */
template<typename T, typename S>
static constexpr inline
T mod(T value, S modulo){
  return value % modulo;
}

/**
 * max double
 */
static constexpr inline
double dmax1(double a, double b){
  return std::max(a, b);
}

/* Functions */
void set_constants();


template<typename T1, typename T2>
using doloop_t = decltype(std::declval<T1>()+std::declval<T2>());

/**
 * For loop abstraction similar to a fortran do loop:
 *
 * doloop(a, b, [](auto i){ std::cout << i; });
 *
 * will print all numbers in the *inclusive* range [a, b].
 */
template<typename UnaryFuncT, typename Size1T, typename Size2T>
void doloop(Size1T from, Size2T to, UnaryFuncT fn)
{
  using SizeT = decltype(from+to);
  for (SizeT i = from; i <= to; ++i) {
    fn(i);
  }
}

#define DO(var, from, to, body) \
  doloop(from, to, [&](doloop_t<decltype(from), decltype(to)> var) body )


void exact_rhs(
  double* forcing,
  size_type nx,
  size_type nxmax,
  size_type ny,
  size_type nz);


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
  int iprn_msg);

#endif // HAVE_HEADER_H
