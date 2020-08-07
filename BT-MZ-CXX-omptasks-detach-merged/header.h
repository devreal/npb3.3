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
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <omp.h>
#include <thread>
#include <chrono>

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

// thread-private global variables
//extern thread_local StaticArrayT3<double, size_type, 5, 5,    problem_size+1, 1, 1, 0> fjac/*(5, 5, problem_size+1)*/;
//extern thread_local StaticArrayT3<double, size_type, 5, 5,    problem_size+1, 1, 1, 0> njac/*(5, 5, problem_size+1)*/;
//extern thread_local StaticArrayT4<double, size_type, 5, 5, 3, problem_size+1, 1, 1, 1, 0> lhs /*(5, 5, 3, problem_size+1)*/;
//extern thread_local StaticArrayT2<double, size_type, 5, problem_size+1, 1, 0> rtmp/*(5, problem_size+1)*/;


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

static inline
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



template<typename ValueT, typename SizeT>
void print_zone(
  ValueT* u_ptr,
  SizeT nx, SizeT nxmax, SizeT ny, SizeT nz, int zone, int step)
{
#if 0
  ArrayViewT4<ValueT, SizeT> u(u_ptr, 5, nxmax, ny, nz);

  std::cout << "------ Printing U zone " << zone << " at step " << step << " -------" << std::fixed << std::endl;
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


bool mpi_poll_service(void* data);

#include <mpi.h>
#include <chrono>

//#define USE_PROGRESS_THREAD

static
void print_binding()
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    printf("Set returned by pthread_getaffinity_np() contained:\n");
    for (int j = 0; j < CPU_SETSIZE; j++)
        if (CPU_ISSET(j, &cpuset))
            printf("  %d", j);
    printf("Running on thread %d\n", omp_get_thread_num());
}

static volatile int do_progress = 1;

static std::thread progress_thread;

static void start_progress_thread()
{
#if defined(USE_PROGRESS_THREAD)
  progress_thread = std::thread([&]()
  {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    for (int i = 0; i < CPU_SETSIZE; ++i) {
      if (CPU_ISSET(i, &cpuset)) {
        CPU_CLR(i, &cpuset);
        break;
      }
    }
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    print_binding();

    while(do_progress) {
      //std::cout << "Calling poll service" << std::endl;
      mpi_poll_service(nullptr);
      std::this_thread::sleep_for(std::chrono::microseconds(100));
    }
  });
#endif // defined(USE_PROGRESS_THREAD)
}

static void stop_progress_thread()
{
#if defined(USE_PROGRESS_THREAD)
  do_progress = 0;
  progress_thread.join();
#endif
}

template<typename UnaryT>
void parallel_master_with_progress(UnaryT fn)
{
  //print_binding();
#pragma omp parallel default(shared)
{
#pragma omp master
{
  volatile int do_progress = 1; // whether to keep triggering progress
  volatile int in_compute  = 0; // whether the compute task has been launched
//#if 0
#if !defined(USE_PROGRESS_THREAD)
  // we cannot safely run this on just one thread
  assert(omp_get_num_threads() > 1);
  #pragma omp task shared(do_progress, in_compute)
  {
    constexpr const auto poll_interval = std::chrono::microseconds{1};
    //using clock = std::chrono::high_resolution_clock;
    using clock = std::chrono::system_clock;
    auto last_progress_ts = clock::now();
    while(do_progress) {
      auto ts = clock::now();
      if (poll_interval < ts - last_progress_ts) 
      {
        //std::cout << "Calling MPI poll service " << std::endl;
        mpi_poll_service(nullptr);
        last_progress_ts = clock::now();
      }
      // it is not safe to yield here because the thread might steal the compute task
      if (in_compute) {
        #pragma omp taskyield
      } else {
        usleep(1);
      }
    }
  }
#endif // USE_PROGRESS_THREAD
//#endif 
#if !defined(USE_PROGRESS_THREAD)
#pragma omp task default(shared)
#endif // USE_PROGRESS_THREAD
{
  in_compute = 1;
  // invoke the function
  fn();

  // wait for all tasks created inside the function
  #pragma omp taskwait

  //MPI_Barrier(MPI_COMM_WORLD);
  do_progress = 0;
}

  #pragma omp taskwait
}
}
}


#endif // HAVE_HEADER_H
