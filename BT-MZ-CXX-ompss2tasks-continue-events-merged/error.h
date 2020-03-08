#ifndef HAVE_ERROR_H
#define HAVE_ERROR_H

#include <cmath>

#include "NArray.h"

#include "npbparams_cc.h"

void error_norm_impl(
  double* rms_ptr,
  double* u_ptr,
  double* rms_loc_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);


void rhs_norm_impl(
  double* rms_ptr,
  double* rhs_ptr,
  double* rms_loc_ptr,
  size_type k,
  size_type nx, size_type nxmax, size_type ny, size_type nz);


void error_norm(
  double* rms,
  double* u,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

void rhs_norm(
  double* rms,
  double* rhs,
  size_type nx, size_type nxmax, size_type ny, size_type nz);

#endif // HAVE_ERROR_H
