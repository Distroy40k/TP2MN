#include "mnblas.h"
#include "utils.h"

CBLAS_INDEX mnblas_isamax(const int N, const float  *X, const int incX) {
  if ((N < 0) || (incX < 0)) return 0;
  register unsigned int i = 0 ;
  CBLAS_INDEX res = i;
  for (; (i < N) ; i += incX) {
    if (X[i] > X[res]) {
      res = (CBLAS_INDEX) i;
    }
  }
  return res;
}

CBLAS_INDEX mnblas_idamax(const int N, const double *X, const int incX) {
  if ((N < 0) || (incX < 0)) return 0;

  register unsigned int i = 0 ;
  CBLAS_INDEX res = i;
  for (; (i < N) ; i += incX) {
    if (X[i] > X[res]) {
      res = (CBLAS_INDEX) i;
    }
  }
  return res;
}

CBLAS_INDEX mnblas_icamax(const int N, const void   *X, const int incX) {
  if ((N < 0) || (incX < 0)) return 0;
  register unsigned int i = 0 ;
  CBLAS_INDEX res = i;
  float * Xp = (float *)Xp;
  float val = absf(Xp[2 * res]) + absf(Xp[2 * res + 1]);
  float tmp_val;
  for (; (i < N) ; i += incX) {
    tmp_val = absf(Xp[2 * i]) + absf(Xp[2 * i + 1]);
    if (tmp_val > val) {
      res = (CBLAS_INDEX) i;
      val = tmp_val;
    }
  }
  return res;
}

CBLAS_INDEX mnblas_izamax(const int N, const void   *X, const int incX) {
  if ((N < 0) || (incX < 0)) return 0;
  register unsigned int i = 0 ;
  CBLAS_INDEX res = i;
  double * Xp = (double *)X;
  double val = absd(Xp[2 * res]) + absd(Xp[2 * res + 1]);
  double tmp_val;
  for (; (i < N) ; i += incX) {
    tmp_val = absd(Xp[2 * i]) + absd(Xp[2 * i + 1]);
    if (tmp_val > val) {
      res = (CBLAS_INDEX) i;
      val = tmp_val;
    }
  }
  return res;
}
