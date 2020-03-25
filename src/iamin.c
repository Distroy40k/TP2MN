#include "mnblas.h"
#include "utils.h"
#include "complexe2.h"

CBLAS_INDEX mnblas_isamin(const int N, const float  *X, const int incX) {
  if ((N < 0) || (incX < 0)) return 0;
  register unsigned int i = 0 ;
  CBLAS_INDEX res = i;
  for (; (i < N) ; i += incX) {
    if (X[i] < X[res]) {
      res = (CBLAS_INDEX) i;
    }
  }
  return res;
}

CBLAS_INDEX mnblas_idamin(const int N, const double *X, const int incX) {
  if ((N < 0) || (incX < 0)) return 0;

  register unsigned int i = 0 ;
  CBLAS_INDEX res = i;
  for (; (i < N) ; i += incX) {
    if (X[i] < X[res]) {
      res = (CBLAS_INDEX) i;
    }
  }
  return res;
}

CBLAS_INDEX mnblas_icamin(const int N, const void   *X, const int incX) {
  if ((N < 0) || (incX < 0)) return 0;
  register unsigned int i = 0 ;
  CBLAS_INDEX res = i;
  complexe_float_t * Xp = (complexe_float_t *)X;
  float val = absf(Xp[i].real) + absf(Xp[i].imaginary);
  float tmp_val;
  for (; (i < N) ; i += incX) {
    tmp_val = absf(Xp[i].real) + absf(Xp[i].imaginary);
    if (tmp_val < val) {
      res = (CBLAS_INDEX) i;
      val = tmp_val;
    }
  }
  return res;
}

CBLAS_INDEX mnblas_izamin(const int N, const void   *X, const int incX) {
  if ((N < 0) || (incX < 0)) return 0;
  register unsigned int i = 0 ;
  CBLAS_INDEX res = i;
  complexe_double_t * Xp = (complexe_double_t *)X;
  double val = absd(Xp[i].real) + absd(Xp[i].imaginary);
  double tmp_val;
  for (; (i < N) ; i += incX) {
    tmp_val = absd(Xp[i].real) + absd(Xp[i].imaginary);
    if (tmp_val < val) {
      res = (CBLAS_INDEX) i;
      val = tmp_val;
    }
  }
  return res;
}
