#include "mnblas.h"
#include "utils.h"
#include "complexe2.h"

float  mnblas_sasum(const int N, const float *X, const int incX) {
  register unsigned int i = 0 ;

  float res = 0;
  for (; (i < N) ; i += incX)
    {
      res += X [i] ;
    }

  return res;
}

double mnblas_dasum(const int N, const double *X, const int incX) {
  register unsigned int i = 0 ;

  double res = 0;
  for (; (i < N) ; i += incX)
    {
      res += X [i] ;
    }

  return res;
}

float  mnblas_scasum(const int N, const void *X, const int incX) {
  register unsigned int i = 0 ;

  float res = 0;
  complexe_float_t *Xp = (complexe_float_t *)X;
  for (; (i < N) ; i += incX)
    {
      res += absf(Xp[i].real) + absf(Xp[i].imaginary);
    }

  return res;
}

double mnblas_dzasum(const int N, const void *X, const int incX) {
  register unsigned int i = 0 ;

  double res = 0;
  complexe_double_t *Xp = (complexe_double_t *) X;
  for (; (i < N) ; i += incX)
    {
      res += absd(Xp[i].real) + absd(Xp[i].imaginary);
    }
  return res;
}
