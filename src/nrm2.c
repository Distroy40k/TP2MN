#include "mnblas.h"
#include "utils.h"

#include <math.h>

float  mnblas_snrm2(const int N, const float *X, const int incX){
    register unsigned int i = 0 ;
    float res = 0;
    for (; (i < N) ; i += incX) {
      res += mn_squaref(X[i]);
    }
    return sqrt(res);
}

double mnblas_dnrm2(const int N, const double *X, const int incX) {
  register unsigned int i = 0 ;
  double res = 0;
  for (; (i < N) ; i += incX) {
    res += mn_squared(X[i]);
  }
  return sqrt(res);
}

float  mnblas_scnrm2(const int N, const void *X, const int incX) {
  register unsigned int i = 0 ;
  float * Xp = (float *)X;
  float res = 0;
  for (; (i < N) ; i += incX) {
    res += mn_squaref(Xp[2 * i]) + mn_squaref(Xp[2 * i + 1]);
  }
  return sqrt(res);
}

double mnblas_dznrm2(const int N, const void *X, const int incX) {
  register unsigned int i = 0 ;
  double * Xp = (double *)X;
  double res = 0;
  for (; (i < N) ; i += incX) {
    res += mn_squaref(Xp[2 * i]) + mn_squaref(Xp[2 * i + 1]);
  }
  return sqrt(res);

}
