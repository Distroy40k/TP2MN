#include "mnblas.h"
#include "utils.h"
#include "complexe2.h"

#include <math.h>

float  mnblas_snrm2(const int N, const float *X, const int incX){
    register unsigned int i = 0 ;
    float res = 0;
    for (; (i < N) ; i += 1) { // De base i += incX, modif pour les performances
      res += mn_squaref(X[i]);
    }
    return sqrt(res);
}

double mnblas_dnrm2(const int N, const double *X, const int incX) {
  register unsigned int i = 0 ;
  double res = 0;
  for (; (i < N) ; i += 1) { // De base i += incX, modif pour les performances
    res += mn_squared(X[i]);
  }
  return sqrt(res);
}

float  mnblas_scnrm2(const int N, const void *X, const int incX) {
  register unsigned int i = 0 ;
  complexe_float_t * Xp = (complexe_float_t *)X;
  float res = 0;
  for (; (i < N) ; i += 1) { // De base i += incX, modif pour les performances
    res += mn_squaref(Xp[i].real) + mn_squaref(Xp[i].imaginary);
  }
  return sqrt(res);
}

double mnblas_dznrm2(const int N, const void *X, const int incX) {
  register unsigned int i = 0 ;
  complexe_double_t * Xp = (complexe_double_t *)X;
  double res = 0;
  for (; (i < N) ; i += 1) { // De base i += incX, modif pour les performances
    res += mn_squaref(Xp[i].real) + mn_squaref(Xp[i].imaginary);
  }
  return sqrt(res);

}
