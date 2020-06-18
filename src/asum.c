#include "mnblas.h"
#include "utils.h"
#include "complexe2.h"

float  mnblas_sasum(const int N, const float *X, const int incX) {
  register unsigned int i;
  float res = 0;
  #pragma omp parallel for private(i) reduction(+:res)
  for (i = 0; i < N ; i += 1) // De base i += incX, modif pour les performances
    {
      res += X [i] ;
    }

  return res;
}

double mnblas_dasum(const int N, const double *X, const int incX) {
  register unsigned int i ;
  double res = 0;
  #pragma omp parallel for private(i) reduction(+:res)
  for (i = 0; i < N ; i += 1) // De base i += incX, modif pour les performances
    {
      res += X [i] ;
    }

  return res;
}

float  mnblas_scasum(const int N, const void *X, const int incX) {

  register unsigned int i ;
  double res = 0;
  complexe_float_t *Xp = (complexe_float_t *)X;
  #pragma omp parallel for private(i) reduction(+:res)
  for (i = 0; i < N ; i += 1) // De base i += incX, modif pour les performances
    {
      res += absf(Xp[i].real) + absf(Xp[i].imaginary);
    }
  return res;
}

double mnblas_dzasum(const int N, const void *X, const int incX) {
  register unsigned int i = 0 ;
  double res = 0;
  complexe_double_t *Xp = (complexe_double_t *) X;
  #pragma omp parallel for private(i) reduction(+:res)
  for (i = 0; i < N ; i += 1) // De base i += incX, modif pour les performances
    {
      res += absd(Xp[i].real) + absd(Xp[i].imaginary);
    }
  return res;
}
