#include "mnblas.h"

double absd (double x);
float absf (float x);


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
  float *Xp = (float *)X;
  for (; (i < N) ; i += incX)
    {
      res += absf(Xp [2 * i]) + absf(Xp [2 * i + 1]);
    }

  return res;
}

double mnblas_dzasum(const int N, const void *X, const int incX) {
  register unsigned int i = 0 ;

  double res = 0;
  double * Xp = (double *) X;
  for (; (i < N) ; i += incX)
    {
      res += absd(Xp [2 * i]) + absd(Xp [2 * i + 1]);
    }

  return res;
}



float absf (float x) {
  if (x < 0) return -1 * x;
  return x;
}
double absd (double x) {
  if (x < 0) return -1 * x;
  return x;
}
