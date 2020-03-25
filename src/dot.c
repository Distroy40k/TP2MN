#include "mnblas.h"
#include <stdio.h>
#include <math.h>
#include "complexe2.h"

/*
float mncblas_sdot(const int N, const float *X, const int incX,
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}
*/

float mncblas_sdot(const int N, const float *X, const int incX,
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  float dot = 0.0 ;


  for (i; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}

double mncblas_ddot(const int N, const double *X, const int incX,
                 const double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  double dot = 0.0 ;


  for (i; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;

}

void   mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  complexe_float_t *res = (complexe_float_t *) dotu;

  complexe_float_t * Xp = (complexe_float_t *) X;
  complexe_float_t * Yp = (complexe_float_t *) Y;

  for (; i < N ; i += incX) {
    *res = add_complexe_float(*res,mult_complexe_float(Xp[i],Yp[j]));
    j+= incY;
  }

}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  complexe_float_t *res = (complexe_float_t *) dotc;

  complexe_float_t * Xp = (complexe_float_t *) X;
  complexe_float_t * Yp = (complexe_float_t *) Y;

  for (; i < N ; i += incX) {
    *res = add_complexe_float(*res,mult_complexe_float(conj_float(Xp[i]),Yp[j]));
    j+= incY;
  }
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  complexe_double_t *res = (complexe_double_t *) dotu;

  complexe_double_t * Xp = (complexe_double_t *) X;
  complexe_double_t * Yp = (complexe_double_t *) Y;

  for (; i < N ; i += incX) {
    *res = add_complexe_double(*res,mult_complexe_double(Xp[i],Yp[j]));
    j+= incY;
  }

}

void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  complexe_double_t *res = (complexe_double_t *) dotc;

  complexe_double_t * Xp = (complexe_double_t *) X;
  complexe_double_t * Yp = (complexe_double_t *) Y;

  for (; i < N ; i += incX) {
    *res = add_complexe_double(*res,mult_complexe_double(conj_double(Xp[i]),Yp[j]));
    j+= incY;
  }
}
