#include "mnblas.h"
#include <stdio.h>

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


  for (i = 0 ; i < N ; i += incX)
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


  for (i = 0 ; i < N ; i += incX)
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

  *((float *) dotu) = 0;
  *(((float *) dotu) + 1) = 0;

  for (i = 0 ; i < N ; i += incX) {

    *((float *) dotu) += (((float *) X) [2 * i] * ((float *) Y) [2 * j]) - (((float *) X) [2 * i + 1] * ((float *) Y) [2 * j + 1]);
    *(((float *) dotu) + 1) += (((float *) X) [2 * i + 1] * ((float *) Y) [2 * j + 1]) - (((float *) X) [2 * i + 1] * ((float *) Y) [2 * j]);
    j+= incY;
  }

}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  *((float *) dotc) = 0;
  *(((float *) dotc) + 1) = 0;

  for (i = 0 ; i < N ; i += incX) {

    *((float *) dotc) += (((float *) X) [2 * i] * ((float *) Y) [2 * j]) + (((float *) X) [2 * i + 1] * ((float *) Y) [2 * j + 1]);
    *(((float *) dotc) + 1) += (-1 * ((float *) X) [2 * i + 1] * ((float *) Y) [2 * j + 1]) - (-1 * ((float *) X) [2 * i + 1] * ((float *) Y) [2 * j]);
    j+= incY;
  }
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  *((double *) dotu) = 0;
  *(((double *) dotu) + 1) = 0;

  for (i = 0 ; i < N ; i += incX) {

    *((double *) dotu) += (((double *) X) [2 * i] * ((double *) Y) [2 * j]) - (((double *) X) [2 * i + 1] * ((double *) Y) [2 * j + 1]);
    *(((double *) dotu) + 1) += (((double *) X) [2 * i + 1] * ((double *) Y) [2 * j + 1]) - (((double *) X) [2 * i + 1] * ((double *) Y) [2 * j]);
    j+= incY;
  }
}

void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  *((double *) dotc) = 0;
  *(((double *) dotc) + 1) = 0;

  for (i = 0 ; i < N ; i += incX) {

    *((double *) dotc) += (((double *) X) [2 * i] * ((double *) Y) [2 * j]) + (((double *) X) [2 * i + 1] * ((double *) Y) [2 * j + 1]);
    *(((double *) dotc) + 1) += (-1 * ((double *) X) [2 * i + 1] * ((double *) Y) [2 * j + 1]) - (-1 * ((double *) X) [2 * i + 1] * ((double *) Y) [2 * j]);
    j+= incY;
  }

  return ;
}
