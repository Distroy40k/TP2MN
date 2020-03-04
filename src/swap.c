#include "mnblas.h"

void mncblas_sswap(const int N, float *X, const int incX,
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX,
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_cswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = ((float *)Y) [2*j] ;
      ((float *)Y) [2*j] = ((float *)X) [2*i] ;
      ((float *)X) [2*i] = save ;
      save = ((float *)Y) [2*j+1] ;
      ((float *)Y) [2*j+1] = ((float *)X) [2*i+1] ;
      ((float *)X) [2*i+1] = save ;
    }
  return ;
}

void mncblas_zswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = ((double *)Y) [2*j] ;
      ((double *)Y) [2*j] = ((double *)X) [2*i] ;
      ((double *)X) [2*i] = save ;
      save = ((double *)Y) [2*j+1] ;
      ((double *)Y) [2*j+1] = ((double *)X) [2*i+1] ;
      ((double *)X) [2*i+1] = save ;
    }
  return ;
}
