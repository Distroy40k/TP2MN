#include "mnblas.h"
#include "complexe2.h"

void mncblas_sswap(const int N, float *X, const int incX,
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;

  for (; ((i < N) && (j < N)) ; i += 1, j+=1) // De base i += incX & j+= incY, modif pour les performances
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

  for (; ((i < N) && (j < N)) ; i += 1, j+=1) // De base i += incX & j+= incY, modif pour les performances
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
  register complexe_float_t save ;
  complexe_float_t * Xp = (complexe_float_t *) X;
  complexe_float_t * Yp = (complexe_float_t *) Y;

  for (; ((i < N) && (j < N)) ; i += 1, j+=1) // De base i += incX & j+= incY, modif pour les performances
    {
      save = Yp[j];
      Yp[j] = Xp[i];
      Xp[i] = save;
    }
  return ;
}

void mncblas_zswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register complexe_double_t save ;
  complexe_double_t * Xp = (complexe_double_t *) X;
  complexe_double_t * Yp = (complexe_double_t *) Y;

  for (; ((i < N) && (j < N)) ; i += 1, j+=1) // De base i += incX & j+= incY, modif pour les performances
    {
      save = Yp[j];
      Yp[j] = Xp[i];
      Xp[i] = save;
    }
  return ;
}
