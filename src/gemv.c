#include "mnblas.h"
#include <stdio.h>
/*
Pour simpplifier on utilise les paramètre suivants :

Layout = MNCblasRowMajor
trans = MNCblasNoTrans
m = n

On ne fais pas cas de lda


Le résultat est donc :
y := alpha*A*x + beta*y

*/

void mncblas_sgemv(const MNCBLAS_LAYOUT layout, const MNCBLAS_TRANSPOSE TransA, const int M, const int N, const float alpha,
  const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY) {

    register int i = 0;
    register int j = 0;
    for (; i < M; i += incX) {
      Y[i] *= (beta/alpha);
      for (j = 0; j < N; j += incX) {
        Y[i] += A[N * i + j] * X[j];
      }
      Y[i] *= alpha;
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A,
  const int lda,const double *X, const int incX, const double beta, double *Y, const int incY) {
    register int i = 0;
    register int j = 0;

    for (; i < M; i += incX) {
      Y[i] *= (beta/alpha);
      for (j = 0; j < N; j += incX) {
        Y[i] += A[N * i + j] * X[j];
      }
      Y[i] *= alpha;
    }

  }

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha,
  const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY) {
    register int i = 0;
    register int j = 0;

    float *Ap = (float *) A;
    float *Xp = (float *) X;
    float *Yp = (float *) Y;
    float *alphap = (float *) alpha;
    float *betap = (float *) beta;

    float tmp_y[2];
    float tmp_y2;
    for (; i < M; i += incX) {
      tmp_y[0] =  Yp[2 * i] * betap[0] - Yp[2 * i + 1] * betap[1];
      tmp_y[1] = Yp[2 * i] * betap[1] + Yp[2 * i + 1] * betap[0];
      Yp[2 * i]=0;
      Yp[2*i+1]=0;
      for (j = 0; j < N; j += incX) {
        Yp[2 * i] += Ap[N * 2 * i + 2 * j] * Xp[2 * j] - Ap[N * 2 * i + 2 * j + 1] * Xp[2 * j + 1];
        Yp[2 * i + 1] += Ap[N * 2 * i + 2 * j] * Xp[2 * j + 1] + Ap[N * 2 * i + 2 * j+ 1] * Xp[2 * j];
      }
      tmp_y2 = Yp[2 *i];
      Yp[2*i] = Yp[2*i] * alphap[0] - Yp[2*i+1] * alphap[1] + tmp_y[0];
      Yp[2*i+1] = tmp_y2 * alphap[1] + Yp[2*i+1] * alphap[0] + tmp_y[1];
    }
  }

void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha, const void *A,
  const int lda, const void *X, const int incX, const void *beta,void *Y, const int incY) {
   register int i = 0;
   register int j = 0;

   double *Ap = (double *) A;
   double *Xp = (double *) X;
   double *Yp = (double *) Y;
   double *alphap = (double *) alpha;
   double *betap = (double *) beta;

   double tmp_y[2];
   double tmp_y2;
   for (; i < M; i += incX) {
     tmp_y[0] =  Yp[2 * i] * betap[0] - Yp[2 * i + 1] * betap[1];
     tmp_y[1] = Yp[2 * i] * betap[1] + Yp[2 * i + 1] * betap[0];
     Yp[2 * i]=0;
     Yp[2*i+1]=0;
     for (j = 0; j < N; j += incX) {
       Yp[2 * i] += Ap[N * 2 * i + 2 * j] * Xp[2 * j] - Ap[N * 2 * i + 2 * j + 1] * Xp[2 * j + 1];
       Yp[2 * i + 1] += Ap[N * 2 * i + 2 * j] * Xp[2 * j + 1] + Ap[N * 2 * i + 2 * j+ 1] * Xp[2 * j];
     }
     tmp_y2 = Yp[2 *i];
     Yp[2*i] = Yp[2*i] * alphap[0] - Yp[2*i+1] * alphap[1] + tmp_y[0];
     Yp[2*i+1] = tmp_y2 * alphap[1] + Yp[2*i+1] * alphap[0] + tmp_y[1];
   }
 }
