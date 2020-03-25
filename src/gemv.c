#include "mnblas.h"
#include <stdio.h>
#include "complexe2.h"
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

    register unsigned int i = 0;
    register unsigned int j = 0;
    register float beta_alpha = beta/alpha;
    register unsigned int in;
    for (; i < M; i += incX) {
      Y[i] *= beta_alpha;
      in = i * N;
      for (j = 0; j < N; j += incX) {
        Y[i] += A[in + j] * X[j];
      }
      Y[i] *= alpha;
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A,
  const int lda,const double *X, const int incX, const double beta, double *Y, const int incY) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    register double beta_alpha = beta/alpha;
    register unsigned int in;

    for (; i < M; i += incX) {
      Y[i] *= beta_alpha;
      in = i * N;
      for (j = 0; j < N; j += incX) {
        Y[i] += A[in + j] * X[j];
      }
      Y[i] *= alpha;
    }

  }

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha,
  const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    complexe_float_t *Ap = (complexe_float_t *) A;
    complexe_float_t *Xp = (complexe_float_t *) X;
    complexe_float_t *Yp = (complexe_float_t *) Y;
    complexe_float_t *alphap = (complexe_float_t *) alpha;
    complexe_float_t *betap = (complexe_float_t *) beta;
    register complexe_float_t beta_alpha = div_complexe_float(*betap, *alphap);

    for (; i < M; i += incX) {
      Yp[i] = mult_complexe_float(Yp[i], beta_alpha);
      for (j = 0; j < N; j += incX) {
        Yp[i] = add_complexe_float(mult_complexe_float(Xp[i], Ap[i]), Yp[i]);
      }
      Yp[i] = mult_complexe_float(Yp[i], *alphap);
    }
  }

void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha, const void *A,
  const int lda, const void *X, const int incX, const void *beta,void *Y, const int incY) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    complexe_double_t *Ap = (complexe_double_t *) A;
    complexe_double_t *Xp = (complexe_double_t *) X;
    complexe_double_t *Yp = (complexe_double_t *) Y;
    complexe_double_t *alphap = (complexe_double_t *) alpha;
    complexe_double_t *betap = (complexe_double_t *) beta;
    register complexe_double_t beta_alpha = div_complexe_double(*betap, *alphap);

    for (; i < M; i += incX) {
      Yp[i] = mult_complexe_double(Yp[i], beta_alpha);
      for (j = 0; j < N; j += incX) {
        Yp[i] = add_complexe_double(mult_complexe_double(Xp[i], Ap[i]), Yp[i]);
      }
      Yp[i] = mult_complexe_double(Yp[i], *alphap);
    }
   }
