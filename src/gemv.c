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

    register unsigned int i;
    register unsigned int j;
    register float beta_alpha = beta/alpha;
    register unsigned int I_N;

    #pragma omp parallel for private(i, j, I_N)
    for (i = 0; i < M; i += 1) { // De base i += incX, modif pour les performances
      Y[i] *= beta_alpha;
      I_N = i * N;
      for (j = 0; j < N; j += 8) {
        Y[i] += A[I_N + j] * X[j];
        Y[i] += A[I_N + j + 1] * X[j + 1];
        Y[i] += A[I_N + j + 2] * X[j + 2];
        Y[i] += A[I_N + j + 3] * X[j + 3];
        Y[i] += A[I_N + j + 4] * X[j + 4];
        Y[i] += A[I_N + j + 5] * X[j + 5];
        Y[i] += A[I_N + j + 6] * X[j + 6];
        Y[i] += A[I_N + j + 7] * X[j + 7];

      }
      Y[i] *= alpha;
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A,
  const int lda,const double *X, const int incX, const double beta, double *Y, const int incY) {
    register unsigned int i;
    register unsigned int j;
    register double beta_alpha = beta/alpha;
    register unsigned int I_N;

    #pragma omp parallel for private(i, j, I_N)
    for (i = 0; i < M; i += 1) { // De base i += incX, modif pour les performances
      Y[i] *= beta_alpha;
      I_N = i * N;
      for (j = 0; j < N; j += 8) {
        Y[i] += A[I_N + j] * X[j];
        Y[i] += A[I_N + j + 1] * X[j + 1];
        Y[i] += A[I_N + j + 2] * X[j + 2];
        Y[i] += A[I_N + j + 3] * X[j + 3];
        Y[i] += A[I_N + j + 4] * X[j + 4];
        Y[i] += A[I_N + j + 5] * X[j + 5];
        Y[i] += A[I_N + j + 6] * X[j + 6];
        Y[i] += A[I_N + j + 7] * X[j + 7];
      }
      Y[i] *= alpha;
    }

  }

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha,
  const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY) {
    register unsigned int i;
    register unsigned int j;
    complexe_float_t *Ap = (complexe_float_t *) A;
    complexe_float_t *Xp = (complexe_float_t *) X;
    complexe_float_t *Yp = (complexe_float_t *) Y;
    complexe_float_t *alphap = (complexe_float_t *) alpha;
    complexe_float_t *betap = (complexe_float_t *) beta;
    register complexe_float_t beta_alpha = div_complexe_float(*betap, *alphap);

    #pragma omp parallel for private(i, j)
    for (i = 0; i < M; i += 1) { // De base i += incX, modif pour les performances
      Yp[i] = mult_complexe_float(Yp[i], beta_alpha);
      for (j = 0; j < N; j += 8) {
        Yp[i] = add_complexe_float(mult_complexe_float(Xp[i], Ap[i]), Yp[i]);
        Yp[i] = add_complexe_float(mult_complexe_float(Xp[i + 1], Ap[i + 1]), Yp[i]);
        Yp[i] = add_complexe_float(mult_complexe_float(Xp[i + 2], Ap[i + 2]), Yp[i]);
        Yp[i] = add_complexe_float(mult_complexe_float(Xp[i + 3], Ap[i + 3]), Yp[i]);
        Yp[i] = add_complexe_float(mult_complexe_float(Xp[i + 4], Ap[i + 4]), Yp[i]);
        Yp[i] = add_complexe_float(mult_complexe_float(Xp[i + 5], Ap[i + 5]), Yp[i]);
        Yp[i] = add_complexe_float(mult_complexe_float(Xp[i + 6], Ap[i + 6]), Yp[i]);
        Yp[i] = add_complexe_float(mult_complexe_float(Xp[i + 7], Ap[i + 7]), Yp[i]);

      }
      Yp[i] = mult_complexe_float(Yp[i], *alphap);
    }
  }

void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha, const void *A,
  const int lda, const void *X, const int incX, const void *beta,void *Y, const int incY) {
    register unsigned int i;
    register unsigned int j;
    complexe_double_t *Ap = (complexe_double_t *) A;
    complexe_double_t *Xp = (complexe_double_t *) X;
    complexe_double_t *Yp = (complexe_double_t *) Y;
    complexe_double_t *alphap = (complexe_double_t *) alpha;
    complexe_double_t *betap = (complexe_double_t *) beta;
    register complexe_double_t beta_alpha = div_complexe_double(*betap, *alphap);

    #pragma omp parallel for private(i, j)

    for (i = 0; i < M; i += 1) { // De base i += incX, modif pour les performances
      Yp[i] = mult_complexe_double(Yp[i], beta_alpha);
      for (j = 0; j < N; j += 8) {
        Yp[i] = add_complexe_double(mult_complexe_double(Xp[j], Ap[j]), Yp[i]);
        Yp[i] = add_complexe_double(mult_complexe_double(Xp[j + 1], Ap[j + 1]), Yp[i]);
        Yp[i] = add_complexe_double(mult_complexe_double(Xp[j + 2], Ap[j + 2]), Yp[i]);
        Yp[i] = add_complexe_double(mult_complexe_double(Xp[j + 3], Ap[j + 3]), Yp[i]);
        Yp[i] = add_complexe_double(mult_complexe_double(Xp[j + 4], Ap[j + 4]), Yp[i]);
        Yp[i] = add_complexe_double(mult_complexe_double(Xp[j + 5], Ap[j + 5]), Yp[i]);
        Yp[i] = add_complexe_double(mult_complexe_double(Xp[j + 6], Ap[j + 6]), Yp[i]);
        Yp[i] = add_complexe_double(mult_complexe_double(Xp[j + 7], Ap[j + 7]), Yp[i]);
      }
      Yp[i] = mult_complexe_double(Yp[i], *alphap);
    }
   }
