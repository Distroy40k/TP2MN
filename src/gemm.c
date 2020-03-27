#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const float alpha,
                   const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)
{

  register unsigned int k;
  register unsigned int i;
  register unsigned int j;
  register unsigned int K_M;
  register float beta_div_alpha = (beta / alpha);

#pragma omp parallel for private(i, j, K_M)
  for (k = 0; k < M; ++k)
  {
    K_M = k * M;
    for (i = 0; i < M; ++i)
    {
      C[i + K_M] *= beta_div_alpha;

      for (j = 0; j < N; j += 8)
      {
        C[i + K_M] += A[j + K_M] * B[i + j * M];
        C[i + K_M] += A[(j + 1) + K_M] * B[i + (j + 1) * M];
        C[i + K_M] += A[(j + 2) + K_M] * B[i + (j + 2) * M];
        C[i + K_M] += A[(j + 3) + K_M] * B[i + (j + 3) * M];
        C[i + K_M] += A[(j + 4) + K_M] * B[i + (j + 4) * M];
        C[i + K_M] += A[(j + 5) + K_M] * B[i + (j + 5) * M];
        C[i + K_M] += A[(j + 6) + K_M] * B[i + (j + 6) * M];
        C[i + K_M] += A[(j + 7) + K_M] * B[i + (j + 7) * M];
      }

      C[i + K_M] *= alpha;
    }
  }
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha,
                   const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{

  register unsigned int k;
  register unsigned int i;
  register unsigned int j;
  register unsigned int K_M;
  register double beta_div_alpha = (beta / alpha);

#pragma omp parallel for private(i, j, K_M)
  for (k = 0; k < M; ++k)
  {
    K_M = k * M;

    for (i = 0; i < M; ++i)
    {
      C[i + K_M] *= beta_div_alpha;

      for (j = 0; j < N; j += 8)
      {
        C[i + K_M] += A[j + K_M] * B[i + j * M];
        C[i + K_M] += A[(j + 1) + K_M] * B[i + (j + 1) * M];
        C[i + K_M] += A[(j + 2) + K_M] * B[i + (j + 2) * M];
        C[i + K_M] += A[(j + 3) + K_M] * B[i + (j + 3) * M];
        C[i + K_M] += A[(j + 4) + K_M] * B[i + (j + 4) * M];
        C[i + K_M] += A[(j + 5) + K_M] * B[i + (j + 5) * M];
        C[i + K_M] += A[(j + 6) + K_M] * B[i + (j + 6) * M];
        C[i + K_M] += A[(j + 7) + K_M] * B[i + (j + 7) * M];
      }

      C[i + K_M] *= alpha;
    }
  }
}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void *alpha,
                   const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc)
{

  complexe_float_t *A_tab = (complexe_float_t *)A;
  complexe_float_t *B_tab = (complexe_float_t *)B;
  complexe_float_t *C_tab = (complexe_float_t *)C;

  register unsigned int k;
  register unsigned int i;
  register unsigned int j;
  register unsigned int K_M;
  register unsigned int I_K_M;

  complexe_float_t *alpha_c = (complexe_float_t *)alpha;
  complexe_float_t *beta_c = (complexe_float_t *)beta;
  complexe_float_t beta_div_alpha = div_complexe_float(*beta_c, *alpha_c);

#pragma omp parallel for private(i, j, K_M)
  for (k = 0; k < M; ++k)
  {
    K_M = k * M;
    for (i = 0; i < M; ++i)
    {
      I_K_M = i + K_M;
      C_tab[I_K_M] = mult_complexe_float(C_tab[I_K_M], beta_div_alpha);
      for (j = 0; j < N; j += 8)
      {
        C_tab[I_K_M] = add_complexe_float(C_tab[I_K_M], mult_complexe_float(A_tab[j + K_M], B_tab[i + j * M]));
        C_tab[I_K_M] = add_complexe_float(C_tab[I_K_M], mult_complexe_float(A_tab[(j + 1) + K_M], B_tab[i + (j + 1) * M]));
        C_tab[I_K_M] = add_complexe_float(C_tab[I_K_M], mult_complexe_float(A_tab[(j + 2) + K_M], B_tab[i + (j + 2) * M]));
        C_tab[I_K_M] = add_complexe_float(C_tab[I_K_M], mult_complexe_float(A_tab[(j + 3) + K_M], B_tab[i + (j + 3) * M]));
        C_tab[I_K_M] = add_complexe_float(C_tab[I_K_M], mult_complexe_float(A_tab[(j + 4) + K_M], B_tab[i + (j + 4) * M]));
        C_tab[I_K_M] = add_complexe_float(C_tab[I_K_M], mult_complexe_float(A_tab[(j + 5) + K_M], B_tab[i + (j + 5) * M]));
        C_tab[I_K_M] = add_complexe_float(C_tab[I_K_M], mult_complexe_float(A_tab[(j + 6) + K_M], B_tab[i + (j + 6) * M]));
        C_tab[I_K_M] = add_complexe_float(C_tab[I_K_M], mult_complexe_float(A_tab[(j + 7) + K_M], B_tab[i + (j + 7) * M]));
      }
      C_tab[I_K_M] = mult_complexe_float(C_tab[I_K_M], *alpha_c);
    }
  }
}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void *alpha,
                   const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc)
{

  complexe_double_t *A_tab = (complexe_double_t *)A;
  complexe_double_t *B_tab = (complexe_double_t *)B;
  complexe_double_t *C_tab = (complexe_double_t *)C;

  register unsigned int k;
  register unsigned int i;
  register unsigned int j;
  register unsigned int K_M;
  register unsigned int I_K_M;

  complexe_double_t *alpha_c = (complexe_double_t *)alpha;
  complexe_double_t *beta_c = (complexe_double_t *)beta;
  complexe_double_t beta_div_alpha = div_complexe_double(*beta_c, *alpha_c);

#pragma omp parallel for private(i, j, K_M)
  for (k = 0; k < M; ++k)
  {
    K_M = k * M;
    for (i = 0; i < M; ++i)
    {
      I_K_M = i + K_M;
      C_tab[I_K_M] = mult_complexe_double(C_tab[I_K_M], beta_div_alpha);

      for (j = 0; j < N; j += 8)
      {
        // C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab [j + K_M], B_tab[i + j*M]));
        C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab[j + K_M], B_tab[i + j * M]));
        C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab[(j + 1) + K_M], B_tab[i + (j + 1) * M]));
        C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab[(j + 2) + K_M], B_tab[i + (j + 2) * M]));
        C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab[(j + 3) + K_M], B_tab[i + (j + 3) * M]));
        C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab[(j + 4) + K_M], B_tab[i + (j + 4) * M]));
        C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab[(j + 5) + K_M], B_tab[i + (j + 5) * M]));
        C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab[(j + 6) + K_M], B_tab[i + (j + 6) * M]));
        C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab[(j + 7) + K_M], B_tab[i + (j + 7) * M]));
      }
      C_tab[I_K_M] = mult_complexe_double(C_tab[I_K_M], *alpha_c);
    }
  }
}
