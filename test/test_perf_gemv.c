#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"
#include "flop.h"

#define N 128
#define NB_FOIS 10

void fvector_init(float *V, float x) {
  for (int i = 0; i < N; i++) {
    V[i] = x;
  }
}
void dvector_init(double *V, double x) {
  for (int i = 0; i < N; i++) {
    V[i] = x;
  }
}
void fmatrice_init( float *A, float x) {
  for (int i = 0; i < N * N; i++) {
    A[i] = x;
  }
}

void dmatrice_init(double *A,double x) {
  for (int i = 0; i < N * N; i++) {
    A[i] = x;
  }
}
void fmatriceC_init(complexe_float_t *A, complexe_float_t x) {
  for (int i = 0; i < N * N; i++) {
    A[i] = x;
  }
}
void dmatriceC_init(complexe_double_t *A, complexe_double_t x) {
  for (int i = 0; i < N * N; i++) {
    A[i] = x;
  }
}

void fvectorC_init(complexe_float_t *V, complexe_float_t x) {
  for (int i = 0; i < N ; i++) {
    V[i] = x;
  }
}

void dvectorC_init(complexe_double_t *V, complexe_double_t x) {
  for (int i = 0; i < N; i++) {
    V[i] = x;
  }
}


int main(int argc, char **argv) {
  int i;
  float fA[N * N];
  double dA[N * N];
  float fX[N];
  double dX[N];
  float fY[N];
  double dY[N];

  complexe_float_t CfA[N * N];
  complexe_double_t CdA[N * N];
  complexe_float_t CfX[N];
  complexe_double_t CdX[N];
  complexe_float_t CfY[N];
  complexe_double_t CdY[N];

  float falpha = (float)1;
  float fbeta = (float)1;
  double dalpha = (double)1;
  double dbeta = (double)1;
  complexe_float_t Cfalpha = init_complexe_float((float) 1, (float) 1);
  complexe_double_t Cdalpha = init_complexe_double((double) 1, (double) 1);
  complexe_float_t Cfbeta = init_complexe_float((float) 1, (float) 1);
  complexe_double_t Cdbeta = init_complexe_double((double) 1, (double) 1);

  unsigned long long int start, end ;

  init_flop () ;
  printf("\nTests des fonctions de gemv\n\n");
  printf("\n##################\nTests de mncblas_sgemv\n##################\n");
  fmatrice_init(fA, (float)1);
  fvector_init(fX, (float)1);
  fvector_init(fY, (float)1);
  mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, falpha, fA, 0, fX, 1, fbeta, fY, 1);
  start = _rdtsc () ;
  for (i = 0; i < NB_FOIS; i ++) {
      mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, falpha, fA, 0, fX, 1, fbeta, fY, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_sgemv", N * (2 * N + 3) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_dgemv\n##################\n");
  dmatrice_init(dA, (double)1);
  dvector_init(dX, (double)1);
  dvector_init(dY, (double)1);
  mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, dalpha, dA, 0, dX, 1, dbeta, dY, 1);
  start = _rdtsc () ;
  for (i = 0; i < NB_FOIS; i ++) {
    mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, dalpha, dA, 0, dX, 1, dbeta, dY, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_dgemv", N * (2 * N +3) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_cgemv\n##################\n");
  fmatriceC_init(CfA, init_complexe_float((float) 1, (float) 1));
  fvectorC_init(CfX, init_complexe_float((float) 1, (float) 1));
  fvectorC_init(CfY, init_complexe_float((float) 1, (float) 1));

  mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, &Cfalpha, CfA, 0, CfX, 1, &Cfbeta, CfY, 1);
  start = _rdtsc () ;
  for (i = 0; i < NB_FOIS; i ++) {
    mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, &Cfalpha, CfA, 0, CfX, 1, &Cfbeta, CfY, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_cgemv", N * (14 + (N * 8)) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_zgemv\n##################\n");
  dmatriceC_init(CdA, init_complexe_double((double) 1, (double) 1));
  dvectorC_init(CdX, init_complexe_double((double) 1, (double) 1));
  dvectorC_init(CdY, init_complexe_double((double) 1, (double) 1));

  mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, &Cdalpha, CdA, 0, CdX, 1, &Cdbeta, CdY, 1);
  start = _rdtsc () ;
  for (i = 0; i < NB_FOIS; i ++) {
    mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, &Cdalpha, CdA, 0, CdX, 1, &Cdbeta, CdY, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_zgemv", N * (14 + (N * 8)) * NB_FOIS, end-start) ;
}
