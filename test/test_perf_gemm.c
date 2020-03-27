#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "flop.h"
#include "complexe2.h"

#define N  48
#define NB_FOIS 1024


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

int main(int argc, char **argv) {

  float fA[N * N];
  double dA[N * N];
  float fB[N * N];
  double dB[N * N];
  float fC[N * N];
  double dC[N * N];

  complexe_float_t CfA[N * N];
  complexe_double_t CdA[N * N];
  complexe_float_t CfB[N * N];
  complexe_double_t CdB[N * N];
  complexe_float_t CfC[N * N];
  complexe_double_t CdC[N * N];

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
  printf("\nTests des fonctions de gemm\n\n");

  printf("\n##################\nTests de mncblas_sgemm\n##################\n");
  fmatrice_init(fA, (float)1);
  fmatrice_init(fB, (float)1);
  fmatrice_init(fC, (float)1);
  mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, falpha, fA, 0, fB, 1, fbeta, fC, 1);
  start = _rdtsc () ;
  for (int k=0; k<NB_FOIS; k++){
    mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, falpha, fA, 0, fB, 1, fbeta, fC, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_sgemm", N * N * ((N * 2) + 3 ) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_dgemm\n##################\n");
  dmatrice_init(dA, (double)1);
  dmatrice_init(dB, (double)1);
  dmatrice_init(dC, (double)1);
  mncblas_dgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, dalpha, dA, 0, dB, 1, dbeta, dC, 1);
  start = _rdtsc () ;
  for (int k=0; k<NB_FOIS; k++){
    mncblas_dgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, dalpha, dA, 0, dB, 1, dbeta, dC, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_dgemm", N * N * ((N * 2) + 3 ) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_cgemm\n##################\n");
  fmatriceC_init(CfA, init_complexe_float((float)1,(float)1));
  fmatriceC_init(CfB, init_complexe_float((float)1,(float)1));
  fmatriceC_init(CfC, init_complexe_float((float)1,(float)1));
  mncblas_cgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, &Cfalpha, CfA, 0, CfB, 1, &Cfbeta, CfC, 1);
  start = _rdtsc () ;
  for (int k=0; k<NB_FOIS; k++){
    mncblas_cgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, &Cfalpha, CfA, 0, CfB, 1, &Cfbeta, CfC, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_cgemm", (N * N * (12 + N * 8) + 11) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_zgemm\n##################\n");
  dmatriceC_init(CdA, init_complexe_double((double)1,(double)1));
  dmatriceC_init(CdB, init_complexe_double((double)1,(double)1));
  dmatriceC_init(CdC, init_complexe_double((double)1,(double)1));
  mncblas_zgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, &Cdalpha, CdA, 0, CdB, 1, &Cdbeta, CdC, 1);
  start = _rdtsc () ;
  for (int k=0; k<NB_FOIS; k++){
    mncblas_zgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, &Cdalpha, CdA, 0, dB, 1, &Cdbeta, CdC, 1);
  }
  end = _rdtsc ();
  calcul_flop ("Tests de mncblas_zgemm",  (N * N * (12 + N * 8) + 11) * (unsigned long) NB_FOIS , end-start) ;
}
