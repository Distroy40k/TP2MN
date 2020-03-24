#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"

#include "flop.h"

#define N 258
#define NB_FOIS 10


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
void fmatriceC_init(float *A, float x) {
  for (int i = 0; i < N * N; i++) {
    A[2 * i] = x;
    A[2 * i + 1] = x;
  }
}
void dmatriceC_init(double *A, double x) {
  for (int i = 0; i < N * N; i++) {
    A[2 * i] = x;
    A[2 * i + 1] = x;
  }
}


int main(int argc, char **argv) {

  float *fA;
  double *dA;
  float *fB;
  double *dB;
  float *fC;
  double *dC;

  float falpha[2] = {(float)1 , (float)1};
  double dalpha[2] = {(double)1, (double)1};
  float fbeta[2] = {(float)1, (float)1};
  double dbeta[2] = {(double)1, (double)1};

  unsigned long long int start, end ;

  init_flop () ;
  printf("\nTests des fonctions de gemm\n\n");

  printf("\n##################\nTests de mncblas_sgemm\n##################\n");
  fA = (float *) malloc (sizeof(float) * N * N);
  fB = (float *) malloc (sizeof(float) * N * N);
  fC = (float *) malloc (sizeof(float) * N * N);
  fmatrice_init(fA, (float)1);
  fmatrice_init(fB, (float)1);
  fmatrice_init(fC, (float)1);
  start = _rdtsc () ;
  for (int k=0; k<NB_FOIS; k++){
    mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, falpha[0], fA, 0, fB, 1, fbeta[0], fC, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_sgemm", N * N * ((N * 2) + 3 ) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_dgemm\n##################\n");
  dA = (double *) malloc (sizeof(double) * N * N);
  dB = (double *) malloc (sizeof(double) * N * N);
  dC = (double *) malloc (sizeof(double) * N * N);
  dmatrice_init(dA, (double)1);
  dmatrice_init(dB, (double)1);
  dmatrice_init(dC, (double)1);
  start = _rdtsc () ;
  for (int k=0; k<NB_FOIS; k++){
    mncblas_dgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, dalpha[0], dA, 0, dB, 1, dbeta[0], dC, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_dgemm", N * N * ((N * 2) + 3 ) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_cgemm\n##################\n");
  free(fA);
  free(fB);
  free(fC);
  fA = (float *) malloc (sizeof(float) * N * N * 2);
  fB = (float *) malloc (sizeof(float) * N * N * 2);
  fC = (float *) malloc (sizeof(float) * N * N * 2);
  fmatriceC_init(fA, (float)1);
  fmatriceC_init(fB, (float)1);
  fmatriceC_init(fC, (float)1);
  start = _rdtsc () ;
  for (int k=0; k<NB_FOIS; k++){
    mncblas_cgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, falpha, fA, 0, fB, 1, fbeta, fC, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_cgemm", N * N * ((N * 8) + 14) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_zgemm\n##################\n");
  free(dA);
  free(dB);
  free(dC);
  dA = (double *) malloc (sizeof(double) * N * N * 2);
  dB = (double *) malloc (sizeof(double) * N * N * 2);
  dC = (double *) malloc (sizeof(double) * N * N * 2);
  dmatriceC_init(dA, (double)1);
  dmatriceC_init(dB, (double)1);
  dmatriceC_init(dC, (double)1);
  start = _rdtsc () ;
  for (int k=0; k<NB_FOIS; k++){
    mncblas_zgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, N, N, N, dalpha, dA, 0, dB, 1, dbeta, dC, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_zgemm",  N * N * ((N * 8) + 14) * NB_FOIS , end-start) ;

  free(fA);
  free(fB);
  free(fC);
  free(dA);
  free(dB);
  free(dC);
}
