#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"

#include "flop.h"

#define N 258
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

void fvectorC_init(float *V, float x) {
  for (int i = 0; i < N ; i++) {
    V[2 * i] = x;
    V[2 * i + 1] = x;
  }
}

void dvectorC_init(double *V, double x) {
  for (int i = 0; i < N; i++) {
    V[2 * i] = x;
    V[2 * i + 1] = x;
  }
}


int main(int argc, char **argv) {
  int i;
  float *fA;
  double *dA;
  float *fX;
  double *dX;
  float *fY;
  double *dY;

  float falpha[2] = {(float)1 , (float)1};
  double dalpha[2] = {(double)1, (double)1};
  float fbeta[2] = {(float)1, (float)1};
  double dbeta[2] = {(double)1, (double)1};

  unsigned long long int start, end ;

  init_flop () ;
  printf("\nTests des fonctions de gemv\n\n");

  printf("\n##################\nTests de mncblas_sgemv\n##################\n");
  fA = (float *) malloc (sizeof(float) * N * N);
  fX = (float *) malloc (sizeof(float) * N);
  fY = (float *) malloc (sizeof(float) * N);
  fmatrice_init(fA, (float)1);
  fvector_init(fX, (float)1);
  fvector_init(fY, (float)1);
  start = _rdtsc () ;
  for (i = 0; i < NB_FOIS; i ++) {
      mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, falpha[0], fA, 0, fX, 1, fbeta[0], fY, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_sgemv", N * (2 * N + 3) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_dgemv\n##################\n");
  dA = (double *) malloc (sizeof(double) * N * N);
  dX = (double *) malloc (sizeof(double) * N);
  dY = (double *) malloc (sizeof(double) * N);
  dmatrice_init(dA, (double)1);
  dvector_init(dX, (double)1);
  dvector_init(dY, (double)1);
  start = _rdtsc () ;
  for (i = 0; i < NB_FOIS; i ++) {
    mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, dalpha[0], dA, 0, dX, 1, dbeta[0], dY, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_dgemv", N * (2 * N +3) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_cgemv\n##################\n");
  free(fA);
  free(fX);
  free(fY);
  fA = (float *) malloc (sizeof(float) * N * N * 2);
  fX = (float *) malloc (sizeof(float) * N * 2);
  fY = (float *) malloc (sizeof(float) * N * 2);
  fmatriceC_init(fA, (float)1);
  fvectorC_init(fX, (float)1);
  fvectorC_init(fY, (float)1);
  start = _rdtsc () ;
  for (i = 0; i < NB_FOIS; i ++) {
    mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, falpha, fA, 0, fX, 1, fbeta, fY, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_cgemv", N * (14 + (N * 8)) * NB_FOIS, end-start) ;

  printf("\n##################\nTests de mncblas_zgemv\n##################\n");
  free(dA);
  free(dX);
  free(dY);
  dA = (double *) malloc (sizeof(double) * N * N * 2);
  dX = (double *) malloc (sizeof(double) * N * 2);
  dY = (double *) malloc (sizeof(double) * N * 2);
  dmatriceC_init(dA, (double)1);
  dvectorC_init(dX, (double)1);
  dvectorC_init(dY, (double)1);
  start = _rdtsc () ;
  for (i = 0; i < NB_FOIS; i ++) {
    mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, N, N, dalpha, dA, 0, dX, 1, dbeta, dY, 1);
  }
  end = _rdtsc () ;
  calcul_flop ("Tests de mncblas_zgemv", N * (14 + (N * 8)) * NB_FOIS, end-start) ;

  free(fA);
  free(fX);
  free(fY);
  free(dA);
  free(dX);
  free(dY);
}
