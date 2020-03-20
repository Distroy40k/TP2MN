#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"

#include "flop.h"

#define NB_FOIS 45126

void fvector_init(FILE *f, float *V, int n) {
  for (int i = 0; i < n; i++) {
    if (fscanf(f, "%f", V + i) != 1){
        printf("Failed to read fvector.\n");
    }
  }
}

void dvector_init(FILE *f, double *V, int n) {
  for (int i = 0; i < n; i++) {
    if (fscanf(f, "%lf", V + i) != 1){
        printf("Failed to read dvector.\n");
    }
  }
}

void fmatrice_init(FILE* f, float *A, int n) {
  for (int i = 0; i < n * n; i++) {
      if (fscanf(f, "%f", &(A[i])) != 1){
        printf("Failed to read fmatrice.\n");
    }
  }
}

void dmatrice_init(FILE* f, double *A, int n) {
  for (int i = 0; i < n * n; i++) {
      if (fscanf(f, "%lf", A + i) != 1){
        printf("Failed to read dmatrice.\n");
    }
  }
}

void fvectorC_init(FILE *f, float *V, int n) {
  for (int i = 0; i < n; i++) {
      if (fscanf(f, "%f", V + 2 * i) != 1){
        printf("Failed to read fvectorC.\n");
    }
    if (fscanf(f, "%f", V + 2 * i + 1) != 1){
        printf("Failed to read fvectorC.\n");
    }
  }
}

void dvectorC_init(FILE *f, double *V, int n) {
  for (int i = 0; i < n ; i++) {
    if (fscanf(f, "%lf", V + 2 * i) != 1){
        printf("Failed to read dvectorC.\n");
    }
    if (fscanf(f, "%lf", V + 2 * i + 1) != 1){
        printf("Failed to read dvectorC.\n");
    }
  }
}

void fmatriceC_init(FILE* f, float *A, int n) {

  for (int i = 0; i < n * n ; i++) {
    if (fscanf(f, "%f", &A[2 * i]) != 1){
        printf("Failed to read fmatriceC.\n");
    }
    if (fscanf(f, "%f", &A[2 * i + 1]) != 1){
        printf("Failed to read fmatriceC.\n");
    }
  }
}

void dmatriceC_init(FILE* f, double *A, int n) {
  for (int i = 0; i < n * n; i++) {
    if (fscanf(f, "%lf", A + 2 * i) != 1){
        printf("Failed to read dmatriceC.\n");
    }
    if (fscanf(f, "%lf", A + 2 * i + 1) != 1){
        printf("Failed to read dmatriceC.\n");
    }
  }
}

void afficher_scalc(float alpha, float *A, float * B, float beta, float *C, int n) {
  printf("n : %d\n", n);
  int i = 0;
  int j=0;
  printf("Matrice A\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", A[i + j*n]);
    }
    printf("\n");
  }
  printf("Matrice B\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", B[i + j*n]);
    }
    printf("\n");
  }
  printf("Matrice C\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", C[i + j*n]);
    }
    printf("\n");
  }
  printf("Alpha %f\n", alpha);
  printf("Beta %f\n", beta);
}

void afficher_ccalc(float *alpha, float *A, float * B, float *beta, float *C, int n) {
  int i = 0;
  int j = 0;
  printf("Matrice A\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", A[2 * i + j*n]);
      printf("%f ", A[2 * i + 1 + j*n]);
    }
    printf("\n");
  }
  printf("Matrice B\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", B[2 * i + j*n]);
      printf("%f ", B[2 * i + 1 + j*n]);
    }
    printf("\n");
  }
  printf("Matrice C\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", C[2 * i + j*n]);
      printf("%f ", C[2 * i + 1 + j*n]);
    }
    printf("\n");
  }
  printf("Alpha %f, %f\n", alpha[0], alpha[1]);
  printf("Beta %f, %f\n", beta[0], beta[1]);
}


void afficher_dcalc(double alpha, double *A, double * B, double beta, double *C, int n) {
  printf("Matrice A\n");
  int i = 0;
  int j=0;
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", A[i + j*n]);
    }
    printf("\n");
  }
  printf("Matrice B\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", B[i + j*n]);
    }
    printf("\n");
  }
  printf("Matrice C\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", C[i + j*n]);
    }
    printf("\n");
  }
  printf("Alpha %f\n", alpha);
  printf("Beta %f\n", beta);
}

void afficher_zcalc(double *alpha, double *A, double * B, double *beta, double *C, int n) {
  int i = 0;
  int j = 0;
  printf("Matrice A\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", A[2 * i + j*n]);
      printf("%f ", A[2 * i + 1 + j*n]);
    }
    printf("\n");
  }
  printf("Matrice B\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", B[2 * i + j*n]);
      printf("%f ", B[2 * i + 1 + j*n]);
    }
    printf("\n");
  }
  printf("Matrice C\n");
  for (j = 0; j < n ; j++){
    for (i = 0; i < n ; i ++) {
      printf("%f ", C[2 * i + j*n]);
      printf("%f ", C[2 * i + 1 + j*n]);
    }
    printf("\n");
  }
  printf("Alpha %f, %f\n", alpha[0], alpha[1]);
  printf("Beta %f, %f\n", beta[0], beta[1]);
}

int main(int argc, char **argv) {

  if (argc != 8) {
    printf("Utilisation:\ntest_gemm.c mode alpha A B beta C N\n");
    printf("Les arguments autres que le mode et N doivent être des fichiers\n");
    printf("Mode :\n s single real\n d double real\n c single complex\n z double complex\n");
    printf("Renseigner la taille du vecteur et de la matrice au début du fichier\n");
    exit(1);
  }
  int n;
  sscanf(argv[7], "%d", &n);
  char mode = argv[1][0];
  float *fA;
  double *dA;
  float *fB;
  double *dB;
  float *fC;
  double *dC;

  switch (mode) {
    case 's':
    case 'd':
    fA = (float *) malloc(sizeof(float) * n * n);
    dA = (double *) malloc(sizeof(double) * n * n);
    fB = (float *) malloc(sizeof(float) * n * n);
    dB = (double *) malloc(sizeof(double) * n * n);
    fC = (float *) malloc(sizeof(float) * n * n);
    dC = (double *) malloc(sizeof(double) * n * n);
    break;
    case 'c':
    case 'z':
    default:
    fA = (float *) malloc(sizeof(float) * n * n * 2);
    dA = (double *) malloc(sizeof(double) * n * n * 2);
    fB = (float *) malloc(sizeof(float) * n * n * 2);
    dB = (double *) malloc(sizeof(double) * n * n * 2);
    fC = (float *) malloc(sizeof(float) * n * n * 2);
    dC = (double *) malloc(sizeof(double) * n * n * 2);
    break;
  }

  float falpha[2];
  double dalpha[2];
  float fbeta[2];
  double dbeta[2];

  if (mode == 's') {
    fmatrice_init(fopen(argv[3], "r"), fA, n);
    fmatrice_init(fopen(argv[4], "r"), fB, n);
    fmatrice_init(fopen(argv[6], "r"), fC, n);
  } else if (mode == 'd') {
    dmatrice_init(fopen(argv[3], "r"), dA, n);
    dmatrice_init(fopen(argv[4], "r"), dB, n);
    dmatrice_init(fopen(argv[6], "r"), dC, n);
  } else if (mode == 'c') {
    fmatriceC_init(fopen(argv[3], "r"), fA, n);
    fmatriceC_init(fopen(argv[4], "r"), fB, n);
    fmatriceC_init(fopen(argv[6], "r"), fC, n);
  } else if (mode == 'z') {
    dmatriceC_init(fopen(argv[3], "r"), dA, n);
    dmatriceC_init(fopen(argv[4], "r"), dB, n);
    dmatriceC_init(fopen(argv[6], "r"), dC, n);
  }

  if (mode == 's') {
    if (fscanf(fopen(argv[2], "r"), "%f", falpha) != 1){
        printf("Failed to read falpha.\n");
    }
    if (fscanf(fopen(argv[5], "r"), "%f", fbeta) != 1){
        printf("Failed to read fbeta.\n");
    }
  } else  if (mode == 'd') {
    if (fscanf(fopen(argv[2], "r"), "%lf", dalpha) != 1){
        printf("Failed to read dalpha.\n");
    }
    if (fscanf(fopen(argv[5], "r"), "%lf", dbeta) != 1){
        printf("Failed to read dbeta.\n");
    }
  } else if (mode == 'c') {
    FILE *Falpha = fopen(argv[2], "r");
    FILE *Fbeta = fopen(argv[5], "r");
    
    if ((fscanf(Falpha, "%f", falpha) != 1) || (fscanf(Falpha, "%f", falpha + 1) != 1)){
        printf("Failed to read falpha.\n");
    }
    if ((fscanf(Fbeta, "%f", fbeta) != 1) || (fscanf(Fbeta, "%f", fbeta + 1) != 1)){
        printf("Failed to read fbeta.\n");
    }
  } else if (mode == 'z') {
    FILE *Falpha = fopen(argv[2], "r");
    FILE *Fbeta = fopen(argv[5], "r");
    if ((fscanf(Falpha, "%lf", dalpha) != 1) || (fscanf(Falpha, "%lf", dalpha + 1) != 1)){
        printf("Failed to read dalpha.\n");
    }
    if ((fscanf(Fbeta, "%lf", dbeta) != 1) || (fscanf(Fbeta, "%lf", dbeta + 1) != 1)){
        printf("Failed to read dbeta.\n");
    }
  }

  unsigned long long int start, end ;
  init_flop () ;
  printf("\nTests des fonctions de gemm\n\n");

  switch (mode) {
    
    case 's':
      printf("\n##################\nTests de mncblas_sgemm\n##################\n");
      afficher_scalc(falpha[0], fA, fB, fbeta[0], fC, n);
      mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, n, n, n, falpha[0], fA, 0, fB, 1, fbeta[0], fC, 1);
      printf("\nResultat:\n");
      for (int j = 0; j < n ; j++){
        for (int i =0; i < n ; i++) {
          printf("%f ", fC[i]);
        }
      printf("\n");
      }
      start = _rdtsc () ;
      for (int k=0; k<NB_FOIS; k++){
        mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, n, n, n, falpha[0], fA, 0, fB, 1, fbeta[0], fC, 1);
      }
      end = _rdtsc () ;
      calcul_flop ("Tests de mncblas_sgemm", n * n * ((n * 2) + 3 ) * NB_FOIS, end-start) ;
      break;

    case 'd':
      printf("\n##################\nTests de mncblas_dgemm\n##################\n");
      afficher_dcalc(dalpha[0], dA, dB, dbeta[0], dC, n);
      mncblas_dgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, n, n, n, dalpha[0], dA, 0, dB, 1, dbeta[0], dC, 1);
      printf("\nResultat:\n");
      for (int j = 0; j < n ; j++){
        for (int i =0; i < n ; i++) {
          printf("%f ", dC[i]);
        }
      printf("\n");
      }
      start = _rdtsc () ;
      for (int k=0; k<NB_FOIS; k++){
        mncblas_dgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, n, n, n, dalpha[0], dA, 0, dB, 1, dbeta[0], dC, 1);
      }
      end = _rdtsc () ;
      calcul_flop ("Tests de mncblas_dgemm", n * n * ((n * 2) + 3 ) * NB_FOIS, end-start) ;
      break;

    case 'c':
      printf("\n##################\nTests de mncblas_cgemm\n##################\n");
      afficher_ccalc(falpha, fA, fB, fbeta, fC,  n);
      mncblas_cgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, n, n, n, falpha, fA, 0, fB, 1, fbeta, fC, 1);
      printf("\nResultat:\n");
      for (int j = 0; j< n ; j++){
        for (int i =0; i <  n ; i++) {
          printf("%f ", fC[2 *i]);
          printf("%f ", fC[2 * i + 1]);
        }
      printf("\n");
      }
      start = _rdtsc () ;
      for (int k=0; k<NB_FOIS; k++){
        mncblas_cgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, n, n, n, falpha, fA, 0, fB, 1, fbeta, fC, 1);
      }
      end = _rdtsc () ;
      calcul_flop ("Tests de mncblas_cgemm", n * n * ((n * 8) + 14) * NB_FOIS, end-start) ;
      break;

    case 'z':
      printf("\n##################\nTests de mncblas_zgemv\n##################\n");
      afficher_zcalc(dalpha, dA, dB, dbeta, dC, n);
      mncblas_zgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, n, n, n, dalpha, dA, 0, dB, 1, dbeta, dC, 1);
      printf("\nResultat:\n");
      for (int j = 0; j< n ; j++){
        for (int i =0; i < n ; i++) {
          printf("%f ", dC[2 *i]);
          printf("%f ", dC[2 * i + 1]);
        }
      printf("\n");
      }
      start = _rdtsc () ;
      for (int k=0; k<NB_FOIS; k++){
        mncblas_zgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, n, n, n, dalpha, dA, 0, dB, 1, dbeta, dC, 1);
      }
      end = _rdtsc () ;
      calcul_flop ("Tests de mncblas_zgemm",  n * n * ((n * 8) + 14) * NB_FOIS , end-start) ;
      break;

    default:
      printf("Mode non reconnu\n");
      exit(1);
  }
}