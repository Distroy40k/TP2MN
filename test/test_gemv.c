#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"

#include "flop.h"

#define NB_FOIS 45126

void fvector_init(FILE *f, float *V, int n) {
  for (int i = 0; i < n; i++) {
    fscanf(f, "%f", V + i);
  }
}

void dvector_init(FILE *f, double *V, int n) {
  for (int i = 0; i < n; i++) {
    fscanf(f, "%lf", V + i);
  }
}

void fmatrice_init(FILE* f, float *A, int n) {
  for (int i = 0; i < n * n; i++) {
    fscanf(f, "%f", &(A[i]));
  }
}

void dmatrice_init(FILE* f, double *A, int n) {
  for (int i = 0; i < n * n; i++) {
    fscanf(f, "%lf", A + i);
  }
}

void afficher_scalc(float alpha, float *A, float * X, float beta, float *Y, int n) {
  printf("n : %d\n", n);
  int i = 0;
  printf("Matrice A\n");

  for (i = 0; i < n * 2; i ++) {
      printf("%f\n", A[i]);
    }
  printf("\nX:\n");
  for (i =0; i < n; i++) {
    printf("%f\n", X[i]);
  }
  printf("\nY:\n");
  for (i =0; i < n; i++) {
    printf("%f\n", Y[i]);
  }
  printf("Alpha %f\n", alpha);
  printf("Beta %f\n", beta);
}

void afficher_ccalc(float *alpha, float *A, float * X, float *beta, float *Y, int n) {
  printf("Matrice A\n");
  int i = 0;
  for (i = 0; i < n * n; i ++) {
      printf("%f\n", A[i]);
  }
  printf("\nX:\n");
  for (i =0; i < n; i++) {
    printf("%f\n", X[i]);
  }
  printf("\nY:\n");
  for (i =0; i < n; i++) {
    printf("%f\n", Y[i]);
  }
  printf("Alpha %f %f\n", alpha[0], alpha[1]);
  printf("Beta %f %f\n", beta[0], beta[1]);

}


void afficher_dcalc(double alpha, double *A, double * X, double beta, double *Y, int n) {
  printf("Matrice A\n");
  int i = 0;
  for (i = 0; i < n * n; i ++) {
      printf("%f\n", A[i]);
  }
  printf("\nX:\n");
  for (i =0; i < n; i++) {
    printf("%f\n", X[i]);
  }
  printf("\nY:\n");
  for (i =0; i < n; i++) {
    printf("%f\n", Y[i]);
  }
  printf("Alpha %f\n", alpha);
  printf("Beta %f\n", beta);
}

void afficher_zcalc(double *alpha, double *A, double * X, double *beta, double *Y, int n) {
  printf("Matrice A\n");
  int i = 0;
  for (i = 0; i < n * n; i ++) {
    printf("%f\n", A[i]);
  }
  printf("\nX:\n");
  for (i =0; i < n; i++) {
    printf("%f\n", X[i]);
  }
  printf("\nY:\n");
  for (i =0; i < n; i++) {
    printf("%f\n", Y[i]);
  }
  printf("Alpha %f, %f\n", alpha[0], alpha[1]);
  printf("Beta %f, %f\n", beta[0], beta[1]);
}

int main(int argc, char **argv) {

  if (argc != 8) {
    printf("Utilisation:\ntest_gemv.c mode alpha A X beta Y N\n");
    printf("Les arguments autres que le mode et N doivent être des fichiers\n");
    printf("Mode :\n s single real\n d double real\n c single complex\n z double complex\n");
    printf("Renseigner la taille du vecteur et de la matrice au début du fichier\n");
    exit(1);
  }
  int n;
  sscanf(argv[7], "%d", &n);
  printf("%d\n", n);
  char mode = argv[1][0];
  float *fA = (float *) malloc(sizeof(float) * n * n);
  double dA[n*n];
  float fX[n*n];
  double dX[n*n];
  float fY[n*n];
  double dY[n*n];
  float falpha[2];
  double dalpha[2];
  float fbeta[2];
  double dbeta[2];

  if (mode == 's' || mode == 'c') {
    fmatrice_init(fopen(argv[3], "r"), fA, n);
    printf("A\n");
    fvector_init(fopen(argv[4], "r"), fX, n);
    fvector_init(fopen(argv[6], "r"), fY, n);
  } else {
    dmatrice_init(fopen(argv[3], "r"), dA, n);
    dvector_init(fopen(argv[4], "r"), dX, n);
    dvector_init(fopen(argv[6], "r"), dY, n);
  }

  if (mode == 's') {
    fscanf(fopen(argv[2], "r"), "%f", falpha);
    fscanf(fopen(argv[5], "r"), "%f", fbeta);
  } else  if (mode == 'd') {
    fscanf(fopen(argv[2], "r"), "%lf", dalpha);
    fscanf(fopen(argv[5], "r"), "%lf", dbeta);
  } else if (mode == 'c') {
    FILE *Falpha = fopen(argv[2], "r");
    FILE *Fbeta = fopen(argv[5], "r");
    fscanf(Falpha, "%f", falpha);
    fscanf(Falpha, "%f", falpha + 1);
    fscanf(Fbeta, "%f", fbeta);
    fscanf(Fbeta, "%f", fbeta + 1);
  } else if (mode == 'z') {
    FILE *Falpha = fopen(argv[2], "r");
    FILE *Fbeta = fopen(argv[5], "r");
    fscanf(Falpha, "%lf", dalpha);
    fscanf(Falpha, "%lf", dalpha + 1);
    fscanf(Fbeta, "%lf", dbeta);
    fscanf(Fbeta, "%lf", dbeta + 1);
  }



  unsigned long long int start, end ;

  init_flop () ;

  printf("\nTests des fonctions de gemv\n\n");

  switch (mode) {
    case 's':
      printf("\n##################\nTests de mncblas_sgemv\n##################\n");
      afficher_scalc(falpha[0], fA, fX, fbeta[0], fY, n);

      start = _rdtsc () ;
      mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, n, n, falpha[0], fA, 0, fX, 1, fbeta[0], fY, 1);
      end = _rdtsc () ;

      printf("\nResultat:\n");
      for (int i =0; i < n; i++) {
        printf("%f\n", fY[i]);
      }

      calcul_flop ("Tests de mncblas_sgemv", n * n * NB_FOIS, end-start) ;
      break;
    case 'd':
      printf("\n##################\nTests de mncblas_dgemv\n##################\n");
      afficher_dcalc(falpha[0], dA, dX, fbeta[0], dY, n);
      start = _rdtsc () ;
      mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, n, n, dalpha[0], dA, 0, dX, 1, dbeta[0], dY, 1);
      end = _rdtsc () ;
      printf("\nResultat:\n");
      for (int i =0; i < n; i++) {
        printf("%f\n", fY[i]);
      }

      calcul_flop ("Tests de mncblas_dgemv", n * n * NB_FOIS, end-start) ;
      break;
    case 'c':
      printf("\n##################\nTests de mncblas_cgemv\n##################\n");
      afficher_ccalc(falpha, fA, fX, fbeta, fY, n);
      start = _rdtsc () ;
      mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, n/2, n/2, falpha, fA, 0, fX, 1, fbeta, fY, 1);
      end = _rdtsc () ;
      printf("\nResultat:\n");
      for (int i =0; i < n; i++) {
        printf("%f\n", fY[i]);
      }

      calcul_flop ("Tests de mncblas_sgemv", n * n * NB_FOIS, end-start) ;
      break;
    case 'z':
      printf("\n##################\nTests de mncblas_zgemv\n##################\n");
      afficher_zcalc(dalpha, dA, dX, dbeta, dY, n);
      start = _rdtsc () ;
      mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, n/2, n/2, dalpha, dA, 0, dX, 1, dbeta, dY, 1);
      end = _rdtsc () ;
      printf("\nResultat:\n");
      for (int i =0; i < n; i++) {
        printf("%f\n", fY[i]);
      }
      calcul_flop ("Tests de mncblas_dgemv", n * n * NB_FOIS, end-start) ;
      break;
  }
}
