#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"

#include "flop.h"

#define NB_FOIS 419430
#define SIZE_VECTOR 100

int is_equal_float(float *f1, float *f2, const int N) {
  for (int i = 0; i < N; i++) {
    if (f1[i] != f2[i]) {
      return 0;
    }
  }
  return 1;
}

int is_equal_double(double *d1, double *d2, const int N) {
  for (int i = 0; i < N; i++) {
    if (d1[i] != d2[i]) {
      return 0;
    }
  }
  return 1;
}

int main (int argc, char **argv)
{

  if (argc != 2) {
    printf("Utilisation: ./test_iamx mode\n");
    exit(1);
  }
  char mode = argv[1][0];

 float f1[SIZE_VECTOR];
 float f2[SIZE_VECTOR];
 double d1[SIZE_VECTOR];
 double d2[SIZE_VECTOR];
 unsigned long long int start, end ;
 int i;
 printf("\nTests des fonctions de copy\n\n");
 printf("####################\nPour tester ces fonctions, on cré un vecteur de taille %d, avec des élements de 0 à %d, et on le copie dans un autre vecteur\nOn verifie ensuite que la copie c'est ben passées####################\n\n", SIZE_VECTOR, SIZE_VECTOR);

switch (mode) {
  case 's':
    printf("\nTests de mncblas_scopy\n\n");
    for (int var = 0; var < SIZE_VECTOR; var += 1) {
      f1[var] = (float) var;
    }
    init_flop () ;
    start = _rdtsc () ;
    for (i = 0; i < NB_FOIS; i ++) {
       mncblas_scopy(SIZE_VECTOR, f1, 1, f2, 1);
    }
    end = _rdtsc () ;
    printf ("La copie d'un tableau de vecteur de float simple précision, de 0 à %d: %lld cycles \n", SIZE_VECTOR, end-start) ;
    calcul_octet ("mncblas_scopy ", NB_FOIS*4*SIZE_VECTOR, end-start) ;
    if (is_equal_float(f1,f2, SIZE_VECTOR) != 1) {
      printf("Erreur copie mncblas_scopy\n");
      exit(1);
    }
    break;
  case 'd':
    printf("\n");
    printf("Tests de mncblas_dcopy\n");
    printf("\n");

    for (int var = 0; var < SIZE_VECTOR; var += 1) {
      d1[var] = (double) var;
    }
    start = _rdtsc () ;
    for (i = 0; i < NB_FOIS; i ++) {
       mncblas_dcopy(SIZE_VECTOR, d1, 1, d2, 1);
    }
    end = _rdtsc () ;
    printf ("La copie d'un tableau de vecteur de double simple précision, de 0 à %d: %lld cycles \n", SIZE_VECTOR, end-start) ;
    calcul_octet ("mncblas_dcopy ", NB_FOIS*4, end-start) ;
    if (is_equal_double(d1,d2, SIZE_VECTOR) != 1) {
      printf("Erreur copie mncblas_dcopy\n");
      exit(1);
    }
    break;
  case 'c':
    printf("\n");
    printf("Tests de mncblas_ccopy\n");
    printf("\n");
    init_flop () ;
    start = _rdtsc () ;
    for (i = 0; i < NB_FOIS; i ++) {
      mncblas_ccopy(SIZE_VECTOR / 2, f1, 1, f2, 1);
    }
    end = _rdtsc () ;
    printf ("La copie d'un tableau de vecteur de float double précision, de 0 à %d:  %lld cycles \n", SIZE_VECTOR, end-start) ;
    calcul_octet ("mncblas_ccopy ", NB_FOIS*4, end-start) ;
    if (is_equal_float(f1,f2, SIZE_VECTOR) != 1) {
      printf("Erreur copie mncblas_ccopy\n");
      exit(1);
    }
    break;
  case 'z':
    printf("\n");
    printf("Tests de mncblas_zcopy\n");
    printf("\n");
    start = _rdtsc () ;
    for (i = 0; i < NB_FOIS; i ++) {
    mncblas_zcopy(SIZE_VECTOR / 2, d1, 1, d2, 1);
    }
    end = _rdtsc () ;
    printf ("La copie d'un tableau de vecteur de double double précision, de 0 à %d: %lld cycles \n", SIZE_VECTOR, end-start) ;
    calcul_octet ("mncblas_zcopy ", NB_FOIS*4, end-start) ;
    if (is_equal_double(d1,d2, SIZE_VECTOR) != 1) {
      printf("Erreur copie mncblas_zcopy\n");
      exit(1);
    }
    break;
  default:
    printf("Mode non reconnu\n");
    exit(1);
}
printf("All test passed\n");
  exit (0) ;
}
