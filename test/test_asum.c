#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"

#include "flop.h"

#define NB_TEST 5
#define NB_FOIS 41943
#define SIZE_VECTOR 426

void init_vecteurf(float* v, float x) {
  for (int i = 0; i < SIZE_VECTOR; i++) {
    v[i] = x;
  }
}

void init_vecteurd(double* v, double x) {
  for (int i = 0; i < SIZE_VECTOR; i++) {
    v[i] = x;
  }
}

int main (int argc, char **argv)
{

  if (argc != 2) {
    printf("Utilisation: ./test_iamx mode\n");
    exit(1);
  }
  char mode = argv[1][0];

 float fvector[SIZE_VECTOR];
 double dvector[SIZE_VECTOR];

 unsigned long long int start, end ;
 int i;
 int test;
 float resf;
 double resd;
 init_flop () ;

 printf("\n");
 printf("Tests des fonctions de asum\n");
 printf("\n#######################\nPour tester ces fonctions, on initialise un vecteur de taille %d avec pour chaque coordonées le numéros du test\net on vérifie que la somme est correctement calculée\n###################", SIZE_VECTOR);

 switch (mode) {
   case 's':
     printf("\n##################\nTests de mnblas_sasum\n##################\n");
      for (test = 0; test < NB_TEST; test ++) {
        printf("\n########### TEST %d############\n", test);
        init_vecteurf(fvector, (float) test);
        start = _rdtsc () ;
        for (i = 0; i < NB_FOIS; i ++) {
           resf = mnblas_sasum(SIZE_VECTOR, fvector, 1);
        }
        end = _rdtsc () ;
        printf ("La somme d'un vecteur de float simple précision: %lld cycles \n", end-start) ;
        calcul_flop ("mnblas_sasum ", NB_FOIS, end-start) ;
        if (resf != test * SIZE_VECTOR) {
          printf("Erreur copie mnblas_sasum\n");
          exit(1);
        }
      }
      break;
    case 'd':
      printf("\n##################\nTests de mnblas_dasum\n##################\n");
      for (test = 0; test < NB_TEST; test ++) {
        printf("\n########### TEST %d ############\n", test);
        init_vecteurd(dvector, (float) test);
        start = _rdtsc () ;
        for (i = 0; i < NB_FOIS; i ++) {
           resd = mnblas_dasum(SIZE_VECTOR, dvector, 1);
        }
        end = _rdtsc () ;
        printf ("La somme d'un vecteur de float double précision: %lld cycles \n", end-start) ;
        calcul_flop ("mnblas_dasum ", NB_FOIS, end-start) ;
        if (resd != test * SIZE_VECTOR) {
          printf("Erreur copie mnblas_dasum\n");
          exit(1);
        }
      }
      break;
    case 'c':
      printf("\n##################\nTests de mnblas_scasum\n##################\n");
      for (test = 0; test < NB_TEST; test ++) {
        printf("\n########### TEST %d############\n",test);
        init_vecteurf(fvector, (float)test);
        start = _rdtsc () ;
        for (i = 0; i < NB_FOIS; i ++) {
          resf = mnblas_scasum(SIZE_VECTOR/2, fvector, 1);
        }
        end = _rdtsc () ;
        printf ("La somme d'un vecteur de float simple précision: %lld cycles \n", end-start) ;
        calcul_flop ("mnblas_scasum ", NB_FOIS, end-start) ;
        if (resf != test * SIZE_VECTOR) {
         printf("Erreur copie mnblas_scasum\n");
         exit(1);
        }
      }
    case 'z':
      printf("\n##################\nTests de mnblas_dzasum\n##################\n");
        for (test = 0; test < NB_TEST; test ++) {
        printf("\n########### TEST %d ############\n", test);
        init_vecteurd(dvector, (float) test);
        start = _rdtsc () ;
        for (i = 0; i < NB_FOIS; i ++) {
          resd = mnblas_dzasum(SIZE_VECTOR/2, dvector, 1);
        }
        end = _rdtsc () ;
        printf ("La somme d'un vecteur de double double précision: %lld cycles \n", end-start) ;
        calcul_flop ("mnblas_dzasum ", NB_FOIS, end-start) ;
        if (resd != test * SIZE_VECTOR) {
         printf("Erreur copie mnblas_dzasum\n");
         exit(1);
        }
      }
      break;
    default:
      printf("Mode incorrect\n");
      exit(1);
  }
  printf("All test passed\n");
  exit (0) ;
}
