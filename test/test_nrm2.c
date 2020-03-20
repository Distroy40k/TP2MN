#include <stdio.h>
#include <x86intrin.h>
#include <math.h>

#include "mnblas.h"
#include "flop.h"

#define NB_TEST 10
#define NB_FOIS 4194
#define SIZE_VECTOR 4264

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
 printf("\n");
 printf("\n#######################\nPour tester ces fonctions, on initialise un vecteur de taille %d avec pour chaque coordonées le numéros du test\net on vérifie que la norme est correctement calculée\n###################\n", SIZE_VECTOR);


 switch (mode) {
   case 's':
     printf("\n##################\nTests de mnblas_snrm2\n##################\n");
      for (test = 0; test < NB_TEST; test ++) {
        printf("\n########### TEST %d############\n", test);
        init_vecteurf(fvector, (float) test);
        start = _rdtsc () ;
        for (i = 0; i < NB_FOIS; i ++) {
           resf = mnblas_snrm2(SIZE_VECTOR, fvector, 1);
        }
        end = _rdtsc () ;
        printf ("La norme euclidienne d'un vecteur: %lld cycles \n", end-start) ;
        calcul_flop ("mnblas_snrm2 ", 2* NB_FOIS, end-start) ;
        if (resf != (float)sqrt(pow(test, 2) * (float) SIZE_VECTOR)) {;
          printf("Erreur mnblas_snrm2\n");
          exit(1);
        }
      }
      break;
    case 'd':
      printf("\n##################\nTests de mnblas_dnrm2\n##################\n");
       for (test = 0; test < NB_TEST; test ++) {
         printf("\n########### TEST %d############\n", test);
         init_vecteurd(dvector, (double) test);
         start = _rdtsc () ;
         for (i = 0; i < NB_FOIS; i ++) {
            resd = mnblas_dnrm2(SIZE_VECTOR, dvector, 1);
         }
         end = _rdtsc () ;
         printf ("La norme euclidienne d'un vecteur: %lld cycles \n", end-start) ;
         calcul_flop ("mnblas_dnrm2 ", 2* NB_FOIS, end-start) ;
         if (resd != (double)sqrt(pow(test, 2) * (double) SIZE_VECTOR)) {;
           printf("Erreur mnblas_dnrm2\n");
           exit(1);
         }
       }
       break;
    case 'c':
      printf("\n##################\nTests de mnblas_scnrm2\n##################\n");
       for (test = 0; test < NB_TEST; test ++) {
         printf("\n########### TEST %d############\n", test);
         init_vecteurf(fvector, (float) test);
         start = _rdtsc () ;
         for (i = 0; i < NB_FOIS; i ++) {
            resf = mnblas_scnrm2(SIZE_VECTOR/2, fvector, 1);
         }
         end = _rdtsc () ;
         printf ("La norme euclidienne d'un vecteur: %lld cycles \n", end-start) ;
         calcul_flop ("mnblas_scnrm2 ", 2* NB_FOIS, end-start) ;
         if (resf != (float)sqrt(pow(test, 2) * (float) SIZE_VECTOR)) {;
           printf("Erreur mnblas_scnrm2\n");
           exit(1);
         }
       }
       break;
    case 'z':
      printf("\n##################\nTests de mnblas_dznrm2\n##################\n");
       for (test = 0; test < NB_TEST; test ++) {
         printf("\n########### TEST %d############\n", test);
         init_vecteurd(dvector, (double) test);
         start = _rdtsc () ;
         for (i = 0; i < NB_FOIS; i ++) {
            resd = mnblas_dznrm2(SIZE_VECTOR/2, dvector, 1);
         }
         end = _rdtsc () ;
         printf ("La norme euclidienne d'un vecteur: %lld cycles \n", end-start) ;
         calcul_flop ("mnblas_dznrm2 ", 2* NB_FOIS, end-start) ;
         if (resd != (double)sqrt(pow(test, 2) * (double) SIZE_VECTOR)) {;
           printf("Erreur mnblas_dznrm2\n");
           exit(1);
         }
       }
       break;
    default:
      printf("Mode non reconnu\n");
      exit(1);
 }

  printf("All test passed\n");
  exit (0) ;
}
