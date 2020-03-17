#include <stdio.h>
#include <x86intrin.h>
#include <time.h>
#include <stdlib.h>

#include "mnblas.h"

#include "flop.h"

#define NB_TEST 10
#define NB_FOIS 45126
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
  srand(time(NULL));

 float fvector[SIZE_VECTOR];
 double dvector[SIZE_VECTOR];

 unsigned long long int start, end ;
 int i;
 int test;
 CBLAS_INDEX lower_indice;
 CBLAS_INDEX res;
 init_flop () ;

 printf("\n");
 printf("Tests des fonctions de iamin\n");
 printf("\n");

 printf("##########################################\nPour tester ces fonctions, on cré un vecteur dont toutes les coordonées sont égales a 1, et on set alétatoirement une coordonées à 0, et on verifie que le résultat est bon\n####################################\n");
  printf("\n##################\nTests de mnblas_isamin\n##################\n");
   for (test = 0; test < NB_TEST; test ++) {
     printf("\n########### TEST %d ############\n", test);
      lower_indice = rand() % SIZE_VECTOR;
     init_vecteurf(fvector, 1.0);
     fvector[lower_indice] = 0.0;
     start = _rdtsc () ;
     for (i = 0; i < NB_FOIS; i ++) {
        res = mnblas_isamin(SIZE_VECTOR, fvector, 1);
     }
     end = _rdtsc () ;
     calcul_flop ("mnblas_isamin ", NB_FOIS, end-start) ;
     if (res != lower_indice) {
       printf("Erreur copie mnblas_isamin\n");
       exit(1);
     }
   }
 printf("\n##################\nTests de mnblas_idamin\n##################\n");
  for (test = 0; test < NB_TEST; test ++) {
    printf("\n########### TEST %d ############\n", test);
     lower_indice = rand() % SIZE_VECTOR;
    init_vecteurd(dvector, 1.0);
    dvector[lower_indice] = 0.0;
    start = _rdtsc () ;
    for (i = 0; i < NB_FOIS; i ++) {
       res = mnblas_idamin(SIZE_VECTOR, dvector, 1);
    }
    end = _rdtsc () ;
    calcul_flop ("mnblas_idamin ", NB_FOIS, end-start) ;
    if (res != lower_indice) {
      printf("Erreur copie mnblas_idamin\n");
      exit(1);
    }
  }
  printf("\n##################\nTests de mnblas_icamin\n##################\n");
   for (test = 0; test < NB_TEST; test ++) {
     printf("\n########### TEST %d ############\n", test);
      lower_indice = rand() % SIZE_VECTOR;
     init_vecteurf(fvector, 1.0);
     fvector[lower_indice] = 0.0;
     start = _rdtsc () ;
     for (i = 0; i < NB_FOIS; i ++) {
        res = mnblas_icamin(SIZE_VECTOR/2, fvector, 1);
     }
     end = _rdtsc () ;
     calcul_flop ("mnblas_icamin ", 3 * NB_FOIS, end-start) ;
     lower_indice = lower_indice / 2;
     if (res != lower_indice) {
       printf("Erreur copie mnblas_icamin\n");
       exit(1);
     }
   }

   printf("\n##################\nTests de mnblas_izamin\n##################\n");
    for (test = 0; test < NB_TEST; test ++) {
      printf("\n########### TEST %d ############\n", test);
       lower_indice = rand() % SIZE_VECTOR;
      init_vecteurd(dvector, 1.0);
      dvector[lower_indice] = 0.0;
      start = _rdtsc () ;
      for (i = 0; i < NB_FOIS; i ++) {
         res = mnblas_izamin(SIZE_VECTOR/2, dvector, 1);
      }
      end = _rdtsc () ;
      calcul_flop ("mnblas_izamin ", 3 * NB_FOIS, end-start) ;
      lower_indice = lower_indice / 2;
      if (res != lower_indice) {
        printf("Erreur copie mnblas_izamin\n");
        exit(1);
      }
    }

  printf("All test passed\n");

  exit (0) ;
}
