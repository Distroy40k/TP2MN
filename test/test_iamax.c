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
 CBLAS_INDEX higher_indice;
 CBLAS_INDEX res;
 init_flop () ;

 printf("\n");
 printf("Tests des fonctions de iamax\n");
 printf("\n");

 printf("##########################################\nPour tester ces fonctions, on cré un vecteur dont toutes les coordonées sont égales a 0, et on set alétatoirement une coordonées à 1, et on verifie que le résultat est bon\n####################################\n");
  printf("\n##################\nTests de mnblas_isamax\n##################\n");
   for (test = 0; test < NB_TEST; test ++) {
     printf("\n########### TEST %d ############\n", test);
      higher_indice = rand() % SIZE_VECTOR;
     init_vecteurf(fvector, 0);
     fvector[higher_indice] = 1;
     start = _rdtsc () ;
     for (i = 0; i < NB_FOIS; i ++) {
        res = mnblas_isamax(SIZE_VECTOR, fvector, 1);
     }
     end = _rdtsc () ;
     calcul_flop ("mnblas_isamax ", NB_FOIS, end-start) ;
     if (res != higher_indice) {
       printf("Erreur copie mnblas_isamax\n");
       exit(1);
     }
   }
 printf("\n##################\nTests de mnblas_idamax\n##################\n");
  for (test = 0; test < NB_TEST; test ++) {
    printf("\n########### TEST %d ############\n", test);
     higher_indice = rand() % SIZE_VECTOR;
    init_vecteurd(dvector, 0);
    dvector[higher_indice] = 1;
    start = _rdtsc () ;
    for (i = 0; i < NB_FOIS; i ++) {
       res = mnblas_idamax(SIZE_VECTOR, dvector, 1);
    }
    end = _rdtsc () ;
    calcul_flop ("mnblas_idamax ", NB_FOIS, end-start) ;
    if (res != higher_indice) {
      printf("Erreur copie mnblas_idamax\n");
      exit(1);
    }
  }
  printf("\n##################\nTests de mnblas_icamax\n##################\n");
   for (test = 0; test < NB_TEST; test ++) {
     printf("\n########### TEST %d ############\n", test);
      higher_indice = rand() % SIZE_VECTOR;
     init_vecteurf(fvector, 0);
     fvector[higher_indice] = 1;
     start = _rdtsc () ;
     for (i = 0; i < NB_FOIS; i ++) {
        res = mnblas_icamax(SIZE_VECTOR/2, fvector, 1);
     }
     end = _rdtsc () ;
     calcul_flop ("mnblas_icamax ", 3 * NB_FOIS, end-start) ;
     higher_indice = higher_indice / 2;
     if (res != higher_indice) {
       printf("Erreur copie mnblas_icamax\n");
       exit(1);
     }
   }

   printf("\n##################\nTests de mnblas_izamax\n##################\n");
    for (test = 0; test < NB_TEST; test ++) {
      printf("\n########### TEST %d ############\n", test);
       higher_indice = rand() % SIZE_VECTOR;
      init_vecteurd(dvector, 0);
      dvector[higher_indice] = 1;
      start = _rdtsc () ;
      for (i = 0; i < NB_FOIS; i ++) {
         res = mnblas_izamax(SIZE_VECTOR/2, dvector, 1);
      }
      end = _rdtsc () ;
      calcul_flop ("mnblas_izamax ", 3 * NB_FOIS, end-start) ;
      higher_indice = higher_indice / 2;
      if (res != higher_indice) {
        printf("Erreur copie mnblas_izamax\n");
        exit(1);
      }
    }

  printf("All test passed\n");

  exit (0) ;
}
