#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "flop.h"

#define NB_TEST 10
#define NB_FOIS 45125
#define SIZE_VECTOR 452

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
int is_equal_float(float *f1, float *f2) {
  for (int i = 0; i < SIZE_VECTOR; i++) {
    if (f1[i] != f2[i]) {
      return 0;
    }
  }
  return 1;
}
int is_equal_double(double *d1, double *d2) {
  for (int i = 0; i < SIZE_VECTOR; i++) {
    if (d1[i] != d2[i]) {
      return 0;
    }
  }
  return 1;
}

int main (int argc, char **argv)
{
 float f1[SIZE_VECTOR];
 float f2[SIZE_VECTOR];
 float f1p[SIZE_VECTOR];
 float f2p[SIZE_VECTOR];
 double d1[SIZE_VECTOR];
 double d2[SIZE_VECTOR];
 double d1p[SIZE_VECTOR];
 double d2p[SIZE_VECTOR];

 unsigned long long int start, end ;
 int i;

 printf("\nTests des fonctions de swap\n\n");

 printf("\nTests de mncblas_sswap\n\n");

 init_vecteurf(f1,  (float) 1);
 init_vecteurf(f1p, (float) 1);
 init_vecteurf(f2,  (float) 2);
 init_vecteurf(f2p, (float) 2);
 init_flop () ;
 start = _rdtsc () ;
 for (i = 0; i < NB_FOIS; i ++) {
    mncblas_sswap(SIZE_VECTOR, f1, 1, f2, 1);
 }
 end = _rdtsc () ;
 printf ("Copy avec sswap: %lld cycles \n", end-start) ;
 calcul_octet ("mncblas_sswap ", NB_FOIS*6*SIZE_VECTOR, end-start) ;
 if (is_equal_float(f1,f2p) != 1) {
   printf("Erreur mncblas_sswap\n");
   exit(1);
 }
 if (is_equal_float(f2,f1p) != 1) {
   printf("Erreur mncblas_sswap\n");
   exit(1);
 }

 printf("\nTests de mncblas_dswap\n\n");

 init_vecteurd(d1, 1);
 init_vecteurd(d1p,1);
 init_vecteurd(d2, 2);
 init_vecteurd(d2p,2);
 init_flop () ;
 start = _rdtsc () ;
 for (i = 0; i < NB_FOIS; i ++) {
    mncblas_dswap(SIZE_VECTOR, d1, 1, d2, 1);
 }
 end = _rdtsc () ;
 printf ("Copy avec dswap: %lld cycles \n", end-start) ;
 calcul_octet ("mncblas_dswap ", NB_FOIS*6*SIZE_VECTOR, end-start) ;
 if (is_equal_double(d1,d2p) != 1) {
   printf("Erreur mncblas_dswap\n");
   exit(1);
 }
 if (is_equal_double(d2,d1p) != 1) {
   printf("Erreur mncblas_dswap\n");
   exit(1);
 }

 printf("\nTests de mncblas_cswap\n\n");

 init_vecteurf(f1, 1);
 init_vecteurf(f1p,1);
 init_vecteurf(f2, 2);
 init_vecteurf(f2p,2);
 init_flop () ;
 start = _rdtsc () ;
 for (i = 0; i < NB_FOIS; i ++) {
    mncblas_cswap(SIZE_VECTOR/2, f1, 1, f2, 1);
 }
 end = _rdtsc () ;
 printf ("Copy avec cswap: %lld cycles \n", end-start) ;
 calcul_octet ("mncblas_sswap ", NB_FOIS*6*SIZE_VECTOR, end-start) ;
 if (is_equal_float(f1,f2p) != 1) {
   printf("Erreur mncblas_cswap\n");
   exit(1);
 }
 if (is_equal_float(f2,f1p) != 1) {
   printf("Erreur mncblas_cswap\n");
   exit(1);
 }

 printf("\nTests de mncblas_zswap\n\n");
 init_vecteurd(d1, 1);
 init_vecteurd(d1p,1);
 init_vecteurd(d2, 2);
 init_vecteurd(d2p,2);
 init_flop () ;
 start = _rdtsc () ;
 for (i = 0; i < NB_FOIS; i ++) {
    mncblas_zswap(SIZE_VECTOR/2, d1, 1, d2, 1);
 }
 end = _rdtsc () ;
 printf ("Copy avec zswap: %lld cycles \n", end-start) ;
 calcul_octet ("mncblas_zswap ", NB_FOIS*6*SIZE_VECTOR, end-start) ;
 if (is_equal_double(d1,d2p) != 1) {
   printf("Erreur mncblas_dswap\n");
   exit(1);
 }
 if (is_equal_double(d2,d1p) != 1) {
   printf("Erreur mncblas_dswap\n");
   exit(1);
 }

}
