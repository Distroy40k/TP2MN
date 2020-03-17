#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"

#include "flop.h"

#define NB_FOIS 419430

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
  static int taille_tab = 100;
 float *f1 = (float *) malloc(sizeof(float) * taille_tab);
 float *f2 = (float *) malloc(sizeof(float) * taille_tab);
 double *d1 = (double *) malloc(sizeof(double) * taille_tab);
 double *d2 = (double *) malloc(sizeof(double) * taille_tab);

 unsigned long long int start, end ;
 int i;

 printf("\n");
 printf("Tests des fonctions de copy\n");
 printf("\n");

 printf("\n");
 printf("Tests de mncblas_scopy\n");
 printf("\n");

 for (int var = 0; var < taille_tab; var += 1) {
   f1[var] = (float) var;
 }
 init_flop () ;
 start = _rdtsc () ;
 for (i = 0; i < NB_FOIS; i ++) {
    mncblas_scopy(taille_tab, f1, 1, f2, 1);
 }
 end = _rdtsc () ;
 printf ("La copie d'un tableau de vecteur de float simple précision, de 0 à %d: %lld cycles \n", taille_tab, end-start) ;
 calcul_octet ("mncblas_scopy ", NB_FOIS*4*taille_tab, end-start) ;
 if (is_equal_float(f1,f2, taille_tab) != 1) {
   printf("Erreur copie mncblas_scopy\n");
   exit(1);
 }

 printf("\n");
 printf("Tests de mncblas_dcopy\n");
 printf("\n");

 for (int var = 0; var < taille_tab; var += 1) {
   d1[var] = (double) var;
 }
 start = _rdtsc () ;
 for (i = 0; i < NB_FOIS; i ++) {
    mncblas_dcopy(taille_tab, d1, 1, d2, 1);
 }
 end = _rdtsc () ;
 printf ("La copie d'un tableau de vecteur de double simple précision, de 0 à %d: %lld cycles \n", taille_tab, end-start) ;
 calcul_octet ("mncblas_dcopy ", NB_FOIS*4, end-start) ;
 if (is_equal_double(d1,d2, taille_tab) != 1) {
   printf("Erreur copie mncblas_dcopy\n");
   exit(1);
 }



 printf("\n");
 printf("Tests de mncblas_ccopy\n");
 printf("\n");

 f2 = (float *) malloc(sizeof(float) * taille_tab);
 init_flop () ;
 start = _rdtsc () ;
 for (i = 0; i < NB_FOIS; i ++) {
   mncblas_ccopy(taille_tab / 2, f1, 1, f2, 1);
 }
 end = _rdtsc () ;
 printf ("La copie d'un tableau de vecteur de float double précision, de 0 à %d:  %lld cycles \n", taille_tab, end-start) ;
 calcul_octet ("mncblas_ccopy ", NB_FOIS*4, end-start) ;
 if (is_equal_float(f1,f2, taille_tab) != 1) {
   printf("Erreur copie mncblas_ccopy\n");
   exit(1);
 }

 printf("\n");
 printf("Tests de mncblas_zcopy\n");
 printf("\n");

 d2 = (double *) malloc(sizeof(double) * taille_tab);
 start = _rdtsc () ;
 for (i = 0; i < NB_FOIS; i ++) {
 mncblas_zcopy(taille_tab / 2, d1, 1, d2, 1);
 }
 end = _rdtsc () ;
 printf ("La copie d'un tableau de vecteur de double double précision, de 0 à %d: %lld cycles \n", taille_tab, end-start) ;
 calcul_octet ("mncblas_zcopy ", NB_FOIS*4, end-start) ;
 if (is_equal_double(d1,d2, taille_tab) != 1) {
   printf("Erreur copie mncblas_zcopy\n");
   exit(1);
 }


 free(f1);
 free(f2);
 free(d1);
 free(d2);



  exit (0) ;
}
