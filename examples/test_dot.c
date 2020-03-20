#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"

#include "flop.h"

#define VECSIZE    256000

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;
typedef double vdouble [VECSIZE] ;

vfloat vec1, vec2 ;
vdouble dvec1, dvec2 ;

void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void dvector_init (vdouble V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;

  return ;
}

int main (int argc, char **argv)
{
 unsigned long long start, end ;
 float res ;
 double dres;
 int i ;
 float tmp_fres[2];
 double tmp_dres[2];

 init_flop () ;

 printf("\n########### TEST #############\nFonction mncblas_sdot\n\n");

     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     res = 0.0 ;

     start = _rdtsc () ;
     for (i = 0 ; i < NB_FOIS; i++)
       {
        res = mncblas_sdot (VECSIZE, vec1, 1, vec2, 1) ;
           }
     end = _rdtsc () ;

     printf ("mncblas_sdot : res = %3.2f nombre de cycles: %Ld \n",  res, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE * NB_FOIS, end-start) ;


   printf("\n########### TEST #############\nFonction mncblas_ddot\n\n");

     dvector_init (dvec1, 1.0) ;
     dvector_init (dvec2, 2.0) ;
     dres = 0.0 ;
     start = _rdtsc () ;
     for (i = 0 ; i < NB_FOIS; i++)
      {
        dres = mncblas_ddot (VECSIZE, dvec1, 1, dvec2, 1) ;
      }
     end = _rdtsc () ;

     printf ("mncblas_ddot : res = %3.2f nombre de cycles: %Ld \n", dres, end-start) ;
     calcul_flop ("ddot ", 2 * VECSIZE * NB_FOIS, end-start) ;

   printf("\n########### TEST #############\nFonction mncblas_cdotu_sub\n\n");

     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;

     start = _rdtsc () ;
     for (i = 0 ; i < NB_FOIS; i++)
      {
        mncblas_cdotu_sub (VECSIZE/2, vec1, 1, vec2, 1, tmp_fres) ;
      }
     end = _rdtsc () ;

     printf ("mncblas_cdotu_sub : res = %3.2f + i %3.2f nombre de cycles: %Ld \n", tmp_fres[0], tmp_fres[1], end-start) ;
     calcul_flop ("cdotu_sub", 2 * VECSIZE * NB_FOIS, end-start) ;

   printf("\n########### TEST #############\nFonction mncblas_cdotc_sub\n\n");

      vector_init (vec1, 1.0) ;
      vector_init (vec2, 2.0) ;
      start = _rdtsc () ;
      for (i = 0 ; i < NB_FOIS; i++)
       {
         mncblas_cdotc_sub (VECSIZE/2, vec1, 1, vec2, 1, tmp_fres) ;
       }
      end = _rdtsc () ;

      printf ("mncblas_cdotc_sub : res = %3.2f + i %3.2f nombre de cycles: %Ld \n", tmp_fres[0], tmp_fres[1], end-start) ;
      calcul_flop ("cdotc_sub", 2 * VECSIZE * NB_FOIS, end-start) ;

    printf("\n########### TEST #############\nFonction mncblas_zdotu_sub\n\n");

       dvector_init (dvec1, 1.0) ;
       dvector_init (dvec2, 2.0) ;

       start = _rdtsc () ;
       for (i = 0 ; i < NB_FOIS; i++)
        {
          mncblas_zdotu_sub (VECSIZE/2, dvec1, 1, dvec2, 1, tmp_dres) ;
        }
       end = _rdtsc () ;

       printf ("mncblas_zdotu_sub : res = %3.2f + i %3.2f nombre de cycles: %Ld \n", tmp_dres[0], tmp_dres[1], end-start) ;
       calcul_flop ("zdotu_sub", 2 * VECSIZE * NB_FOIS, end-start) ;


     printf("\n########### TEST #############\nFonction mncblas_zdotc_sub\n\n");

        dvector_init (dvec1, 1.0) ;
        dvector_init (dvec2, 2.0) ;

        start = _rdtsc () ;
        for (i = 0 ; i < NB_FOIS; i++)
         {
           mncblas_zdotc_sub (VECSIZE/2, dvec1, 1, dvec2, 1, tmp_dres) ;
         }
        end = _rdtsc () ;

        printf ("mncblas_zdotc_sub : res = %3.2f + i %3.2f nombre de cycles: %Ld \n", tmp_dres[0], tmp_dres[1], end-start) ;
        calcul_flop ("zdotc_sub",2 * VECSIZE * NB_FOIS, end-start) ;

}
