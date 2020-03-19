#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    65536

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

     printf ("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;


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

     printf ("mncblas_ddot %d : res = %3.2f nombre de cycles: %Ld \n", i, dres, end-start) ;
     calcul_flop ("ddot ", 2 * VECSIZE, end-start) ;

   printf("\n########### TEST #############\nFonction mncblas_cdotu_sub\n\n");

     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     float tmp_res[2];

     start = _rdtsc () ;
     for (i = 0 ; i < NB_FOIS; i++)
      {
        mncblas_cdotu_sub (VECSIZE/2, vec1, 1, vec2, 1, tmp_res) ;
      }
     end = _rdtsc () ;

     printf ("mncblas_cdotu_sub %d : res = %3.2f + i %3.2f nombre de cycles: %Ld \n", i, tmp_res[0], tmp_res[1], end-start) ;
     calcul_flop ("cdotu_sub", 4 * VECSIZE, end-start) ;

   printf("\n########### TEST #############\nFonction mncblas_cdotc_sub\n\n");

      vector_init (vec1, 1.0) ;
      vector_init (vec2, 2.0) ;
      float tmp_res[2];

      start = _rdtsc () ;
      for (i = 0 ; i < NB_FOIS; i++)
       {
         mncblas_cdotc_sub (VECSIZE/2, vec1, 1, vec2, 1, tmp_res) ;
       }
      end = _rdtsc () ;

      printf ("mncblas_cdotc_sub %d : res = %3.2f + i %3.2f nombre de cycles: %Ld \n", i, tmp_res[0], tmp_res[1], end-start) ;
      calcul_flop ("cdotc_sub", 4 * VECSIZE, end-start) ;

    printf("\n########### TEST #############\nFonction mncblas_zdotu_sub\n\n");

       dvector_init (dvec1, 1.0) ;
       dvector_init (dvec2, 2.0) ;
       double tmp_res[2];

       start = _rdtsc () ;
       for (i = 0 ; i < NB_FOIS; i++)
        {
          mncblas_zdotu_sub (VECSIZE/2, dvec1, 1, dvec2, 1, tmp_res) ;
        }
       end = _rdtsc () ;

       printf ("mncblas_zdotu_sub %d : res = %3.2f + i %3.2f nombre de cycles: %Ld \n", i, tmp_res[0], tmp_res[1], end-start) ;
       calcul_flop ("zdotu_sub", 4 * VECSIZE, end-start) ;


     printf("\n########### TEST #############\nFonction mncblas_zdotc_sub\n\n");

        dvector_init (dvec1, 1.0) ;
        dvector_init (dvec2, 2.0) ;
        double tmp_res[2];

        start = _rdtsc () ;
        for (i = 0 ; i < NB_FOIS; i++)
         {
           mncblas_zdotc_sub (VECSIZE/2, dvec1, 1, dvec2, 1, tmp_res) ;
         }
        end = _rdtsc () ;

        printf ("mncblas_zdotc_sub %d : res = %3.2f + i %3.2f nombre de cycles: %Ld \n", i, tmp_res[0], tmp_res[1], end-start) ;
        calcul_flop ("zdotc_sub", 4 * VECSIZE, end-start) ;

}
