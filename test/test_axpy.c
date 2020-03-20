#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"

#include "flop.h"

#define NB_FOIS 4000
#define VECSIZE 4

void fvector_init(float *V)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = i;

  return;
}

void dvector_init(double *V)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = i;

  return;
}

void fvector_print(float *V)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("%f ", V[i]);
  printf("\n");

  return;
}

void dvector_print(double *V)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("%f ", V[i]);
  printf("\n");

  return;
}

int max(int a, int b)
{
  if (a > b)
  {
    return a;
  }
  return b;
}

int main(int argc, char **argv)
{

  if (argc != 2) {
    printf("Utilisation: ./test_iamx mode\n");
    exit(1);
  }
  char mode = argv[1][0];

  float fX[VECSIZE];
  float fY[VECSIZE];
  double dX[VECSIZE];
  double dY[VECSIZE];
  double dalpha[2];
  float falpha[2];

  falpha[0] = 2.0;
  falpha[1] = 2.0;
  dalpha[0] = 2.0;
  dalpha[1] = 2.0;
  float alpha = 2.0;

  int incX = 1;
  int incY = 1;
  int nb_boucles;

  fvector_init(fX);
  fvector_init(fY);
  dvector_init(dX);
  dvector_init(dY);

  unsigned long long int start, end;
  int i;

  init_flop();

  switch (mode) {
    case 's':
      printf("\n########### TEST #############\nFonction mnblas_saxpy\n\n");
      printf("%f * X + Y\n", alpha);
      printf("X :\n");
      fvector_print(fX);
      printf("Y :\n");
      fvector_print(fY);
      mnblas_saxpy(VECSIZE, alpha, fX, incX, fY, incY);
      printf("\n");
      printf("Resultat :\n");
      fvector_print(fY);
      fvector_init(fY);
      mnblas_saxpy(VECSIZE, alpha, fX, incX + 1, fY, incY);
      printf("Resultat avec incX = 2 :\n");
      fvector_print(fY);
      fvector_init(fY);
      start = _rdtsc();
      for (i = 0; i < NB_FOIS; i++)
      {
        mnblas_saxpy(VECSIZE, alpha, fX, incX, fY, incY);
      }
      end = _rdtsc();
      printf("mnblas_saxpy %d : nombre de cycles: %Ld \n", i, end - start);
      nb_boucles = VECSIZE / max(incX, incY);
      calcul_flop("mnblas_saxpy ", nb_boucles * 2, end - start);
      printf("\n");
      break;
    case 'd':
      printf("\n########### TEST #############\nFonction mnblas_daxpy\n\n");
      printf("%lf * X + Y\n", alpha);
      printf("X :\n");
      dvector_print(dX);
      printf("Y :\n");
      dvector_print(dY);
      mnblas_daxpy(VECSIZE, alpha, dX, incX, dY, incY);
      printf("\n");
      printf("Resultat :\n");
      dvector_print(dY);
      dvector_init(dY);
      mnblas_daxpy(VECSIZE, alpha, dX, incX + 1, dY, incY);
      printf("Resultat avec incX = 2 :\n");
      dvector_print(dY);
      dvector_init(dY);
      start = _rdtsc();
      for (i = 0; i < NB_FOIS; i++)
      {
        mnblas_daxpy(VECSIZE, alpha, dX, incX, dY, incY);
      }
      end = _rdtsc();
      printf("mnblas_daxpy %d : nombre de cycles: %Ld \n", i, end - start);
      nb_boucles = VECSIZE / max(incX, incY);
      calcul_flop("mnblas_daxpy ", nb_boucles * 2, end - start);
      printf("\n");
      break;
    case 'c':
      printf("\n########### TEST #############\nFonction mnblas_caxpy\n\n");
      fvector_init(fY);
      printf("%f * X + Y\n", falpha[0]);
      printf("X :\n");
      fvector_print(fX);
      printf("Y :\n");
      fvector_print(fY);
      mnblas_caxpy(VECSIZE, falpha, fX, incX, fY, incY);
      printf("\n");
      printf("Resultat :\n");
      fvector_print(fY);

      fvector_init(fY);
      mnblas_caxpy(VECSIZE, falpha, fX, incX + 1, fY, incY);
      printf("Resultat avec incX = 2 :\n");
      fvector_print(fY);
      fvector_init(fY);
      start = _rdtsc();
      for (i = 0; i < NB_FOIS; i++)
      {
        mnblas_caxpy(VECSIZE, falpha, fX, incX, fY, incY);
      }
      end = _rdtsc();
      printf("mnblas_caxpy %d : nombre de cycles: %Ld \n", i, end - start);
      nb_boucles = VECSIZE / max(incX, incY);
      calcul_flop("mnblas_caxpy ", nb_boucles * 4, end - start);
      printf("\n");
      break;
    case 'z':
      printf("\n########### TEST #############\nFonction mnblas_zaxpy\n\n");
      dvector_init(dY);
      printf("%lf * X + Y\n", dalpha[0]);
      printf("X :\n");
      dvector_print(dX);
      printf("Y :\n");
      dvector_print(dY);
      mnblas_zaxpy(VECSIZE, dalpha, dX, incX, dY, incY);
      printf("\n");
      printf("Resultat :\n");
      dvector_print(dY);
      dvector_init(dY);
      mnblas_zaxpy(VECSIZE, dalpha, dX, incX + 1, dY, incY);
      printf("Resultat avec incX = 2 :\n");
      dvector_print(dY);
      dvector_init(dY);
      start = _rdtsc();
      for (i = 0; i < NB_FOIS; i++)
      {
        mnblas_zaxpy(VECSIZE, dalpha, dX, incX, dY, incY);
      }
      end = _rdtsc();
      printf("mnblas_zaxpy %d : nombre de cycles: %Ld \n", i, end - start);
      nb_boucles = VECSIZE / max(incX, incY);
      calcul_flop("mnblas_zaxpy ", nb_boucles * 4, end - start);
      printf("\n");
      break;
    default:
      printf("Mode non reconnu\n");
      exit(1);
  }
  exit(0);
}
