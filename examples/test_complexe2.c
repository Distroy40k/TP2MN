#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#define NB_FOIS 4194304
#include "flop.h"

int main(int argc, char **argv)
{
  //Pour vérifier les résultats, vous pouvez décommenter les parties situées avant les boucles 
  //d'évaluation des fonctions

  complexe_float_t c1 = {10.0, 7.0};
  complexe_float_t c2 = {25.0, 32.0};
  volatile complexe_float_t c3;

  complexe_double_t cd1;
  complexe_double_t cd2;
  volatile complexe_double_t cd3;

  unsigned long long int start, end;
  int i;

  // Partie addition

  printf("C1 : %f + i x %f , C2 : %f + i x %f\n", c1.real, c1.imaginary, c2.real, c2.imaginary);
  printf("\n");
  printf("Tests des fonctions d'addition de complexes\n");
  printf("\n");
  printf("Float\n");
  printf("\n");

  init_flop();

  printf("C1 : %f + i x %f , C2 : %f + i x %f\n", c1.real, c1.imaginary, c2.real, c2.imaginary);
  c3 = add_complexe_float(c1, c2);
  printf("Evaluation des performances : \n");

  start = _rdtsc();
  for (i = 0; i < NB_FOIS; i++)
  {
    c3 = add_complexe_float(c1, c2);
  }
  end = _rdtsc();

  printf("Nombre de cycles :  %lld cycles \n", end - start);
  calcul_flop("Addition de complexes Float", NB_FOIS * 2, end - start);
  printf("Resultat : %f + i x %f\n", c3.real, c3.imaginary);


  cd1 = (complexe_double_t){10.0, 7.0};
  cd2 = (complexe_double_t){25.0, 32.0};


  printf("\n");
  printf("Double\n");
  printf("\n");

  cd3 = add_complexe_double(cd1, cd2);
  printf("Evaluation des performances : \n");

  start = _rdtsc();
  for (i = 0; i < NB_FOIS; i++)
  {
    cd3 = add_complexe_double(cd1, cd2);
  }
  end = _rdtsc();


  printf("Nombre de cycles :  %lld cycles \n", end - start);
  calcul_flop("Addition de complexes Double", NB_FOIS * 2, end - start);
  printf("Resultat : %f + i x %f\n", cd3.real, cd3.imaginary);

  // Partie multiplication

  printf("\n");
  printf("Tests des fonctions de multiplication de complexes\n");
  printf("\n");
  printf("Float\n");
  printf("\n");

  c3 = mult_complexe_float(c1, c2);
  printf("Evaluation des performances : \n");

  start = _rdtsc();
  for (i = 0; i < NB_FOIS; i++)
  {
    c3 = mult_complexe_float(c1, c2);
  }
  end = _rdtsc();

  printf("Nombre de cycles :  %lld cycles \n", end - start);
  calcul_flop("Multiplication de complexes Float", NB_FOIS * 6, end - start);
  printf("Resultat : %f + i x %f\n", c3.real, c3.imaginary);

  cd1 = (complexe_double_t){10.0, 7.0};
  cd2 = (complexe_double_t){25.0, 32.0};

  printf("\n");
  printf("Double\n");
  printf("\n");

  cd3 = mult_complexe_double(cd1, cd2);
  printf("Evaluation des performances : \n");

  start = _rdtsc();

  for (i = 0; i < NB_FOIS; i++)
  {
    cd3 = mult_complexe_double(cd1, cd2);
  }

  end = _rdtsc();


  printf("Nombre de cycles :  %lld cycles \n", end - start);
  calcul_flop("Addition de complexes Double", NB_FOIS * 6, end - start);
  printf("Resultat : %f + i x %f\n", cd3.real, cd3.imaginary);

  // Partie division

  printf("\n");
  printf("Tests des fonctions de division de complexes\n");
  printf("\n");
  printf("Float\n");
  printf("\n");

  c3 = div_complexe_float(c1, c2);
  printf("Evaluation des performances : \n");

  start = _rdtsc();
  for (i = 0; i < NB_FOIS; i++)
  {
    c3 = div_complexe_float(c1, c2);
  }
  end = _rdtsc();

  printf("Nombre de cycles :  %lld cycles \n", end - start);
  calcul_flop("Division de complexes Float", NB_FOIS * 11, end - start);
  printf("Resultat : %f + i x %f\n", c3.real, c3.imaginary);

  cd1 = (complexe_double_t){10.0, 7.0};
  cd2 = (complexe_double_t){25.0, 32.0};

  printf("\n");
  printf("Double\n");
  printf("\n");

  cd3 = div_complexe_double(cd1, cd2);
  printf("Evaluation des performances : \n");

  start = _rdtsc();

  for (i = 0; i < NB_FOIS; i++)
  {
    cd3 = div_complexe_double(cd1, cd2);
  }

  end = _rdtsc();


  printf("Nombre de cycles :  %lld cycles \n", end - start);
  calcul_flop("Division de complexes Double", NB_FOIS * 11, end - start);
  printf("Resultat : %f + i x %f\n", cd3.real, cd3.imaginary);

  exit(0);
}
