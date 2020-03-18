#include <stdio.h>

void remplir_matrice (int n,float x,FILE *f) {
    for(int i =0; i<n*n; i++){
        fprintf(f,"%f\n",x);
    }
}
void remplir_vecteur (int n,float x,FILE *f) {
    for(int i =0; i<n; i++){
        fprintf(f,"%f\n",x);
    }
}

int main()  {
  FILE *f = fopen("m_2_2_0.txt", "w");
  remplir_matrice(2, 0, f);
  f = fopen("m_2_2_1.txt", "w");
  remplir_matrice(2, 1, f);
  f = fopen("m_256_256_1.txt", "w");
  remplir_matrice(256, 1, f);
  f = fopen("m_256_256_0.txt", "w");
  remplir_matrice(256, 0, f);
  f = fopen("m_128_128_1.txt", "w");
  remplir_matrice(256, 1, f);
  f = fopen("m_128_128_0.txt", "w");
  remplir_matrice(256, 0, f);

  f = fopen("x_2_0.txt", "w");
  remplir_vecteur(2,0,f);
  f = fopen("x_2_1.txt", "w");
  remplir_vecteur(2,1,f);
  f = fopen("x_256_0.txt", "w");
  remplir_vecteur(256,0,f);
  f = fopen("x_256_1.txt", "w");
  remplir_vecteur(256,1,f);
  f = fopen("x_128_0.txt", "w");
  remplir_vecteur(128,0,f);
  f = fopen("x_128_1.txt", "w");
  remplir_vecteur(128,1,f);

  f = fopen("y_2_0.txt", "w");
  remplir_vecteur(2,0,f);
  f = fopen("y_2_1.txt", "w");
  remplir_vecteur(2,1,f);
  f = fopen("y_256_0.txt", "w");
  remplir_vecteur(256,0,f);
  f = fopen("y_256_1.txt", "w");
  remplir_vecteur(256,1,f);
  f = fopen("y_128_0.txt", "w");
  remplir_vecteur(128,0,f);
  f = fopen("y_128_1.txt", "w");
  remplir_vecteur(128,1,f);
}
