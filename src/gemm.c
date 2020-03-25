#include "mnblas.h"
#include "complexe2.h"

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const float alpha, 
const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc){

                    register unsigned int k = 0;
                    register unsigned int i;
                    register unsigned int j;
                    register unsigned int K_M;

                    for (; k<M ; k++){
                        K_M = k*M;
                        for (i = 0 ; i<M ; i++) {
                            C [i + K_M] *= (beta/alpha);
                            for (j = 0; j<N ; j++) {
                                C [i + K_M] += A [j + K_M] * B [i + j*M]; 
                            }
                            C [i + K_M] *= alpha;
                        }
                    }
                 }

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, 
const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc){

                    register unsigned int k = 0;
                    register unsigned int i;
                    register unsigned int j;
                    register unsigned int K_M;

                    for (; k<M ; k++){
                        K_M = k*M;
                        for (i = 0 ; i<M ; i++) {
                            C [i + K_M] *= (beta/alpha);
                            for (j = 0; j<N ; j++) {
                                C [i + K_M] += A [j + K_M] * B [i + j*M]; 
                            }
                            C [i + K_M] *= alpha;
                        }
                    }
                 }

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void *alpha, 
const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc){

                    float * alpha_tab = (float *)alpha;
                    float * beta_tab = (float *)beta;
                    complexe_float_t * A_tab = (complexe_float_t *)A;
                    complexe_float_t * B_tab = (complexe_float_t *)B;
                    complexe_float_t * C_tab = (complexe_float_t *)C;

                    register unsigned int k = 0;
                    register unsigned int i;
                    register unsigned int j;
                    register unsigned int K_M;

                    complexe_float_t alpha_c;
                    complexe_float_t beta_c;
                    complexe_float_t beta_div_alpha;
                    complexe_float_t A_B;
                    alpha_c = (complexe_float_t){alpha_tab[0], alpha_tab[1]};
                    beta_c = (complexe_float_t){beta_tab[0], beta_tab[1]};
                    beta_div_alpha = div_complexe_float(beta_c, alpha_c);


                    for (; k<M ; k++) {
                        K_M = k*M;
                        for (i = 0; i<M ; i++) { 
                            C_tab [i + K_M] = mult_complexe_float(C_tab [i + K_M], beta_div_alpha);
                            for (j = 0; j<N ; j++){
                                A_B = mult_complexe_float(A_tab [j + K_M], B_tab[i + j*M]);
                                C_tab[i + K_M] = add_complexe_float(C_tab[i + K_M], A_B);
                            }
                            C_tab[i + K_M] = mult_complexe_float(C_tab[i + K_M], alpha_c);
                        }
                    }
                 }

void mncblas_zgemm (MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void *alpha, 
const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc){

                    complexe_double_t * A_tab = (complexe_double_t *)A;
                    complexe_double_t * B_tab = (complexe_double_t *)B;
                    complexe_double_t * C_tab = (complexe_double_t *)C;

                    register unsigned int k = 0;
                    register unsigned int i;
                    register unsigned int j;
                    register unsigned int K_M;
                    register unsigned int I_K_M;

                    complexe_double_t *alpha_c = (complexe_double_t *)alpha;
                    complexe_double_t *beta_c = (complexe_double_t *)beta;
                    complexe_double_t beta_div_alpha = div_complexe_double(*beta_c, *alpha_c);

                    for (; k<M ; k++) {
                        K_M = k*M;
                        for (i = 0; i<M ; i++) { 
                            I_K_M = i + K_M;
                            C_tab [I_K_M] = mult_complexe_double(C_tab [I_K_M], beta_div_alpha);
                            for (j = 0; j<N ; j++){
                                C_tab[I_K_M] = add_complexe_double(C_tab[I_K_M], mult_complexe_double(A_tab [j + K_M], B_tab[i + j*M]));
                            }
                            C_tab[I_K_M] = mult_complexe_double(C_tab[I_K_M], *alpha_c);
                        }
                    }
                 }