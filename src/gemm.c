#include "mnblas.h"

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const float alpha, 
const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc){

                    for (int k = 0; k<M ; k++){
                        for (int i = 0; i<M ; i++){
                            C [i + k*M] *= (beta/alpha);
                            for (int j = 0; j<N ; j++){
                                C [i + k*M] += A [j + k*M] * B [i + j*M]; 
                            }
                            C [i + k*M] *= alpha;
                        }
                    }
                 }

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, 
const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc){

                    for (int k = 0; k<M ; k++){
                        for (int i = 0; i<M ; i++){
                            C [i + k*M] *= (beta/alpha);
                            for (int j = 0; j<N ; j++){
                                C [i + k*M] += A [j + k*M] * B [i + j*M]; 
                            }
                            C [i + k*M] *= alpha;
                        }
                    }
                 }

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void *alpha, 
const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc){

                    float * alpha_tab = (float *)alpha;
                    float * beta_tab = (float *)beta;
                    float * A_tab = (float *)A;
                    float * B_tab = (float *)B;
                    float * C_tab = (float *)C;

                    float tmp_C1[2];
                    float tmp_C2;

                    for (int k = 0; k<M ; k++){
                        for (int i = 0; i<M ; i++){ 
                            tmp_C1[0] = beta_tab[0] * C_tab[2*i + 2*k*M] - (beta_tab[1] * C_tab[2*i + 1 + 2*k*M]);
                            tmp_C1[1] = beta_tab[1] * C_tab[2*i + 2*k*M] + (beta_tab[0] * C_tab[2*i + 1 + 2*k*M]);
                            C_tab [2*i + 2*k*M] = 0;
                            C_tab [2*i + 1 + 2*k*M] = 0; 
                            for (int j = 0; j<N ; j++){
                                C_tab [2*i + 2*k*M] += (A_tab [2*j + 2*k*M] * B_tab [2*i + 2*j*M]) - (A_tab [2*j + 2*k*M + 1] * B_tab [2*i + 2*j*M + 1]); 
                                C_tab [2*i + 2*k*M + 1] += (A_tab [2*j + 2*k*M ] * B_tab [2*i + 2*j*M + 1]) + (A_tab [2*j + 2*k*M + 1 ] * B_tab [2*i + 2*j*M]); 
                            }
                            tmp_C2 = C_tab [2*i + 2*k*M];
                            C_tab [2*i + 2*k*M] = alpha_tab [0] * C_tab [2*i + 2*k*M] - alpha_tab [1] * C_tab [2*i + 1 + 2*k*M] + tmp_C1[0];
                            C_tab [2*i + 2*k*M + 1] = alpha_tab [1] * tmp_C2 + alpha_tab[0] * C_tab [2*i + 1 + 2*k*M ] + tmp_C1[1];
                        }
                    }
                 }

void mncblas_zgemm (MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void *alpha, 
const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc){

                    double * alpha_tab = (double *)alpha;
                    double * beta_tab = (double *)beta;
                    double * A_tab = (double *)A;
                    double * B_tab = (double *)B;
                    double * C_tab = (double *)C;

                    double tmp_C1[2];
                    double tmp_C2;

                    for (int k = 0; k<M ; k++){
                        for (int i = 0; i<M ; i++){ 
                            tmp_C1[0] = beta_tab[0] * C_tab[2*i + 2*k*M] - (beta_tab[1] * C_tab[2*i + 1 + 2*k*M]);
                            tmp_C1[1] = beta_tab[1] * C_tab[2*i + 2*k*M] + (beta_tab[0] * C_tab[2*i + 1 + 2*k*M]);
                            C_tab [2*i + 2*k*M] = 0;
                            C_tab [2*i + 1 + 2*k*M] = 0; 
                            for (int j = 0; j<N ; j++){
                                C_tab [2*i + 2*k*M] += (A_tab [2*j + 2*k*M] * B_tab [2*i + 2*j*M]) - (A_tab [2*j + 2*k*M + 1] * B_tab [2*i + 2*j*M + 1]); 
                                C_tab [2*i + 2*k*M + 1] += (A_tab [2*j + 2*k*M ] * B_tab [2*i + 2*j*M + 1]) + (A_tab [2*j + 2*k*M + 1 ] * B_tab [2*i + 2*j*M]); 
                            }
                            tmp_C2 = C_tab [2*i + 2*k*M];
                            C_tab [2*i + 2*k*M] = alpha_tab [0] * C_tab [2*i + 2*k*M] - alpha_tab [1] * C_tab [2*i + 1 + 2*k*M] + tmp_C1[0];
                            C_tab [2*i + 2*k*M + 1] = alpha_tab [1] * tmp_C2 + alpha_tab[0] * C_tab [2*i + 1 + 2*k*M ] + tmp_C1[1];
                        }
                    }
                 }