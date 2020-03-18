#include "mnblas.h"

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const float alpha, 
const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc){

                    if (alpha == 0){
                        for(; i< M * M){
                            C [i] *= beta;
                        }
                    }else{
                        for (int k = 0; k<M ; k++){
                            for (int i = 0; i<M ; i++){
                                C [i + k*M] *= (alpha/beta);
                                for (int j = 0; j<N ; j++){
                                    C [i + k*M] += A [j + k*M] * B [i + j*M]; 
                                }
                                C [i + k*M] *= alpha;
                            }
                        }
                    }
                 }

void mncblas_dgemm(MNCBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, 
const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc){

                    if (alpha_tab [0] == 0){
                        for(int i = 0; i< M * M){
                            C [i] *= beta;
                        }
                    }else{
                        for (int k = 0; k<M ; k++){
                            for (int i = 0; i<M ; i++){
                                C [i + k*M] *= (alpha/beta);
                                for (int j = 0; j<N ; j++){
                                    C [i + k*M] += A [j + k*M] * B [i + j*M]; 
                                }
                                C [i + k*M] *= alpha;
                            }
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

                    if (alpha_tab == [0]){
                        for(int i = 0; i< M * M){
                            C_tab [i] *= beta_tab [0];
                        }
                    }else{
                        for (int k = 0; k<M ; k++){
                            for (int i = 0; i<M ; i++){
                                C_tab [2*i + 2*k*M] *= (alpha_tab [0]/beta_tab [0]);
                                C_tab [2*i + 1 + 2*k*M] *= (alpha_tab [0]/beta_tab [0]);
                                for (int j = 0; j<N ; j++){
                                    C_tab [2*i + 2*k*M] += (A_tab [2*j + 2*k*M] * B_tab [2*i + 2*j*M]) - (A_tab [2*j + 2*k*M + 1] * B_tab [2*i + 2*j*M + 1]); 
                                    C_tab [2*i + 2*k*M + 1] += (A_tab [2*j + 2*k*M ] * B_tab [2*i + 2*j*M + 1]) + (A_tab [2*j + 2*k*M + 1 ] * B_tab [2*i + 2*j*M]); 
                                }
                                C_tab [2*i + 2*k*M] *= alpha_tab [0];
                                C_tab [2*i + 2*k*M + 1] *= alpha_tab [0];
                            }
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

                    if (alpha_tab == [0]){
                        for(int i = 0; i< M * M){
                            C_tab [i] *= beta_tab [0];
                        }
                    }else{
                        for (int k = 0; k<M ; k++){
                            for (int i = 0; i<M ; i++){
                                C_tab [2*i + 2*k*M] *= (alpha_tab [0]/beta_tab [0]);
                                C_tab [2*i + 1 + 2*k*M] *= (alpha_tab [0]/beta_tab [0]);
                                for (int j = 0; j<N ; j++){
                                    C_tab [2*i + 2*k*M] += (A_tab [2*j + 2*k*M] * B_tab [2*i + 2*j*M]) - (A_tab [2*j + 2*k*M + 1] * B_tab [2*i + 2*j*M + 1]); 
                                    C_tab [2*i + 2*k*M + 1] += (A_tab [2*j + 2*k*M ] * B_tab [2*i + 2*j*M + 1]) + (A_tab [2*j + 2*k*M + 1 ] * B_tab [2*i + 2*j*M]); 
                                }
                                C_tab [2*i + 2*k*M] *= alpha_tab [0];
                                C_tab [2*i + 2*k*M + 1] *= alpha_tab [0];
                            }
                        }
                    }
                 }