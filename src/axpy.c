#include "mnblas.h"

void mnblas_saxpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY){
    int idx_Y = 0;
    int idx_X = 0;
    while ((idx_Y < N) && (idx_X < N)){
        Y [idx_Y] = alpha * X [idx_X] + Y [idx_Y];
        idx_X += incX;
        idx_Y += incY;
    }
}

void mnblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY){
    int idx_Y = 0;
    int idx_X = 0;
    while ((idx_Y < N) && (idx_X < N)){
        Y [idx_Y] = alpha * X [idx_X] + Y [idx_Y];
        idx_X += incX;
        idx_Y += incY;
    }
}

void mnblas_caxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY){
    int idx_Y = 0;
    int idx_X = 0;
    float alpha_val = *((float *) alpha);
    float * X_val = (float *) X;
    float * Y_val = (float *) Y; 
    while ((idx_Y < N) && (idx_X < N)){
        Y_val [idx_Y] = alpha_val * X_val [idx_X] + Y_val [idx_Y];
        Y_val [idx_Y+1] = alpha_val * X_val [idx_X+1] + Y_val [idx_Y+1];
        idx_X += 2 * incX;
        idx_Y += 2 * incY;
    }
}


void mnblas_zaxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY){
    int idx_Y = 0;
    int idx_X = 0;
    double alpha_val = *((double *) alpha);
    double * X_val = (double *) X;
    double * Y_val = (double *) Y;
    while ((idx_Y < N) && (idx_X < N)){
        Y_val [idx_Y] = alpha_val * X_val [idx_X] + Y_val [idx_Y];
        Y_val [idx_Y+1] = alpha_val * X_val [idx_X+1] + Y_val [idx_Y+1];
        idx_X += 2 * incX;
        idx_Y += 2 * incY;
    }
}
