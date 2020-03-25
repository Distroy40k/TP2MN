#include "mnblas.h"
#include "complexe2.h"

void mnblas_saxpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY){
    register int idx_Y = 0;
    register int idx_X = 0;
    while ((idx_Y < N) && (idx_X < N)){
        Y [idx_Y] = alpha * X [idx_X] + Y [idx_Y];
        idx_X += incX;
        idx_Y += incY;
    }
}

void mnblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY){
    register int idx_Y = 0;
    register int idx_X = 0;
    while ((idx_Y < N) && (idx_X < N)){
        Y [idx_Y] = alpha * X [idx_X] + Y [idx_Y];
        idx_X += incX;
        idx_Y += incY;
    }
}

void mnblas_caxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY){
    register int idx_Y = 0;
    register int idx_X = 0;
    complexe_float_t * alphap = (complexe_float_t *) alpha;
    complexe_float_t * X_val = (complexe_float_t *) X;
    complexe_float_t * Y_val = (complexe_float_t *) Y;
    while ((idx_Y < N) && (idx_X < N)){
        Y_val[idx_Y] = add_complexe_float(mult_complexe_float(X_val[idx_X],*alphap),Y_val[idx_Y]);
        idx_X += incX;
        idx_Y += incY;
    }
}


void mnblas_zaxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY){
  register int idx_Y = 0;
  register int idx_X = 0;
  complexe_double_t * alphap = (complexe_double_t *) alpha;
  complexe_double_t * X_val = (complexe_double_t *) X;
  complexe_double_t * Y_val = (complexe_double_t *) Y;
  while ((idx_Y < N) && (idx_X < N)){
      Y_val[idx_Y] = add_complexe_double(mult_complexe_double(X_val[idx_X],*alphap),Y_val[idx_Y]);
      idx_X += incX;
      idx_Y += incY;
  }
}
