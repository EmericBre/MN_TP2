#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"
#include "math.h"

#define max(a, b) ((a) > (b) ? (a) : (b))

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
                   const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const float alpha, const float *A, const int lda,
                   const float *X, const int incX, const float beta,
                   float *Y, const int incY)
{
   register unsigned int i;
   register unsigned int j;
   register float save;

   for (i = 0; i < M; i += incX) {
      save = beta*Y[i];
      for (j = 0; j < N; j += incY) {
         save += alpha*A[i*N+j]*X[j]; 
      }
      Y[i] = save;
   }

}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const double alpha, const double *A, const int lda,
                   const double *X, const int incX, const double beta,
                   double *Y, const int incY)
{
   register unsigned int i;
   register unsigned int j;
   register double save;

   for (i = 0; i < M; i += incX) {
      save = beta*Y[i];
      for (j = 0; j < N; j += incY) {
         save += alpha*A[i*N+j]*X[j]; 
      }
      Y[i] = save;
   }
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta,
                   void *Y, const int incY)
{
   register unsigned int i;
   register unsigned int j;
   register complexe_float_t save ;

   for (i = 0; i < M; i += incX) {
      save = mult_complexe_float(*(complexe_float_t*)beta, ((complexe_float_t*)Y)[i]);
      for (j = 0; j < N; j += incY) {
         save = add_complexe_float(save, mult_complexe_float(mult_complexe_float(*(complexe_float_t*)alpha, ((complexe_float_t*)A)[i*N+j]),((complexe_float_t*)X)[j]));
      }
      ((complexe_float_t*)Y)[i] = save;
   }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta,
                   void *Y, const int incY)
{
   register unsigned int i;
   register unsigned int j;
   register complexe_double_t save;

   for (i = 0; i < M; i += incX) {
      save = mult_complexe_double(*(complexe_double_t*)beta, ((complexe_double_t*)Y)[i]);
      for (j = 0; j < N; j += incY) {
         save = add_complexe_double(save, mult_complexe_double(mult_complexe_double(*(complexe_double_t*)alpha, ((complexe_double_t*)A)[i*N+j]),((complexe_double_t*)X)[j]));
      }
      ((complexe_double_t*)Y)[i] = save;
   }
}
