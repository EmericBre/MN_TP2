#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"
#include "math.h"

float mnblas_sasum(const int N, const float *X, const int incX)
{
  register unsigned int i;
  register float sum = 0.0;

  for (i = 0; i < N; i += incX)
  {
    sum += fabs(X[i]);
  }

  return sum;
}

double mnblas_dasum(const int N, const double *X, const int incX)
{
  register unsigned int i;
  register float sum = 0.0;

  for (i = 0; i < N; i += incX)
  {
    sum += fabs(X[i]);
  }

  return sum;
}

float mnblas_scasum(const int N, const void *X, const int incX)
{
  register unsigned int i;
  register float sum = 0.0;

  for (i = 0; i < N; i += incX)
  {
    sum += fabs(((complexe_float_t *)X)[i].real) + fabs(((complexe_float_t *)X)[i].imaginary);
  }

  return sum;
}

double mnblas_dzasum(const int N, const void *X, const int incX)
{
  register unsigned int i;
  register double sum = 0.0;

  for (i = 0; i < N; i += incX)
  {
    sum += fabs(((complexe_double_t *)X)[i].real) + fabs(((complexe_double_t *)X)[i].imaginary);
  }

  return sum;
}
