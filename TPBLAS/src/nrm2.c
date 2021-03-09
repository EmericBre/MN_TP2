#include <stdio.h>
#include <stdlib.h>
#include "mnblas.h"
#include "complexe.h"
#include "math.h"

float mnblas_snrm2(const int N, const float *X, const int incX)
{
  register unsigned int i;
  register float nrm1 = 0.0;
  register float nrm2 = 0.0;

  for (i = 0; i < N; i += incX)
  {
    nrm1 +=pow(X[i], 2);
  }
  nrm2 = sqrtf(nrm1);

  return nrm2;
}

double mnblas_dnrm2(const int N, const double *X, const int incX)
{
  register unsigned int i;
  double nrm1 = 0.0;
  double nrm2 = 0.0;

  for (i = 0; i < N; i += incX)
  {
    nrm1 += pow(X[i], 2);
  }
  nrm2 = sqrt(nrm1);

  return nrm2;
}

float mnblas_scnrm2(const int N, const void *X, const int incX)
{
  register unsigned int i;
  float nrm1 = 0.0;
  float nrm2 = 0.0;

  for (i = 0; i < N; i += incX)
  {
    nrm1 += pow(((complexe_float_t *)X)[i].real, 2) + pow(((complexe_float_t *)X)[i].imaginary, 2);
  }
  nrm2 = sqrt(nrm1);

  return nrm2;
}

double mnblas_dznrm2(const int N, const void *X, const int incX)
{
  register unsigned int i;
  double nrm1 = 0.0;
  double nrm2 = 0.0;

  for (i = 0; i < N; i += incX)
  {
    nrm1 += pow(((complexe_double_t *)X)[i].real, 2) + pow(((complexe_double_t *)X)[i].imaginary, 2);
  }
  nrm2 = sqrt(nrm1);

  return nrm2;
}
