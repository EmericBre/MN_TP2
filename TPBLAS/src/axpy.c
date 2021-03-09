#include "mnblas.h"
#include "complexe.h"

void mnblas_saxpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY)
{
    register unsigned int i;
    register unsigned int j;

    for (i = 0, j = 0; ((i < N) && (j < N)); i += incX, j += incY)
    {
        Y[j] = alpha * X[i]+Y[j];
    }

    return;
}

void mnblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY)
{
    register unsigned int i;
    register unsigned int j;

    for (i = 0, j = 0; ((i < N) && (j < N)); i += incX, j += incY)
    {
        Y[j] = alpha * X[i]+Y[j];
    }

    return;
}

void mnblas_caxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY)
{
    // register unsigned int i;
    // register unsigned int j;

    // for (i = 0, j = 0; ((i < N) && (j < N)); i += incX, j += incY)
    // {
    //     ((complexe_float_t *)Y)[j].real = ((complexe_float_t *)alpha)->real * ((complexe_float_t *)X)[i].real - ((complexe_float_t *)alpha)->imaginary * ((complexe_float_t *)X)[i].imaginary + ((complexe_float_t *)Y)[i].real;
    //     ((complexe_float_t *)Y)[j].imaginary = ((complexe_float_t *)alpha)->real * ((complexe_float_t *)X)[i].imaginary + ((complexe_float_t *)alpha)->imaginary * ((complexe_float_t *)X)[i].real +((complexe_float_t *)Y)[i].imaginary; 
    // }

    return;
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY)
{
    // register unsigned int i;
    // register unsigned int j;

    // for (i = 0, j = 0; ((i < N) && (j < N)); i += incX, j += incY)
    // {

    //     ((complexe_double_t *)Y)[j].real = 
    //     ((complexe_double_t *)Y)[j].imaginary = 
    // }

    return;
}