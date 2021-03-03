#include "mnblas.h"
#include "complexe.h"
#include <limits.h>
#include <math.h>
#include <stdlib.h>

CBLAS_INDEX mnblas_isamin(const int N, const float *X, const int incX)
{
    if ((N <= 0) || (incX <= 0))
    {
        return 0;
    }
    register unsigned int j;
    register float i = X[0];

    for (j = incX; j < N; j += incX)
    {
        if (fabs(X[j]) < i)
        {
            i = fabs(X[j]);
        }
    }

    return i;
}

CBLAS_INDEX mnblas_idamin(const int N, const double *X, const int incX)
{
    if ((N <= 0) || (incX <= 0))
    {
        return 0;
    }
    register unsigned int j;
    register double i = X[0];

    for (j = incX; j < N; j += incX)
    {
        if (fabs(X[j]) < i)
        {
            i = fabs(X[j]);
        }
    }

    return i;
}

CBLAS_INDEX mnblas_icamin(const int N, const void *X, const int incX)
{
    if ((N <= 0) || (incX <= 0))
    {
        return 0;
    }
    register unsigned int j;
    register float i = fabs(((complexe_float_t *)X)[0].real) + fabs(((complexe_float_t *)X)[0].imaginary);

    for (j = incX; j < N; j += incX)
    {
        if ((fabs(((complexe_float_t *)X)[j].real) + fabs(((complexe_float_t *)X)[j].imaginary)) < i)
        {
            i = fabs(((complexe_float_t *)X)[j].real) + fabs(((complexe_float_t *)X)[j].real);
        }
    }

    return i;
}

CBLAS_INDEX mnblas_izamin(const int N, const void *X, const int incX)
{
    if ((N <= 0) || (incX <= 0))
    {
        return 0;
    }
    register unsigned int j;
    register float i = fabs(((complexe_double_t *)X)[0].real) + fabs(((complexe_double_t *)X)[0].imaginary);

    for (j = incX; j < N; j += incX)
    {
        if (fabs(((complexe_double_t *)X)[j].real) + fabs(((complexe_double_t *)X)[j].imaginary) < i)
        {
            i =fabs(((complexe_double_t *)X)[j].real) + fabs(((complexe_double_t *)X)[j].imaginary);
        }
    }

    return i;
}
