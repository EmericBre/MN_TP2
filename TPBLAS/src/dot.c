#include "../include/mnblas.h"
#include "../include/complexe.h"
#include <stdio.h>


float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  float dot = 0.0 ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}

double mncblas_ddot(const int N, const double *X, const int incX, 
                 const double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  double dot = 0.0 ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}

void   mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  ((complexe_float_t *)dotu) [i].real = 0.0;
  ((complexe_float_t *)dotu) [i].imaginary = 0.0;

  
  for (i = 0 ; i < N ; i += incX)
    {
      ((complexe_float_t *)dotu) [i].real += ((complexe_float_t *)X) [i].real * ((complexe_float_t *)Y) [j].real ;
      ((complexe_float_t *)dotu) [i].real -= ((complexe_float_t *)X) [i].imaginary * ((complexe_float_t *)Y) [j].imaginary ;
      ((complexe_float_t *)dotu) [i].imaginary += ((complexe_float_t *)X) [i].real * ((complexe_float_t *)Y) [j].imaginary ;
      ((complexe_float_t *)dotu) [i].imaginary += ((complexe_float_t *)X) [i].imaginary * ((complexe_float_t *)Y) [j].real ;
      j+=incY ;
    }

  return;
}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  ((complexe_float_t *)dotc) [i].real = 0.0;
  ((complexe_float_t *)dotc) [i].imaginary = 0.0;

  
  for (i = 0 ; i < N ; i += incX)
    {
      ((complexe_float_t *)dotc) [i].real += ((complexe_float_t *)X) [i].real * ((complexe_float_t *)Y) [j].real ;
      ((complexe_float_t *)dotc) [i].real += ((complexe_float_t *)X) [i].imaginary * ((complexe_float_t *)Y) [j].imaginary ;
      ((complexe_float_t *)dotc) [i].imaginary += ((complexe_float_t *)X) [i].real * ((complexe_float_t *)Y) [j].imaginary ;
      ((complexe_float_t *)dotc) [i].imaginary -= ((complexe_float_t *)X) [i].imaginary * ((complexe_float_t *)Y) [j].real ;
      j+=incY ;
    }

  return;
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  ((complexe_double_t *)dotu) [i].real = 0.0;
  ((complexe_double_t *)dotu) [i].imaginary = 0.0;

  
  for (i = 0 ; i < N ; i += incX)
    {
      ((complexe_double_t *)dotu) [i].real += ((complexe_double_t *)X) [i].real * ((complexe_double_t *)Y) [j].real ;
      ((complexe_double_t *)dotu) [i].real -= ((complexe_double_t *)X) [i].imaginary * ((complexe_double_t *)Y) [j].imaginary ;
      ((complexe_double_t *)dotu) [i].imaginary += ((complexe_double_t *)X) [i].real * ((complexe_double_t *)Y) [j].imaginary ;
      ((complexe_double_t *)dotu) [i].imaginary += ((complexe_double_t *)X) [i].imaginary * ((complexe_double_t *)Y) [j].real ;
      j+=incY ;
    }

  return;
}
  
void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  ((complexe_double_t *)dotc) [i].real = 0.0;
  ((complexe_double_t *)dotc) [i].imaginary = 0.0;

  
  for (i = 0 ; i < N ; i += incX)
    {
      ((complexe_double_t *)dotc) [i].real += ((complexe_double_t *)X) [i].real * ((complexe_double_t *)Y) [j].real ;
      ((complexe_double_t *)dotc) [i].real += ((complexe_double_t *)X) [i].imaginary * ((complexe_double_t *)Y) [j].imaginary ;
      ((complexe_double_t *)dotc) [i].imaginary += ((complexe_double_t *)X) [i].real * ((complexe_double_t *)Y) [j].imaginary ;
      ((complexe_double_t *)dotc) [i].imaginary -= ((complexe_double_t *)X) [i].imaginary * ((complexe_double_t *)Y) [j].real ;
      j+=incY ;
    }

  return;
}




