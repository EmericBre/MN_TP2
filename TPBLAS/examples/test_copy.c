#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define NB_FOIS 10

float *vec1;
float *vec2;

double *vec1d;
double *vec2d;

complexe_float_t *vec1cf;
complexe_float_t *vec2cf;

complexe_double_t *vec1cd;
complexe_double_t *vec2cd;

complexe_float_t *copycf;
complexe_double_t *copycd;

void init()
{
    vec1 = malloc(sizeof(float *) * VECSIZE);
    vec2 = malloc(sizeof(float *) * VECSIZE);

    vec1d = malloc(sizeof(double *) * VECSIZE);
    vec2d = malloc(sizeof(double *) * VECSIZE);

    vec1cf = malloc(sizeof(complexe_float_t *) * VECSIZE);
    vec2cf = malloc(sizeof(complexe_float_t *) * VECSIZE);

    vec1cd = malloc(sizeof(float *) * (VECSIZE * 2));
    vec2cd = malloc(sizeof(float *) * (VECSIZE * 2));

    copycf = malloc(sizeof(complexe_float_t *));
    copycd = malloc(sizeof(double *) * 2);
}

int main(int argc, char **argv)
{

    float *f1 = malloc(sizeof(float *) * 6);
    float *f2 = malloc(sizeof(float *) * 6);

    double *d1 = malloc(sizeof(double *) * 6);
    double *d2 = malloc(sizeof(double *) * 6);

    complexe_float_t *cf1 = malloc(sizeof(complexe_float_t *) * 6);
    complexe_float_t *cf2 = malloc(sizeof(complexe_float_t *) * 6);

    complexe_double_t *cd1 = malloc(sizeof(double *) * 12);
    complexe_double_t *cd2 = malloc(sizeof(double *) * 12);

    float enter = 1.0;

    for (int i = 0; i < 6; i++)
    {
        f1[i] = enter;
        enter++;
    }

    for (int i = 0; i < 6; i++)
    {
        f2[i] = enter;
        enter++;
    }

    printf("\n\n\nFLOAT\n\nAvant copie :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%f, ", f1[i]);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%f, ", f2[i]);
    }
    printf("\n");

    mncblas_scopy(6, f1, 1, f2, 1);

    printf("Après copie :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%f, ", f1[i]);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%f, ", f2[i]);
    }
    printf("\n");

    double enterd = 1.0;

    for (int i = 0; i < 6; i++)
    {
        d1[i] = enterd;
        enterd++;
    }

    for (int i = 0; i < 6; i++)
    {
        d2[i] = enterd;
        enterd++;
    }

    printf("\n\n\nDOUBLE\n\nAvant copie :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, ", d1[i]);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, ", d2[i]);
    }
    printf("\n");

    mncblas_dcopy(6, d1, 1, d2, 1);

    printf("Après copie :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, ", d1[i]);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, ", d2[i]);
    }
    printf("\n");

    complexe_float_t entercf = {1.0, 1.0};

    for (int i = 0; i < 6; i++)
    {
        cf1[i].real = entercf.real;
        cf1[i].imaginary = entercf.imaginary;
        entercf.real++;
        entercf.imaginary++;
    }

    for (int i = 0; i < 6; i++)
    {
        cf2[i].real = entercf.real;
        cf2[i].imaginary = entercf.imaginary;
        entercf.real++;
        entercf.imaginary++;
    }

    printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant copie :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cf1[i].real, cf1[i].imaginary);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cf2[i].real, cf2[i].imaginary);
    }
    printf("\n");

    mncblas_ccopy(6, cf1, 1, cf2, 1);

    printf("\nAprès copie :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cf1[i].real, cf1[i].imaginary);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cf2[i].real, cf2[i].imaginary);
    }
    printf("\n");

    complexe_double_t entercd = {1.0, 1.0};

    for (int i = 0; i < 6; i++)
    {
        cd1[i].real = entercd.real;
        cd1[i].imaginary = entercd.imaginary;
        entercd.real++;
        entercd.imaginary++;
    }

    for (int i = 0; i < 6; i++)
    {
        cd2[i].real = entercd.real;
        cd2[i].imaginary = entercd.imaginary;
        entercd.real++;
        entercd.imaginary++;
    }

    printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant copie :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cd1[i].real, cd1[i].imaginary);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cd2[i].real, cd2[i].imaginary);
    }
    printf("\n");

    mncblas_zcopy(6, cd1, 1, cd2, 1);

    printf("Après copie :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cd1[i].real, cd1[i].imaginary);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cd2[i].real, cd2[i].imaginary);
    }
    printf("\n");



    printf("\n\n=========================================================\n");
    printf("PERFORMANCES");
    printf("\n=========================================================\n\n\n");

    unsigned long long start, end;
    int i;

    init();

    copycf[0].real = 0.0;
    copycf[0].imaginary = 0.0;
    copycd[0].real = 0.0;
    copycd[0].imaginary = 0.0;

    printf("\n\n\nFLOAT\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_init(vec1, 1.0, VECSIZE);
        vector_init(vec2, 2.0, VECSIZE);

        start = _rdtsc();
        mncblas_scopy(6, vec1, 1, vec2, 1);
        end = _rdtsc();

        printf("mncblas_scopy %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("scopy ", 2 * VECSIZE, end - start);
    }

    printf("\n\n\nDOUBLE\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initd(vec1d, 1.0, VECSIZE);
        vector_initd(vec2d, 2.0, VECSIZE);

        start = _rdtsc();
        mncblas_dcopy(6, vec1d, 1, vec2d, 1);
        end = _rdtsc();

        printf("mncblas_dcopy %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("dswap ", 2 * VECSIZE, end - start);
    }

    printf("\n\n\nCOMPLEXE FLOAT\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcf(vec1cf, copycf[0], VECSIZE);
        vector_initcf(vec2cf, copycf[0], VECSIZE);

        start = _rdtsc();
        mncblas_ccopy(6, vec1cf, 1, vec2cf, 1);
        end = _rdtsc();

        printf("mncblas_ccopy %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("ccopy ", 2 * VECSIZE, end - start);
    }

    printf("\n\n\nCOMPLEXE DOUBLE\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcd(vec1cd, copycd[0], VECSIZE);
        vector_initcd(vec2cd, copycd[0], VECSIZE);

        start = _rdtsc();
        mncblas_zcopy(6, vec1cd, 1, vec2cd, 1);
        end = _rdtsc();

        printf("mncblas_zcopy %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("zcopy ", 2 * VECSIZE, end - start);
    }

    printf("\n");

    exit(0);
}