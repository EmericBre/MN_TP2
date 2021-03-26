#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define NB_FOIS 10

float *vec1;

double *vec1d;

complexe_float_t *vec1cf;

complexe_double_t *vec1cd;

complexe_float_t *iamincf;
complexe_double_t *iamincd;

void init()
{
    vec1 = malloc(sizeof(float *) * VECSIZE);

    vec1d = malloc(sizeof(double *) * VECSIZE);

    vec1cf = malloc(sizeof(complexe_float_t *) * VECSIZE);

    vec1cd = malloc(sizeof(float *) * (VECSIZE * 2));

    iamincf = malloc(sizeof(complexe_float_t *));
    iamincd = malloc(sizeof(double *) * 2);
}

int main(int argc, char **argv)
{

    float *f1 = malloc(sizeof(float *) * 6);
    f1[0] = 3.0;
    f1[1] = 2.0;
    f1[2] = 4.0;
    f1[3] = -1.0;
    f1[4] = -5.0;
    f1[5] = 2.0;
    float f2;

    printf("\n");

    double *d1 = malloc(sizeof(double *) * 6);
    d1[0] = 3.0;
    d1[1] = 2.0;
    d1[2] = 4.0;
    d1[3] = -1.0;
    d1[4] = -5.0;
    d1[5] = 2.0;
    float d2;

    printf("\n");

    complexe_float_t *cf1 = malloc(sizeof(complexe_float_t *) * 6);
    cf1[0].real = 3.0;
    cf1[0].imaginary = -2.0;
    cf1[1].real = 2.0;
    cf1[1].imaginary = -6.0;
    cf1[2].real = 4.0;
    cf1[2].imaginary = 1.0;
    cf1[3].real = -1.0;
    cf1[3].imaginary = 3.0;
    cf1[4].real = -5.0;
    cf1[4].imaginary = 0.0;
    cf1[5].real = 2.0;
    cf1[5].imaginary = -2.0;
    float cf2;

    printf("\n");

    complexe_double_t *cd1 = malloc(sizeof(complexe_double_t *) * 6);
    cd1[0].real = 3.0;
    cd1[0].imaginary = -2.0;
    cd1[1].real = 2.0;
    cd1[1].imaginary = -6.0;
    cd1[2].real = 4.0;
    cd1[2].imaginary = 1.0;
    cd1[3].real = -1.0;
    cd1[3].imaginary = 3.0;
    cd1[4].real = -5.0;
    cd1[4].imaginary = 0.0;
    cd1[5].real = 2.0;
    cd1[5].imaginary = -2.0;
    float cd2;

    printf("\n");

    printf("FLOAT\n\nAvant iamin :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%f, ", f1[i]);
    }
    printf("\n");

    f2 = mnblas_isamin(6, f1, 1);

    printf("Après iamin :\nValeur absolue min : ");
    printf("%f\n", f2);

    printf("\n\n\nDOUBLE\n\nAvant iamin :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, ", d1[i]);
    }
    printf("\n");

    d2 = mnblas_idamin(6, d1, 1);

    printf("Après iamin :\nValeur absolue min : ");
    printf("%f\n", d2);

    printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant iamin :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cf1[i].real, cf1[i].imaginary);
    }
    printf("\n");

    cf2 = mnblas_icamin(6, cf1, 1);

    printf("\nAprès iamin :\nValeur absolue min : ");
    printf("%f\n", cf2);

    printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant iamin :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cd1[i].real, cd1[i].imaginary);
    }
    printf("\n");

    cd2 = mnblas_izamin(6, cd1, 1);

    printf("Après iamin :\nValeur absolue min : ");
    printf("%f\n", cd2);

    printf("\n\n=========================================================\n");
    printf("PERFORMANCES");
    printf("\n=========================================================\n\n\n");

    unsigned long long start, end;
    int i;

    init();

    iamincf[0].real = 0.0;
    iamincf[0].imaginary = 0.0;
    iamincd[0].real = 0.0;
    iamincd[0].imaginary = 0.0;

    printf("\n\n\nFLOAT\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_init(vec1, 1.0, VECSIZE);

        start = _rdtsc();
        f2 = mnblas_isamin(6, vec1, 1);
        end = _rdtsc();

        printf("mnblas_isamin %d : res = %3.2f nombre de cycles: %Ld \n", i, f2, end - start);
        calcul_flop("isamin ", VECSIZE, end - start);
    }

    printf("\n\n\nDOUBLE\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initd(vec1d, 1.0, VECSIZE);

        start = _rdtsc();
        d2 = mnblas_idamin(6, vec1d, 1);
        end = _rdtsc();

        printf("mnblas_idamin %d : res = %3.2f nombre de cycles: %Ld \n", i, d2, end - start);
        calcul_flop("idamin ", VECSIZE, end - start);
    }

    printf("\n\n\nCOMPLEXE FLOAT\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcf(vec1cf, iamincf[0], VECSIZE);

        start = _rdtsc();
        cf2 = mnblas_icamin(6, vec1cf, 1);
        end = _rdtsc();

        printf("mnblas_icamin %d : res = %3.2f nombre de cycles: %Ld \n", i, cf2, end - start);
        calcul_flop("icamin ", VECSIZE, end - start);
    }

    printf("\n\n\nCOMPLEXE DOUBLE\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcd(vec1cd, iamincd[0], VECSIZE);

        start = _rdtsc();
        cd2 = mnblas_izamin(6, vec1cd, 1);
        end = _rdtsc();

        printf("mnblas_izamin %d : res = %3.2f nombre de cycles: %Ld \n", i, cd2, end - start);
        calcul_flop("izamin ", VECSIZE, end - start);
    }

    printf("\n");

    exit(0);
}