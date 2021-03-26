#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define NB_FOIS 10

float *vec1;
float *vec2;
float *vec3;

double *vec1d;
double *vec2d;
double *vec3d;

complexe_float_t *vec1cf;
complexe_float_t *vec2cf;
complexe_float_t *vec3cf;

complexe_double_t *vec1cd;
complexe_double_t *vec2cd;
complexe_double_t *vec3cd;

complexe_float_t *gemmcf;
complexe_double_t *gemmcd;

void init()
{
    vec1 = malloc(sizeof(float *) * 900);
    vec2 = malloc(sizeof(float *) * 900);
    vec3 = malloc(sizeof(float *) * 900);

    vec1d = malloc(sizeof(double *) * 900);
    vec2d = malloc(sizeof(double *) * 900);
    vec3d = malloc(sizeof(double *) * 900);

    vec1cf = malloc(sizeof(complexe_float_t *) * 900);
    vec2cf = malloc(sizeof(complexe_float_t *) * 900);
    vec3cf = malloc(sizeof(complexe_float_t *) * 900);

    vec1cd = malloc(sizeof(float *) * (900 * 2));
    vec2cd = malloc(sizeof(float *) * (900 * 2));
    vec3cd = malloc(sizeof(float *) * (900 * 2));

    gemmcf = malloc(sizeof(complexe_float_t *));
    gemmcd = malloc(sizeof(double *) * 2);
}

int main(int argc, char **argv)
{

    float *fA = malloc(sizeof(float *) * 12);
    float *fB = malloc(sizeof(float *) * 12);
    float *fC = malloc(sizeof(float *) * 9);

    double *dA = malloc(sizeof(double *) * 12);
    double *dB = malloc(sizeof(double *) * 12);
    double *dC = malloc(sizeof(double *) * 12);

    complexe_float_t *cfA = malloc(sizeof(complexe_float_t *) * 12);
    complexe_float_t *cfB = malloc(sizeof(complexe_float_t *) * 12);
    complexe_float_t *cfC = malloc(sizeof(complexe_float_t *) * 12);

    complexe_double_t *cdA = malloc(sizeof(double *) * 24);
    complexe_double_t *cdB = malloc(sizeof(double *) * 24);
    complexe_double_t *cdC = malloc(sizeof(double *) * 18);

    int M = 3;
    int N = 3;
    int K = 4;
    complexe_float_t *alphacf = malloc(sizeof(complexe_float_t *));
    complexe_float_t *betacf = malloc(sizeof(complexe_float_t *));
    complexe_double_t *alphacd = malloc(sizeof(double *) * 2);
    complexe_double_t *betacd = malloc(sizeof(double *) * 2);

    float enter = 1.0;

    for (int i = 0; i < 12; i++)
    {
        fA[i] = enter;
        fB[i] = enter;
        enter++;
    }

    enter = 1.0;

    for (int i = 0; i < 9; i++)
    {
        fC[i] = enter;
        enter++;
    }

    float alphaf = 2.0;
    float betaf = 2.0;

    printf("\n\n\nFLOAT\n\nAvant gemm :\nMatrice A : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 4 == 0)
        {
            printf("\n");
        }
        printf("%f ", fA[i]);
    }
    printf("\nMatrice B : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%f ", fB[i]);
    }
    printf("\nMatrice C : ");
    for (int i = 0; i < 9; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%f ", fC[i]);
    }
    printf("\n");

    mncblas_sgemm(101, 111, 111, M, N, K, alphaf, fA, 0, fB, 0, betaf, fC, 0);

    printf("Après gemm :\nMatrice C : ");
    for (int i = 0; i < 9; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%f ", fC[i]);
    }
    printf("\n");

    double enterd = 1.0;

    for (int i = 0; i < 12; i++)
    {
        dA[i] = enterd;
        dB[i] = enterd;
        enterd++;
    }

    enterd = 1.0;

    for (int i = 0; i < 9; i++)
    {
        dC[i] = enterd;
        enterd++;
    }

    double alphad = 2.0;
    double betad = 2.0;

    printf("\n\n\nDOUBLE\n\nAvant gemm :\nMatrice A : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 4 == 0)
        {
            printf("\n");
        }
        printf("%f ", dA[i]);
    }
    printf("\nMatrice B : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%f ", dB[i]);
    }
    printf("\nMatrice C : ");
    for (int i = 0; i < 9; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%f ", dC[i]);
    }
    printf("\n");

    mncblas_dgemm(101, 111, 111, M, N, K, alphad, dA, 0, dB, 0, betad, dC, 0);

    printf("Après gemm :\nMatrice C : ");
    for (int i = 0; i < 9; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%f ", dC[i]);
    }
    printf("\n");

    complexe_float_t entercf = {1.0, 1.0};

    for (int i = 0; i < 12; i++)
    {
        cfA[i].real = entercf.real;
        cfA[i].imaginary = entercf.imaginary;
        cfB[i].real = entercf.real;
        cfB[i].imaginary = entercf.imaginary;
        entercf.real++;
        entercf.imaginary++;
    }

    entercf.real = 1.0;
    entercf.imaginary = 1.0;

    for (int i = 0; i < 9; i++)
    {
        cfC[i].real = entercf.real;
        cfC[i].imaginary = entercf.imaginary;
        entercf.real++;
        entercf.imaginary++;
    }

    alphacf[0].real = 2.0;
    alphacf[0].imaginary = 2.0;
    betacf[0].real = 2.0;
    betacf[0].imaginary = 2.0;

    printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant gemm :\nMatrice A : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 4 == 0)
        {
            printf("\n");
        }
        printf("%lf %lf        ", cfA[i].real, cfA[i].imaginary);
    }
    printf("\nMatrice B : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }

        printf("%lf %lf        ", cfB[i].real, cfB[i].imaginary);
    }
    printf("\nMatrice C : ");
    for (int i = 0; i < 9; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%lf %lf        ", cfC[i].real, cfC[i].imaginary);
    }
    printf("\n");

    mncblas_cgemm(101, 111, 111, M, N, K, alphacf, cfA, 0, cfB, 0, betacf, cfC, 0);

    printf("Après gemm :\nMatrice C : ");
    for (int i = 0; i < 9; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%lf %lf        ", cfC[i].real, cfC[i].imaginary);
    }
    printf("\n");

    complexe_double_t entercd = {1.0, 1.0};

    for (int i = 0; i < 12; i++)
    {
        cdA[i].real = entercd.real;
        cdA[i].imaginary = entercd.imaginary;
        cdB[i].real = entercd.real;
        cdB[i].imaginary = entercd.imaginary;
        entercd.real++;
        entercd.imaginary++;
    }

    entercd.real = 1.0;
    entercd.imaginary = 1.0;

    for (int i = 0; i < 9; i++)
    {
        cdC[i].real = entercd.real;
        cdC[i].imaginary = entercd.imaginary;
        entercd.real++;
        entercd.imaginary++;
    }

    alphacd[0].real = 2.0;
    alphacd[0].imaginary = 2.0;
    betacd[0].real = 2.0;
    betacd[0].imaginary = 2.0;

    printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant gemm :\nMatrice A : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 4 == 0)
        {
            printf("\n");
        }
        printf("%lf %lf        ", cdA[i].real, cdA[i].imaginary);
    }
    printf("\nMatrice B : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%lf %lf        ", cdB[i].real, cdB[i].imaginary);
    }
    printf("\nMatrice C : ");
    for (int i = 0; i < 9; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%lf %lf        ", cdC[i].real, cdC[i].imaginary);
    }
    printf("\n");

    mncblas_zgemm(101, 111, 111, M, N, K, alphacd, cdA, 0, cdB, 0, betacd, cdC, 0);

    printf("Après gemm :\nMatrice C : ");
    for (int i = 0; i < 9; i++)
    {
        if (i % 3 == 0)
        {
            printf("\n");
        }
        printf("%lf %lf        ", cdC[i].real, cdC[i].imaginary);
    }
    printf("\n");

    printf("\n\n=========================================================\n");
    printf("PERFORMANCES");
    printf("\n=========================================================\n\n\n");

    unsigned long long start, end;
    int i;

    init();

    gemmcf[0].real = 0.0;
    gemmcf[0].imaginary = 0.0;
    gemmcd[0].real = 0.0;
    gemmcd[0].imaginary = 0.0;

    printf("\n\n\nFLOAT\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_init(vec1, 1.0, 900);
        vector_init(vec2, 2.0, 900);
        vector_init(vec3, 3.0, 900);

        start = _rdtsc();
        mncblas_sgemm(101, 111, 111, 30, 30, 30, alphaf, vec1, 0, vec2, 0, betaf, vec3, 0);
        end = _rdtsc();

        printf("mncblas_sgemm %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("sgemm ", 3 * 900, end - start);
    }

    printf("\n\n\nDOUBLE\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initd(vec1d, 1.0, 900);
        vector_initd(vec2d, 2.0, 900);
        vector_initd(vec3d, 3.0, 900);

        start = _rdtsc();
        mncblas_dgemm(101, 111, 111, 30, 30, 30, alphad, vec1d, 0, vec2d, 0, betad, vec3d, 0);
        end = _rdtsc();

        printf("mncblas_dgemm %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("dgemm ", 3 * 900, end - start);
    }

    printf("\n\n\nCOMPLEXE FLOAT\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcf(vec1cf, gemmcf[0], 900);
        vector_initcf(vec2cf, gemmcf[0], 900);
        vector_initcf(vec3cf, gemmcf[0], 900);

        start = _rdtsc();
        mncblas_cgemm(101, 111, 111, 30, 30, 30, alphacf, vec1cf, 0, vec2cf, 0, betacf, vec3cf, 0);
        end = _rdtsc();

        printf("mncblas_cgemm %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("cgemm ", 3 * 900, end - start);
    }

    printf("\n\n\nCOMPLEXE DOUBLE\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcd(vec1cd, gemmcd[0], 900);
        vector_initcd(vec2cd, gemmcd[0], 900);
        vector_initcd(vec3cd, gemmcd[0], 900);

        start = _rdtsc();
        mncblas_zgemm(101, 111, 111, 30, 30, 30, alphacd, vec1cd, 0, vec2cd, 0, betacd, vec3cd, 0);
        end = _rdtsc();

        printf("mncblas_zgemm %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("zgemm ", 3 * 900, end - start);
    }

    printf("\n");

    exit(0);
}