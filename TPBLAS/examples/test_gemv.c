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

complexe_float_t *gemvcf;
complexe_double_t *gemvcd;

void init()
{
    vec1 = malloc(sizeof(float *) * 900);
    vec2 = malloc(sizeof(float *) * 30);
    vec3 = malloc(sizeof(float *) * 30);

    vec1d = malloc(sizeof(double *) * 900);
    vec2d = malloc(sizeof(double *) * 30);
    vec3d = malloc(sizeof(double *) * 30);

    vec1cf = malloc(sizeof(complexe_float_t *) * 900);
    vec2cf = malloc(sizeof(complexe_float_t *) * 30);
    vec3cf = malloc(sizeof(complexe_float_t *) * 30);

    vec1cd = malloc(sizeof(float *) * (900 * 2));
    vec2cd = malloc(sizeof(float *) * (30 * 2));
    vec3cd = malloc(sizeof(float *) * (30 * 2));

    gemvcf = malloc(sizeof(complexe_float_t *));
    gemvcd = malloc(sizeof(double *) * 2);
}

int main(int argc, char **argv)
{

    float *fA = malloc(sizeof(float *) * 12);
    float *fx = malloc(sizeof(float *) * 4);
    float *fy = malloc(sizeof(float *) * 3);

    double *dA = malloc(sizeof(double *) * 12);
    double *dx = malloc(sizeof(double *) * 4);
    double *dy = malloc(sizeof(double *) * 3);

    complexe_float_t *cfA = malloc(sizeof(complexe_float_t *) * 12);
    complexe_float_t *cfx = malloc(sizeof(complexe_float_t *) * 4);
    complexe_float_t *cfy = malloc(sizeof(complexe_float_t *) * 3);

    complexe_double_t *cdA = malloc(sizeof(double *) * 24);
    complexe_double_t *cdx = malloc(sizeof(double *) * 8);
    complexe_double_t *cdy = malloc(sizeof(double *) * 6);

    int M = 3;
    int N = 4;
    complexe_float_t *alphacf = malloc(sizeof(complexe_float_t *));
    complexe_float_t *betacf = malloc(sizeof(complexe_float_t *));
    complexe_double_t *alphacd = malloc(sizeof(double *) * 2);
    complexe_double_t *betacd = malloc(sizeof(double *) * 2);

    float enter = 1.0;

    for (int i = 0; i < 12; i++)
    {
        fA[i] = enter;
        enter++;
    }

    enter = 1.0;

    for (int i = 0; i < 4; i++)
    {
        fx[i] = enter;
        enter++;
    }

    enter = 1.0;

    for (int i = 0; i < 3; i++)
    {
        fy[i] = enter;
        enter++;
    }

    float alphaf = 2.0;
    float betaf = 2.0;

    printf("\n\n\nFLOAT\n\nAvant gemv :\nMatrice A : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 4 == 0)
        {
            printf("\n");
        }
        printf("%f ", fA[i]);
    }
    printf("\nVecteur X : ");
    for (int i = 0; i < 4; i++)
    {
        printf("%f ", fx[i]);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 3; i++)
    {
        printf("%f ", fy[i]);
    }
    printf("\n");

    mncblas_sgemv(101, 111, M, N, alphaf, fA, 0, fx, 1, betaf, fy, 1);

    printf("Après gemv :\nVecteur Y : ");
    for (int i = 0; i < 3; i++)
    {
        printf("%f ", fy[i]);
    }
    printf("\n");

    double enterd = 1.0;

    for (int i = 0; i < 12; i++)
    {
        dA[i] = enterd;
        enterd++;
    }

    enterd = 1.0;

    for (int i = 0; i < 4; i++)
    {
        dx[i] = enterd;
        enterd++;
    }

    enterd = 1.0;

    for (int i = 0; i < 3; i++)
    {
        dy[i] = enterd;
        enterd++;
    }

    float alphad = 2.0;
    float betad = 2.0;

    printf("\n\n\nDOUBLE\n\nAvant gemv :\nMatrice A : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 4 == 0)
        {
            printf("\n");
        }
        printf("%f ", dA[i]);
    }
    printf("\nVecteur X : ");
    for (int i = 0; i < 4; i++)
    {
        printf("%f ", dx[i]);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 3; i++)
    {
        printf("%f ", dy[i]);
    }
    printf("\n");

    mncblas_dgemv(101, 111, M, N, alphad, dA, 0, dx, 1, betad, dy, 1);

    printf("Après gemv :\nVecteur Y : ");
    for (int i = 0; i < 3; i++)
    {
        printf("%f ", dy[i]);
    }
    printf("\n");

    complexe_float_t entercf = {1.0, 1.0};

    for (int i = 0; i < 12; i++)
    {
        cfA[i].real = entercf.real;
        cfA[i].imaginary = entercf.imaginary;
        entercf.real++;
        entercf.imaginary++;
    }

    entercf.real = 1.0;
    entercf.imaginary = 1.0;

    for (int i = 0; i < 4; i++)
    {
        cfx[i].real = entercf.real;
        cfx[i].imaginary = entercf.imaginary;
        entercf.real++;
        entercf.imaginary++;
    }

    entercf.real = 1.0;
    entercf.imaginary = 1.0;

    for (int i = 0; i < 3; i++)
    {
        cfy[i].real = entercf.real;
        cfy[i].imaginary = entercf.imaginary;
        entercf.real++;
        entercf.imaginary++;
    }

    alphacf[0].real = 2.0;
    alphacf[0].imaginary = 2.0;
    betacf[0].real = 2.0;
    betacf[0].imaginary = 2.0;

    printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant gemv :\nMatrice A : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 4 == 0)
        {
            printf("\n");
        }
        printf("%lf %lf        ", cfA[i].real, cfA[i].imaginary);
    }
    printf("\nVecteur X : ");
    for (int i = 0; i < 4; i++)
    {
        printf("%lf %lf        ", cfx[i].real, cfx[i].imaginary);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 3; i++)
    {
        printf("%lf %lf        ", cfy[i].real, cfy[i].imaginary);
    }
    printf("\n");

    mncblas_cgemv(101, 111, M, N, alphacf, cfA, 0, cfx, 1, betacf, cfy, 1);

    printf("Après gemm :\nVecteur Y : ");
    for (int i = 0; i < 3; i++)
    {
        printf("%lf %lf        ", cfy[i].real, cfy[i].imaginary);
    }
    printf("\n");

    complexe_double_t entercd = {1.0, 1.0};

    for (int i = 0; i < 12; i++)
    {
        cdA[i].real = entercd.real;
        cdA[i].imaginary = entercd.imaginary;
        entercd.real++;
        entercd.imaginary++;
    }

    entercd.real = 1.0;
    entercd.imaginary = 1.0;

    for (int i = 0; i < 4; i++)
    {
        cdx[i].real = entercd.real;
        cdx[i].imaginary = entercd.imaginary;
        entercd.real++;
        entercd.imaginary++;
    }

    entercd.real = 1.0;
    entercd.imaginary = 1.0;

    for (int i = 0; i < 3; i++)
    {
        cdy[i].real = entercd.real;
        cdy[i].imaginary = entercd.imaginary;
        entercd.real++;
        entercd.imaginary++;
    }

    alphacd[0].real = 2.0;
    alphacd[0].imaginary = 2.0;
    betacd[0].real = 2.0;
    betacd[0].imaginary = 2.0;

    printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant gemv :\nMatrice A : ");
    for (int i = 0; i < 12; i++)
    {
        if (i % 4 == 0)
        {
            printf("\n");
        }
        printf("%lf %lf        ", cdA[i].real, cdA[i].imaginary);
    }
    printf("\nVecteur X : ");
    for (int i = 0; i < 4; i++)
    {
        printf("%lf %lf        ", cdx[i].real, cdx[i].imaginary);
    }
    printf("\nVecteur Y : ");
    for (int i = 0; i < 3; i++)
    {
        printf("%lf %lf        ", cdy[i].real, cdy[i].imaginary);
    }
    printf("\n");

    mncblas_zgemv(101, 111, M, N, alphacd, cdA, 0, cdx, 1, betacd, cdy, 1);

    printf("Après gemm :\nVecteur Y : ");
    for (int i = 0; i < 3; i++)
    {
        printf("%lf %lf        ", cdy[i].real, cdy[i].imaginary);
    }
    printf("\n");

    printf("\n\n=========================================================\n");
    printf("PERFORMANCES");
    printf("\n=========================================================\n\n\n");

    unsigned long long start, end;
    int i;

    init();

    gemvcf[0].real = 0.0;
    gemvcf[0].imaginary = 0.0;
    gemvcd[0].real = 0.0;
    gemvcd[0].imaginary = 0.0;

    printf("\n\n\nFLOAT\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_init(vec1, 1.0, 900);
        vector_init(vec2, 2.0, 30);
        vector_init(vec3, 3.0, 30);

        start = _rdtsc();
        mncblas_sgemv(101, 111, 30, 30, alphaf, vec1, 0, vec2, 1, betaf, vec3, 1);
        end = _rdtsc();

        printf("mncblas_sgemv %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("sgemv ", 960, end - start);
    }

    printf("\n\n\nDOUBLE\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initd(vec1d, 1.0, 900);
        vector_initd(vec2d, 2.0, 30);
        vector_initd(vec3d, 3.0, 30);

        start = _rdtsc();
        mncblas_dgemv(101, 111, 30, 30, alphad, vec1d, 0, vec2d, 1, betad, vec3d, 1);
        end = _rdtsc();

        printf("mncblas_dgemv %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("dgemv ", 960, end - start);
    }

    printf("\n\n\nCOMPLEXE FLOAT\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcf(vec1cf, gemvcf[0], 900);
        vector_initcf(vec2cf, gemvcf[0], 30);
        vector_initcf(vec3cf, gemvcf[0], 30);

        start = _rdtsc();
        mncblas_cgemv(101, 111, 30, 30, alphacf, vec1cf, 0, vec2cf, 1, betacf, vec3cf, 1);
        end = _rdtsc();

        printf("mncblas_cgemv %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("cgemv ", 960, end - start);
    }

    printf("\n\n\nCOMPLEXE DOUBLE\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcd(vec1cd, gemvcd[0], 900);
        vector_initcd(vec2cd, gemvcd[0], 30);
        vector_initcd(vec3cd, gemvcd[0], 30);

        start = _rdtsc();
        mncblas_zgemv(101, 111, 30, 30, alphacd, vec1cd, 0, vec2cd, 1, betacd, vec3cd, 1);
        end = _rdtsc();

        printf("mncblas_zgemv %d : nombre de cycles: %Ld \n", i, end - start);
        calcul_flop("zgemv ", 960, end - start);
    }

    printf("\n");

    free(fA);
    free(fx);
    free(fy);
    free(dA);
    free(dx);
    free(dy);
    free(cfA);
    free(cfx);
    free(cfy);
    free(cdA);
    free(cdx);
    free(cdy);
    free(alphacf);
    free(alphacd);
    free(betacf);
    free(betacd);
    free(vec1);
    free(vec2);
    free(vec3);
    free(vec1d);
    free(vec2d);
    free(vec3d);
    free(vec1cf);
    free(vec2cf);
    free(vec3cf);
    free(vec1cd);
    free(vec2cd);
    free(vec3cd);
    free(gemvcf);
    free(gemvcd);

    exit(0);
}