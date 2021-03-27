#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define NB_FOIS 10

int main(int argc, char **argv)
{

    float *f1 = malloc(sizeof(float *) * VECSIZE);
    f1[0] = 3.0;
    f1[1] = 2.0;
    f1[2] = 4.0;
    f1[3] = -1.0;
    f1[4] = -5.0;
    f1[5] = 2.0;
    float f2;

    printf("\n");

    double *d1 = malloc(sizeof(double *) * VECSIZE);
    d1[0] = 3.0;
    d1[1] = 2.0;
    d1[2] = 4.0;
    d1[3] = -1.0;
    d1[4] = -5.0;
    d1[5] = 2.0;
    float d2;

    printf("\n");

    complexe_float_t *cf1 = malloc(sizeof(complexe_float_t *) * VECSIZE);
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

    complexe_double_t *cd1 = malloc(sizeof(complexe_double_t *) * VECSIZE);
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

    complexe_float_t* asumcf = malloc(sizeof(complexe_float_t *));
    complexe_double_t* asumcd = malloc(sizeof(double *) * 2);

    printf("\n");

    printf("FLOAT\n\nAvant asum :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%f, ", f1[i]);
    }
    printf("\n");

    f2 = mnblas_sasum(6, f1, 1);

    printf("Après asum :\nSomme : ");
    printf("%f\n", f2);

    printf("\n\n\nDOUBLE\n\nAvant asum :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, ", d1[i]);
    }
    printf("\n");

    d2 = mnblas_dasum(6, d1, 1);

    printf("Après asum :\nSomme : ");
    printf("%f\n", d2);

    printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant asum :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cf1[i].real, cf1[i].imaginary);
    }
    printf("\n");

    cf2 = mnblas_scasum(6, cf1, 1);

    printf("\nAprès asum :\nSomme : ");
    printf("%f\n", cf2);

    printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant asum :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cd1[i].real, cd1[i].imaginary);
    }
    printf("\n");

    cd2 = mnblas_dzasum(6, cd1, 1);

    printf("Après asum :\nSomme : ");
    printf("%f\n", cd2);


    printf("\n\n=========================================================\n");
    printf("PERFORMANCES");
    printf("\n=========================================================\n\n\n");

    unsigned long long start, end;
    int i;

    asumcf[0].real = 0.0;
    asumcf[0].imaginary = 0.0;
    asumcd[0].real = 0.0;
    asumcd[0].imaginary = 0.0;

    printf("\n\n\nFLOAT\n\n");

    float res;

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_init(f1, 1.0, VECSIZE);

        start = _rdtsc();
        res = mnblas_sasum(VECSIZE, f1, 1);
        end = _rdtsc();

        printf("mncblas_sasum %d : res = %3.2f nombre de cycles: %Ld \n", i, res, end - start);
        calcul_flop("sasum ", VECSIZE, end - start);
    }

    printf("\n\n\nDOUBLE\n\n");

    float resd;

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initd(d1, 1.0, VECSIZE);

        start = _rdtsc();
        resd = mnblas_dasum(VECSIZE, d1, 1);
        end = _rdtsc();

        printf("mncblas_dasum %d : res = %3.2f nombre de cycles: %Ld \n", i, resd, end - start);
        calcul_flop("dasum ", VECSIZE, end - start);
    }

    printf("\n\n\nCOMPLEXE FLOAT\n\n");

    float rescf;

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcf(cf1, asumcf[0], VECSIZE);

        start = _rdtsc();
        rescf = mnblas_scasum(VECSIZE, cf1, 1);
        end = _rdtsc();

        printf("mncblas_scasum %d : res = %3.2f nombre de cycles: %Ld \n", i, rescf, end - start);
        calcul_flop("scasum ", VECSIZE, end - start);
    }

    printf("\n\n\nCOMPLEXE DOUBLE\n\n");

    float rescd;

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcd(cd1, asumcd[0], VECSIZE);

        start = _rdtsc();
        rescd = mnblas_dzasum(VECSIZE, cd1, 1);
        end = _rdtsc();

        printf("mncblas_dzasum %d : res = %3.2f nombre de cycles: %Ld \n", i, rescd, end - start);
        calcul_flop("dzasum ", VECSIZE, end - start);
    }

    printf("\n");

    free(f1);
    free(d1);
    free(cd1);
    free(asumcf);
    free(asumcd);

    exit(0);
}