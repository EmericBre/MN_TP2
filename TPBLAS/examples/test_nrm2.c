#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define NB_FOIS 4194

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

    complexe_float_t* nrm2cf = malloc(sizeof(complexe_float_t *));
    complexe_double_t* nrm2cd = malloc(sizeof(double *) * 2);

    printf("\n");

    printf("FLOAT\n\nAvant nrm2 :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%f, ", f1[i]);
    }
    printf("\n");

    f2 = mnblas_snrm2(6, f1, 1);

    printf("Après nrm2 :\nRésultat : ");
    printf("%f\n", f2);

    printf("\n\n\nDOUBLE\n\nAvant nrm2 :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, ", d1[i]);
    }
    printf("\n");

    d2 = mnblas_dnrm2(6, d1, 1);

    printf("Après nrm2 :\nRésultat : ");
    printf("%f\n", d2);

    printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant nrm2 :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cf1[i].real, cf1[i].imaginary);
    }
    printf("\n");

    cf2 = mnblas_scnrm2(6, cf1, 1);

    printf("\nAprès nrm2 :\nRésultat : ");
    printf("%f\n", cf2);

    printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant nrm2 :\nVecteur X : ");
    for (int i = 0; i < 6; i++)
    {
        printf("%lf, %lf; ", cd1[i].real, cd1[i].imaginary);
    }
    printf("\n");

    cd2 = mnblas_dznrm2(6, cd1, 1);

    printf("Après nrm2 :\nRésultat : ");
    printf("%f\n", cd2);

    printf("\n\n=========================================================\n");
    printf("PERFORMANCES");
    printf("\n=========================================================\n\n\n");

    unsigned long long start, end;
    int i;
    long long int moyenne = 0;

    nrm2cf[0].real = 0.0;
    nrm2cf[0].imaginary = 0.0;
    nrm2cd[0].real = 0.0;
    nrm2cd[0].imaginary = 0.0;

    printf("\n\n\nFLOAT\n\n");

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_init(f1, 1.0, VECSIZE);

        start = _rdtsc();
        f2 = mnblas_snrm2(VECSIZE, f1, 1);
        end = _rdtsc();

        moyenne += end - start;
    }

    printf("mnblas_snrm2 moyenne : nombre de cycles: %Ld \n", moyenne/NB_FOIS);
    calcul_flop("snrm2 ", VECSIZE * NB_FOIS, moyenne);

    printf("\n\n\nDOUBLE\n\n");

    moyenne = 0;

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initd(d1, 1.0, VECSIZE);

        start = _rdtsc();
        d2 = mnblas_dnrm2(VECSIZE, d1, 1);
        end = _rdtsc();

        moyenne += end - start;
    }

    printf("mnblas_dnrm2 moyenne : nombre de cycles: %Ld \n", moyenne/NB_FOIS);
    calcul_flop("dnrm2 ", VECSIZE * NB_FOIS, moyenne);

    printf("\n\n\nCOMPLEXE FLOAT\n\n");

    moyenne = 0;

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcf(cf1, nrm2cf[0], VECSIZE);

        start = _rdtsc();
        cf2 = mnblas_scnrm2(VECSIZE, cf1, 1);
        end = _rdtsc();

        moyenne += end - start;
    }

    printf("mnblas_scnrm2 moyenne : nombre de cycles: %Ld \n", moyenne/NB_FOIS);
    calcul_flop("scnrm2 ", VECSIZE * NB_FOIS, moyenne);

    printf("\n\n\nCOMPLEXE DOUBLE\n\n");

    moyenne = 0;

    init_flop();

    for (i = 0; i < NB_FOIS; i++)
    {
        vector_initcd(cd1, nrm2cd[0], VECSIZE);

        start = _rdtsc();
        cd2 = mnblas_dznrm2(VECSIZE, cd1, 1);
        end = _rdtsc();

        moyenne += end - start;
    }

    printf("mnblas_dznrm2 moyenne : nombre de cycles: %Ld \n", moyenne/NB_FOIS);
    calcul_flop("dznrm2 ", VECSIZE * NB_FOIS, moyenne);

    printf("\n");

    free(f1);
    free(d1);
    free(cd1);
    free(nrm2cf);
    free(nrm2cd);

    exit(0);
}