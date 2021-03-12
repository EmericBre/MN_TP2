#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        4194304

int main (int argc, char **argv)
{

 float* fA= malloc (sizeof(float*)*12);
 float* fB= malloc (sizeof(float*)*12);
 float* fC= malloc (sizeof(float*)*9);

 double* dA= malloc (sizeof(double*)*12);
 double* dB= malloc (sizeof(double*)*12);
 double* dC= malloc (sizeof(double*)*12);

 complexe_float_t* cfA = malloc (sizeof(complexe_float_t*)*12);
 complexe_float_t* cfB = malloc (sizeof(complexe_float_t*)*12);
 complexe_float_t* cfC = malloc (sizeof(complexe_float_t*)*12);

 complexe_double_t* cdA = malloc (sizeof(double*)*24);
 complexe_double_t* cdB = malloc (sizeof(double*)*24);
 complexe_double_t* cdC = malloc (sizeof(double*)*18);

 int M = 3;
 int N = 3;
 int K = 4;
 complexe_float_t* alphacf = malloc (sizeof(complexe_float_t*));
 complexe_float_t* betacf = malloc (sizeof(complexe_float_t*));
 complexe_double_t* alphacd = malloc (sizeof(double*)*2);
 complexe_double_t* betacd = malloc (sizeof(double*)*2);

//  unsigned long long int start, end ;

 init_flop () ;

 float enter = 1.0;

 for (int i = 0; i < 12; i++) {
     fA[i] = enter;
     fB[i] = enter;
     enter++;
 }

 enter = 1.0;

 for (int i = 0; i < 9; i++) {
     fC[i] = enter;
     enter++;
 }

 float alphaf = 2.0;
 float betaf = 2.0;

 printf("\n\n\nFLOAT\n\nAvant gemm :\nMatrice A : ");
 for (int i = 0; i < 12; i++) {
     if (i%4 == 0) {
        printf("\n");
     }
     printf("%f ", fA[i]);
 }
 printf("\nMatrice B : ");
 for (int i = 0; i < 12; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%f ", fB[i]);
 }
 printf("\nMatrice C : ");
 for (int i = 0; i < 9; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%f ", fC[i]);
 }
 printf("\n");

 mncblas_sgemm(101, 111, 111, M, N, K, alphaf, fA, 0, fB, 0, betaf, fC, 0);

 printf("Après gemm :\nMatrice C : ");
 for (int i = 0; i < 9; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%f ", fC[i]);
 }
 printf("\n");
 




 double enterd = 1.0;

 for (int i = 0; i < 12; i++) {
     dA[i] = enterd;
     dB[i] = enterd;
     enterd++;
 }

 enterd = 1.0;

 for (int i = 0; i < 9; i++) {
     dC[i] = enterd;
     enterd++;
 }

 double alphad = 2.0;
 double betad = 2.0;

 printf("\n\n\nDOUBLE\n\nAvant gemm :\nMatrice A : ");
 for (int i = 0; i < 12; i++) {
     if (i%4 == 0) {
        printf("\n");
     }
     printf("%f ", dA[i]);
 }
 printf("\nMatrice B : ");
 for (int i = 0; i < 12; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%f ", dB[i]);
 }
 printf("\nMatrice C : ");
 for (int i = 0; i < 9; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%f ", dC[i]);
 }
 printf("\n");

 mncblas_dgemm(101, 111, 111, M, N, K, alphad, dA, 0, dB, 0, betad, dC, 0);

 printf("Après gemm :\nMatrice C : ");
 for (int i = 0; i < 9; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%f ", dC[i]);
 }
 printf("\n");





 complexe_float_t entercf = {1.0, 1.0};

 for (int i = 0; i < 12; i++) {
     cfA[i].real = entercf.real;
     cfA[i].imaginary = entercf.imaginary;
     cfB[i].real = entercf.real;
     cfB[i].imaginary = entercf.imaginary;
     entercf.real++;
     entercf.imaginary++;
 }

 entercf.real = 1.0; entercf.imaginary = 1.0;

 for (int i = 0; i < 9; i++) {
     cfC[i].real = entercf.real;
     cfC[i].imaginary = entercf.imaginary;
     entercf.real++;
     entercf.imaginary++;
 }

 alphacf[0].real = 2.0; alphacf[0].imaginary = 2.0;
 betacf[0].real = 2.0; betacf[0].imaginary = 2.0;

 printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant gemm :\nMatrice A : ");
 for (int i = 0; i < 12; i++) {
     if (i%4 == 0) {
        printf("\n");
     }
     printf("%lf %lf        ", cfA[i].real, cfA[i].imaginary);
 }
 printf("\nMatrice B : ");
 for (int i = 0; i < 12; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     
     printf("%lf %lf        ", cfB[i].real, cfB[i].imaginary);
 }
 printf("\nMatrice C : ");
 for (int i = 0; i < 9; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%lf %lf        ", cfC[i].real, cfC[i].imaginary);
 }
 printf("\n");

 mncblas_cgemm(101, 111, 111, M, N, K, alphacf, cfA, 0, cfB, 0, betacf, cfC, 0);

 printf("Après gemm :\nMatrice C : ");
 for (int i = 0; i < 9; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%lf %lf        ", cfC[i].real, cfC[i].imaginary);
 }
 printf("\n");





 complexe_double_t entercd = {1.0, 1.0};

 for (int i = 0; i < 12; i++) {
     cdA[i].real = entercd.real;
     cdA[i].imaginary = entercd.imaginary;
     cdB[i].real = entercd.real;
     cdB[i].imaginary = entercd.imaginary;
     entercd.real++;
     entercd.imaginary++;
 }

 entercd.real = 1.0; entercd.imaginary = 1.0;

 for (int i = 0; i < 9; i++) {
     cdC[i].real = entercd.real;
     cdC[i].imaginary = entercd.imaginary;
     entercd.real++;
     entercd.imaginary++;
 }

 alphacd[0].real = 2.0; alphacd[0].imaginary = 2.0;
 betacd[0].real = 2.0; betacd[0].imaginary = 2.0;

 printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant gemm :\nMatrice A : ");
 for (int i = 0; i < 12; i++) {
     if (i%4 == 0) {
        printf("\n");
     }
     printf("%lf %lf        ", cdA[i].real, cdA[i].imaginary);
 }
 printf("\nMatrice B : ");
 for (int i = 0; i < 12; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%lf %lf        ", cdB[i].real, cdB[i].imaginary);
 }
 printf("\nMatrice C : ");
 for (int i = 0; i < 9; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%lf %lf        ", cdC[i].real, cdC[i].imaginary);
 }
 printf("\n");

 mncblas_zgemm(101, 111, 111, M, N, K, alphacd, cdA, 0, cdB, 0, betacd, cdC, 0);

 printf("Après gemm :\nMatrice C : ");
 for (int i = 0; i < 9; i++) {
     if (i%3 == 0) {
        printf("\n");
     }
     printf("%lf %lf        ", cdC[i].real, cdC[i].imaginary);
 }
 printf("\n");

//  printf ("Addition de c1 et c2 : c3.r %f c3.i %f\n", c3.real, c3.imaginary) ;

//  start =_rdtsc () ;
 
//  for (i = 0 ; i < NB_FOIS; i++)
//    {
//      cd3 = add_complexe_double (cd1, cd2) ;
//    }

//  end = _rdtsc () ;

//   printf ("apres boucle cd3.real %f cd3.imaginary %f %lld cycles \n", cd3.real, cd3.imaginary, end-start) ;

//   calcul_flop ("calcul complexe ", NB_FOIS*4, end-start) ;
  exit (0) ;


}