#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        4194304

int main (int argc, char **argv)
{

 float* fA= malloc (sizeof(float*)*12);
 float* fx= malloc (sizeof(float*)*4);
 float* fy= malloc (sizeof(float*)*3);

 double* dA= malloc (sizeof(double*)*12);
 double* dx= malloc (sizeof(double*)*4);
 double* dy= malloc (sizeof(double*)*3);

 complexe_float_t* cfA = malloc (sizeof(complexe_float_t*)*12);
 complexe_float_t* cfx = malloc (sizeof(complexe_float_t*)*4);
 complexe_float_t* cfy = malloc (sizeof(complexe_float_t*)*3);

 complexe_double_t* cdA = malloc (sizeof(double*)*24);
 complexe_double_t* cdx = malloc (sizeof(double*)*8);
 complexe_double_t* cdy = malloc (sizeof(double*)*6);

 int M = 3;
 int N = 4;
 complexe_float_t* alphacf = malloc (sizeof(complexe_float_t*));
 complexe_float_t* betacf = malloc (sizeof(complexe_float_t*));
 complexe_double_t* alphacd = malloc (sizeof(double*)*2);
 complexe_double_t* betacd = malloc (sizeof(double*)*2);

//  unsigned long long int start, end ;

 init_flop () ;

 float enter = 1.0;

 for (int i = 0; i < 12; i++) {
     fA[i] = enter;
     enter++;
 }

 enter = 1.0;

 for (int i = 0; i < 4; i++) {
     fx[i] = enter;
     enter++;
 }

 enter = 1.0;

 for (int i = 0; i < 3; i++) {
     fy[i] = enter;
     enter++;
 }

 float alphaf = 2.0;
 float betaf = 2.0;

 printf("\n\n\nFLOAT\n\nAvant gemv :\nMatrice A : ");
 for (int i = 0; i < 12; i++) {
     if (i%4 == 0) {
        printf("\n");
     }
     printf("%f ", fA[i]);
 }
 printf("\nVecteur X : ");
 for (int i = 0; i < 4; i++) {
     printf("%f ", fx[i]);
 }
 printf("\nVecteur Y : ");
 for (int i = 0; i < 3; i++) {
     printf("%f ", fy[i]);
 }
 printf("\n");

 mncblas_sgemv(101, 111, M, N, alphaf, fA, 0, fx, 1, betaf, fy, 1);

 printf("Après gemv :\nVecteur Y : ");
 for (int i = 0; i < 3; i++) {
     printf("%f ", fy[i]);
 }
 printf("\n");
 




 double enterd = 1.0;

 for (int i = 0; i < 12; i++) {
     dA[i] = enterd;
     enterd++;
 }

 enterd = 1.0;

 for (int i = 0; i < 4; i++) {
     dx[i] = enterd;
     enterd++;
 }

 enterd = 1.0;

 for (int i = 0; i < 3; i++) {
     dy[i] = enterd;
     enterd++;
 }

 float alphad = 2.0;
 float betad = 2.0;

 printf("\n\n\nDOUBLE\n\nAvant gemv :\nMatrice A : ");
 for (int i = 0; i < 12; i++) {
     if (i%4 == 0) {
        printf("\n");
     }
     printf("%f ", dA[i]);
 }
 printf("\nVecteur X : ");
 for (int i = 0; i < 4; i++) {
     printf("%f ", dx[i]);
 }
 printf("\nVecteur Y : ");
 for (int i = 0; i < 3; i++) {
     printf("%f ", dy[i]);
 }
 printf("\n");

 mncblas_dgemv(101, 111, M, N, alphad, dA, 0, dx, 1, betad, dy, 1);

 printf("Après gemv :\nVecteur Y : ");
 for (int i = 0; i < 3; i++) {
     printf("%f ", dy[i]);
 }
 printf("\n");





 complexe_float_t entercf = {1.0, 1.0};

 for (int i = 0; i < 12; i++) {
     cfA[i].real = entercf.real;
     cfA[i].imaginary = entercf.imaginary;
     entercf.real++;
     entercf.imaginary++;
 }

 entercf.real = 1.0; entercf.imaginary = 1.0;

 for (int i =0; i < 4; i++) {
     cfx[i].real = entercf.real;
     cfx[i].imaginary = entercf.imaginary;
     entercf.real++;
     entercf.imaginary++;
 }

 entercf.real = 1.0; entercf.imaginary = 1.0;

 for (int i = 0; i < 3; i++) {
     cfy[i].real = entercf.real;
     cfy[i].imaginary = entercf.imaginary;
     entercf.real++;
     entercf.imaginary++;
 }

 alphacf[0].real = 2.0; alphacf[0].imaginary = 2.0;
 betacf[0].real = 2.0; betacf[0].imaginary = 2.0;

 printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant gemv :\nMatrice A : ");
 for (int i = 0; i < 12; i++) {
     if (i%4 == 0) {
        printf("\n");
     }
     printf("%lf %lf        ", cfA[i].real, cfA[i].imaginary);
 }
 printf("\nVecteur X : ");
 for (int i = 0; i < 4; i++) {
     printf("%lf %lf        ", cfx[i].real, cfx[i].imaginary);
 }
 printf("\nVecteur Y : ");
 for (int i = 0; i < 3; i++) {
     printf("%lf %lf        ", cfy[i].real, cfy[i].imaginary);
 }
 printf("\n");

 mncblas_cgemv(101, 111, M, N, alphacf, cfA, 0, cfx, 1, betacf, cfy, 1);

 printf("Après gemm :\nVecteur Y : ");
 for (int i = 0; i < 3; i++) {
     printf("%lf %lf        ", cfy[i].real, cfy[i].imaginary);
 }
 printf("\n");





 complexe_double_t entercd = {1.0, 1.0};

 for (int i = 0; i < 12; i++) {
     cdA[i].real = entercd.real;
     cdA[i].imaginary = entercd.imaginary;
     entercd.real++;
     entercd.imaginary++;
 }

 entercd.real = 1.0; entercd.imaginary = 1.0;

 for (int i =0; i < 4; i++) {
     cdx[i].real = entercd.real;
     cdx[i].imaginary = entercd.imaginary;
     entercd.real++;
     entercd.imaginary++;
 }

 entercd.real = 1.0; entercd.imaginary = 1.0;

 for (int i = 0; i < 3; i++) {
     cdy[i].real = entercd.real;
     cdy[i].imaginary = entercd.imaginary;
     entercd.real++;
     entercd.imaginary++;
 }

 alphacd[0].real = 2.0; alphacd[0].imaginary = 2.0;
 betacd[0].real = 2.0; betacd[0].imaginary = 2.0;

 printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant gemv :\nMatrice A : ");
 for (int i = 0; i < 12; i++) {
     if (i%4 == 0) {
        printf("\n");
     }
     printf("%lf %lf        ", cdA[i].real, cdA[i].imaginary);
 }
 printf("\nVecteur X : ");
 for (int i = 0; i < 4; i++) {
     printf("%lf %lf        ", cdx[i].real, cdx[i].imaginary);
 }
 printf("\nVecteur Y : ");
 for (int i = 0; i < 3; i++) {
     printf("%lf %lf        ", cdy[i].real, cdy[i].imaginary);
 }
 printf("\n");

 mncblas_zgemv(101, 111, M, N, alphacd, cdA, 0, cdx, 1, betacd, cdy, 1);

 printf("Après gemm :\nVecteur Y : ");
 for (int i = 0; i < 3; i++) {
     printf("%lf %lf        ", cdy[i].real, cdy[i].imaginary);
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