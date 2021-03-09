#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        4194304

int main (int argc, char **argv)
{

 float* f1= malloc (sizeof(float*)*6);

 double* d1= malloc (sizeof(double*)*6) ;

 complexe_float_t* cf1 = malloc (sizeof(complexe_float_t*)*6);

 complexe_double_t* cd1 = malloc (sizeof(complexe_double_t*)*6);

//  unsigned long long int start, end ;

 init_flop () ;

 float enter = 1.0;

 for (int i = 0; i < 6; i++) {
     f1[i] = enter;
     enter++;
 }

 printf("\n\n\nFLOAT\n\nAvant iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%f, ", *f1+i);
 }
 printf("\n");
 
 mnblas_isamin (6, f1, 1) ;

 printf("Après iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%f, ", *f1+i);
 }
 printf("\n");





 double enterd = 1.0;

 for (int i = 0; i < 6; i++) {
     d1[i] = enterd;
     enterd++;
 }

 printf("\n\n\nDOUBLE\n\nAvant iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%lf, ", *d1+i);
 }
 printf("\n");
 
 mnblas_idamin (6, d1, 1) ;

 printf("Après iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%lf, ", *d1+i);
 }
 printf("\n");





 complexe_float_t entercf = {1.0, 1.0};

 for (int i = 0; i < 6; i++) {
     cf1[i].real = entercf.real;
     cf1[i].imaginary = entercf.imaginary;
     entercf.real++;
     entercf.imaginary++;
 }

 printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%lf, %lf; ", cf1[i].real, cf1[i].imaginary);
 }
 printf("\n");
 
 mnblas_icamin (6, cf1, 1) ;

 printf("\nAprès iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%lf, %lf; ", cf1[i].real, cf1[i].imaginary);
 }
 printf("\n");





 complexe_double_t entercd = {1.0, 1.0};

 for (int i = 0; i < 6; i++) {
     cd1[i].real = entercd.real;
     cd1[i].imaginary = entercd.imaginary;
     entercd.real++;
     entercd.imaginary++;
 }

 printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%lf, %lf; ", cd1[i].real, cd1[i].imaginary);
 }
 printf("\n");
 
 mnblas_izamin (6, cd1, 1) ;

 printf("Après iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%lf, %lf; ", cd1[i].real, cd1[i].imaginary);
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