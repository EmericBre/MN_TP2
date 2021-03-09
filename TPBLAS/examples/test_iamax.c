#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        4194304

int main (int argc, char **argv)
{

 float* f1= malloc (sizeof(float*)*6);
 f1[0] = 3.0;
 f1[1] = 2.0;
 f1[2] = 4.0;
 f1[3] = -1.0;
 f1[4] = -5.0;
 f1[5] = 2.0;
 float f2;

 printf("\n");

 double* d1= malloc (sizeof(double*)*6);
 d1[0] = 3.0;
 d1[1] = 2.0;
 d1[2] = 4.0;
 d1[3] = -1.0;
 d1[4] = -5.0;
 d1[5] = 2.0;
 float d2;

 printf("\n");

 complexe_float_t* cf1 = malloc (sizeof(complexe_float_t*)*6);
 cf1[0].real = 3.0;  cf1[0].imaginary = -2.0;
 cf1[1].real = 2.0; cf1[1].imaginary = -6.0;
 cf1[2].real = 4.0; cf1[2].imaginary = 1.0;
 cf1[3].real = -1.0; cf1[3].imaginary = 3.0;
 cf1[4].real = -5.0; cf1[4].imaginary = 0.0;
 cf1[5].real = 2.0; cf1[5].imaginary = -2.0;
 float cf2;

 printf("\n");

 complexe_double_t* cd1 = malloc (sizeof(complexe_double_t*)*6);
 cd1[0].real = 3.0;  cd1[0].imaginary = -2.0;
 cd1[1].real = 2.0; cd1[1].imaginary = -6.0;
 cd1[2].real = 4.0; cd1[2].imaginary = 1.0;
 cd1[3].real = -1.0; cd1[3].imaginary = 3.0;
 cd1[4].real = -5.0; cd1[4].imaginary = 0.0;
 cd1[5].real = 2.0; cd1[5].imaginary = -2.0;
 float cd2;

 printf("\n");

//  unsigned long long int start, end ;

 init_flop () ;

 printf("FLOAT\n\nAvant iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%f, ", f1[i]);
 }
 printf("\n");
 
 f2 = mnblas_isamax (6, f1, 1) ;

 printf("Après iamax :\nValeur absolue max : ");
 printf("%f\n", f2);






 printf("\n\n\nDOUBLE\n\nAvant iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%lf, ", d1[i]);
 }
 printf("\n");
 
 d2 = mnblas_idamax (6, d1, 1) ;

 printf("Après iamax :\nValeur absolue max : ");
 printf("%f\n", d2);






 printf("\n\n\nCOMPLEXE_FLOAT\n\nAvant iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%lf, %lf; ", cf1[i].real, cf1[i].imaginary);
 }
 printf("\n");
 
 cf2 = mnblas_icamax (6, cf1, 1) ;

 printf("\nAprès iamax :\nValeur absolue max : ");
 printf("%f\n", cf2);






 printf("\n\n\nCOMPLEXE_DOUBLE\n\nAvant iamax :\nVecteur X : ");
 for (int i = 0; i < 6; i++) {
     printf("%lf, %lf; ", cd1[i].real, cd1[i].imaginary);
 }
 printf("\n");
 
 cd2 = mnblas_izamax (6, cd1, 1) ;

 printf("Après iamax :\nValeur absolue max : ");
 printf("%f\n", cd2);

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