#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "../include/complexe.h"


#include "flop.h"

#define    NB_FOIS        4194304

int main (int argc, char **argv)
{
 complexe_float_t c1= {1.0, 2.0} ;
 complexe_float_t c2= {3.0, 6.0} ;
 complexe_float_t c3 ;

 complexe_double_t cd1 = {10.0, 7.0}  ;
 complexe_double_t cd2 = {25.0, 32.0};
 complexe_double_t cd3 ;

 unsigned long long int start, end ;
 int i ;

 init_flop () ;

 printf("c1 = c1.r -> %f, c1.i -> %f\n", c1.real, c1.imaginary);
 printf("c2 = c2.r -> %f, c2.i -> %f\n", c2.real, c2.imaginary);

 printf("cd1 = cd1.r -> %f, cd1.i -> %f\n", cd1.real, cd1.imaginary);
 printf("cd2 = cd2.r -> %f, cd2.i -> %f\n", cd2.real, cd2.imaginary);
 
 c3 = add_complexe_float (c1, c2) ;

 printf ("Addition de c1 et c2 : c3.r %f c3.i %f\n", c3.real, c3.imaginary) ;


 cd3 = add_complexe_double (cd1, cd2) ;

 printf ("Addition double de cd1 et cd2 : cd3.r %f cd3.i %f\n", cd3.real, cd3.imaginary) ;

 c3 = mult_complexe_float(c1,c2);

 printf ("Multiplication de c1 et c2 : c3.r %f c3.i %f\n", c3.real, c3.imaginary) ;

 cd3 = mult_complexe_double (cd1, cd2) ;

 printf ("Multiplication double de cd1 et cd2 : cd3.r %f cd3.i %f\n", cd3.real, cd3.imaginary) ;

 c3 = div_complexe_float(c1,c2) ;

 printf ("Division de c1 et c2 : c3.r %f c3.i %f\n", c3.real, c3.imaginary) ;

 cd3 = div_complexe_double(cd1,cd2) ;

 printf ("Division double de cd1 et cd2 : cd3.r %f cd3.i %f\n", cd3.real, cd3.imaginary) ;

 start =_rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cd3 = add_complexe_double (cd1, cd2) ;
   }

 end = _rdtsc () ;

  printf ("apres boucle cd3.real %f cd3.imaginary %f %lld cycles \n", cd3.real, cd3.imaginary, end-start) ;

  calcul_flop ("calcul complexe ", NB_FOIS*4, end-start) ;
  exit (0) ;


}
