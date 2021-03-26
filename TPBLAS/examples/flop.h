#define VECSIZE 65536

void init_flop () ;

void calcul_flop (char *message, int nb_operations_flottantes, unsigned long long int cycles) ;

void vector_init(float* V, float x, int size)
{
  register unsigned int i;

  for (i = 0; i < size; i++)
    V[i] = x;

  return;
}

void vector_initd(double* V, double x, int size)
{
  register unsigned int i;

  for (i = 0; i < size; i++)
    V[i] = x;

  return;
}

void vector_initcf(complexe_float_t* V, complexe_float_t x, int size)
{
  register unsigned int i;

  for (i = 0; i < size; i++)
    V[i] = x;

  return;
}

void vector_initcd(complexe_double_t* V, complexe_double_t x, int size)
{
  register unsigned int i;

  for (i = 0; i < size; i++)
    V[i] = x;

  return;
}

void vector_print(float* V)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("%f ", V[i]);
  printf("\n");

  return;
}

void vector_printd(double* V)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("%lf ", V[i]);
  printf("\n");

  return;
}

void vector_printcf(complexe_float_t* V)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("%f %f ", V[i].real, V[i].imaginary);
  printf("\n");

  return;
}

void vector_printcd(complexe_double_t* V)
{
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("%lf %lf ", V[i].real, V[i].imaginary);
  printf("\n");

  return;
}