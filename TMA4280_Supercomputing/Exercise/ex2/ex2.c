#include <stdio.h>
#include <stdlib.h>

const double A[][3] = {{0.3, 0.4, 0.3},
                         {0.7, 0.1, 0.2},
                        {0.5, 0.5, 0.0}};

const double b[3] = {1.0, 1.0, 1.0};
const double a[3] = {0.1, 0.2, 0.3};

int main(int argc, char** argv)
{
  double y[3] = {0}; //sum vector y
  double Ab[3]={0}; //A multiplied by b
  double gb[3]={0}; //b multiplied by gamma
  int i=0, j=0;
  double gamma=0.0;
  double gammaVector[3]={gamma};
  
  if(argc<2)
  {
    fprintf(stdout, "Usage: ex2 <gamma>\n");
    return EXIT_FAILURE;
  }

  gamma = atoi(argv[1]);
  
/*---------------------------------------------------------------------------*/
/* Ab                                                                        */
/*---------------------------------------------------------------------------*/
  for(i=0; i<3; i++)
  {
    Ab[i] = multiply(A[i], b, 3);
  }

/*---------------------------------------------------------------------------*/
/* y = a + Ab                                                                */
/*---------------------------------------------------------------------------*/
  y = add(Ab, a, 3);

/*---------------------------------------------------------------------------*/
/* gb                                                                        */
/*---------------------------------------------------------------------------*/
  gb = multiply(gammaVector, b, 3);

/*---------------------------------------------------------------------------*/
/* x = a + gb                                                                */
/*---------------------------------------------------------------------------*/
  y = add(gb, a, 3);

/*---------------------------------------------------------------------------*/
/* alpha = xT*y                                                              */
/*---------------------------------------------------------------------------*/
  alpha = multiply(x, y, 3);


  fprintf(stdout, "result: y = %f %f %f\n", y[0], y[1], y[2]);
  return EXIT_SUCCESS;
}

double* add(double* vector1, double* vector2, int size)
{
  double* result = (double*)malloc(size*sizeof(double));
  int i=0;

  for(i=0; i<size; i++)
  {
    result[i] = vector1[i] + vector2[i];
  }

  return result;
}

double* multiply(double* vector1, double* vector2, int size)
{
  double* result = (double*)malloc(size*sizeof(double));
  int i=0;

  for(i=0; i<size; i++)
  {
    result[i] = vector1[i]*vector2[i];
  }
  return result;
}
