#include <stdio.h>
#include <stdlib.h>

#define VSIZE 3

//Function Prototypes
double* add(const double* vector1, const double* vector2, int size);
double* multiply(const double* vector1, const double* vector2, int size);

//const static declarations
const double A[][3] = {{0.3, 0.4, 0.3},
                         {0.7, 0.1, 0.2},
                        {0.5, 0.5, 0.0}};

const double b[3] = {1.0, 1.0, 1.0};
const double a[3] = {0.1, 0.2, 0.3};

//main
int main(int argc, char** argv)
{
  double *y = malloc(VSIZE*sizeof(double)); //sum vector y
  double *x = malloc(VSIZE*sizeof(double)); //sum vector x
  double *Ab = malloc(VSIZE*sizeof(double)); //A multiplied by b
  double *gb = malloc(VSIZE*sizeof(double)); //b multiplied by gamma
  double *gammaVector = malloc(VSIZE*sizeof(double));
  int i=0, j=0;
  double gamma=0.0;
  double alpha=0.0;
  
  if(argc<2)
  {
    fprintf(stdout, "Usage: ex2 <gamma>\n");
    return EXIT_FAILURE;
  }

  gamma = atoi(argv[1]);
  for(i=0; i<VSIZE; i++)
  {
    gammaVector[i]=gamma;
  }
  
/*---------------------------------------------------------------------------*/
/* Ab                                                                        */
/*---------------------------------------------------------------------------*/
  for(i=0; i<VSIZE; i++)
  {
    for(j=0; j<VSIZE; j++)
    {
      Ab[i]+=A[i][j]*x[j];
    }
  }
/*---------------------------------------------------------------------------*/
/* y = a + Ab                                                                */
/*---------------------------------------------------------------------------*/
  y = add(Ab, a, 3);

/*---------------------------------------------------------------------------*/
/* gb                                                                        */
/*---------------------------------------------------------------------------*/
  gb = multiply(gammaVector, b, VSIZE);

/*---------------------------------------------------------------------------*/
/* x = a + gb                                                                */
/*---------------------------------------------------------------------------*/
  y = add(gb, a, VSIZE);

/*---------------------------------------------------------------------------*/
/* alpha = xT*y                                                              */
/*---------------------------------------------------------------------------*/
  alpha = multiply(x, y, VSIZE);


  fprintf(stdout, "result: y = %f %f %f\n", y[0], y[1], y[2]);
  return EXIT_SUCCESS;
}

double* add(const double* vector1, const double* vector2, int size)
{
  double* result = (double*)malloc(size*sizeof(double));
  int i=0;

  for(i=0; i<size; i++)
  {
    result[i] = vector1[i] + vector2[i];
  }

  return result;
}

double* multiply(const double* vector1, const double* vector2, int size)
{
  double* result = (double*)malloc(size*sizeof(double));
  int i=0;

  for(i=0; i<size; i++)
  {
    result[i] = vector1[i]*vector2[i];
  }
  return result;
}
