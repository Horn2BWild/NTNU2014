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
  double *y = NULL;  //sum vector y
  double *x = NULL;  //sum vector x
  double *Ab = NULL; //A multiplied by b
  double *gb = NULL; //b multiplied by gamma
  double *gammaVector = NULL; 
  int i=0, j=0;
  double gamma=0.0;
 // double alpha=0.0;
  double* alpha = NULL; 
  
  if(argc<2)
  {
    fprintf(stdout, "Usage: ex2 <gamma>\n");
    return EXIT_FAILURE;
  }

  *gammaVector = (double*)malloc(VSIZE*sizeof(double));
  gamma = atoi(argv[1]);
  for(i=0; i<VSIZE; i++)
  {
    gammaVector[i]=gamma;
  }
  
/*---------------------------------------------------------------------------*/
/* Ab                                                                        */
/*---------------------------------------------------------------------------*/
  *Ab = (double*)malloc(VSIZE*sizeof(double));

  for(i=0; i<VSIZE; i++)
  {
    for(j=0; j<VSIZE; j++)
    {
      Ab[i]+=A[i][j]*b[j];
    }
  }
/*---------------------------------------------------------------------------*/
/* y = a + Ab                                                                */
/*---------------------------------------------------------------------------*/
  *y = (double*)malloc(VSIZE*sizeof(double));
  y = add(Ab, a, 3);

/*---------------------------------------------------------------------------*/
/* gb                                                                        */
/*---------------------------------------------------------------------------*/
  *gb = (double*)malloc(VSIZE*sizeof(double));
  gb = multiply(gammaVector, b, VSIZE);

/*---------------------------------------------------------------------------*/
/* x = a + gb                                                                */
/*---------------------------------------------------------------------------*/
  *x = (double*)malloc(VSIZE*sizeof(double));
  x = add(gb, a, VSIZE);

/*---------------------------------------------------------------------------*/
/* alpha = xT*y                                                              */
/*---------------------------------------------------------------------------*/
  *alpha = (double*)malloc(sizeof(double));
  alpha = multiply(x, y, VSIZE);


  fprintf(stdout, "result: y = %f %f %f\n", y[0], y[1], y[2]);

  free(y);
  free(x);
  free(Ab);
  free(gb);
  free(gammaVector);
  free(alpha);

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
