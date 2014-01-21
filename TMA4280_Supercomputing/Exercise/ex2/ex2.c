#include <stdio.h>
#include <stdlib.h>

// defining the vector size
#define VSIZE 3

// Function Prototypes
double* add(const double* vector1, const double* vector2, int size);
double* multiply(const double* vector1, const double* vector2, int size);

// const static declarations
const double A[][3] = {{0.3, 0.4, 0.3},
                         {0.7, 0.1, 0.2},
                        {0.5, 0.5, 0.0}};

const double b[3] = {1.0, 1.0, 1.0};
const double a[3] = {0.1, 0.2, 0.3};

// main
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
  
// commandline argument
  if(argc<2)
  {
    fprintf(stdout, "Usage: ex2 <gamma>\n");
    return EXIT_FAILURE;
  }

// creating gammaVector from commandline input
  gammaVector = (double*)malloc(VSIZE*sizeof(double));
  gamma = atof(argv[1]);
  for(i=0; i<VSIZE; i++)
  {
    gammaVector[i]=gamma;
/*---DECOMMENT FOR DEBUGGING PURPOSES---
    printf("gamma: %f", gamma);
*/
  }
  
// calculation of A*b
  Ab = (double*)malloc(VSIZE*sizeof(double));
  for(i=0; i<VSIZE; i++)
  {
    for(j=0; j<VSIZE; j++)
    {
      Ab[i]+=A[i][j]*b[j];
    }
  }

//  y = a + Ab   
  y = (double*)malloc(VSIZE*sizeof(double));
  y = add(Ab, a, 3);

// gamma*b     
  gb = (double*)malloc(VSIZE*sizeof(double));
  gb = multiply(gammaVector, b, VSIZE);

/*---DECOMMENT FOR DEBUGGING PURPOSES---
  for(i=0; i<VSIZE; i++)
    printf("gb[%d] = %f", i, gb[i]);
*/

// x = a + gb             
  x = (double*)malloc(VSIZE*sizeof(double));
  x = add(gb, a, VSIZE);

// alpha = xT*y        
  alpha = (double*)malloc(sizeof(double));
  alpha = multiply(x, y, VSIZE);

// Outputs           
  fprintf(stdout, "-----------------Outputs-----------------\n");
  fprintf(stdout, "  x = [%f %f %f]\n", x[0], x[1], x[2]);
  fprintf(stdout, "  y = [%f %f %f]\n", y[0], y[1], y[2]);
  fprintf(stdout, "  alpha = %f\n", *alpha);
  fprintf(stdout, "-------------------End-------------------\n");

// freeing memory
  free(y);
  free(x);
  free(Ab);
  free(gb);
  free(gammaVector);
  free(alpha);

  return EXIT_SUCCESS;
}

/*----------------------------------------------------------------------------*/
/* Function add                                                               */
/* Parameters:                                                                */
/*   const double* vector1......first vector to be added                      */
/*   const double* vector2......second vector to be added                     */
/*   int size...................vector size                                   */
/* Return value:                                                              */
/*   double*....................result of the addition                        */
/*----------------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------------*/
/* Function multiply                                                          */
/* Parameters:                                                                */
/*   const double* vector1......first vector to be added                      */
/*   const double* vector2......second vector to be added                     */
/*   int size...................vector size                                   */
/* Return value:                                                              */
/*   double*....................result of the multiplication                  */
/*----------------------------------------------------------------------------*/
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
