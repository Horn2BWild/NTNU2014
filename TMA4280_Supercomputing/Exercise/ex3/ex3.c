#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

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
//MPI init
  MPI_Status status;
  MPI_Init(&argc, &argv);

  int proc=0;
  double *y = NULL; //sum vector y
  double *x = NULL; //sum vector x
  double *Ab = NULL; //A multiplied by b
  double *gb = NULL; //b multiplied by gamma
  double *gammaVector = NULL;
  int i=0, j=0;
  int n=0; 
  double gamma=0.0;
 // double alpha=0.0;
  double* alpha = NULL;
  int rank=0, size=0, tag=0; //MPI flags
  int numProcesses=0;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  tag = 100;

  double start_time = MPI_Wtime();

// commandline argument
  if(argc<3)
  {
    fprintf(stdout, "Usage: ex2 <gamma> <numProc>\n");
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

  numProcesses=atoi(argv[2]);
  
// calculation of A*b
  //#pragma Omp parallel for reduce(+:Ab[i]) schedule(guided,1)
  Ab = (double*)malloc(VSIZE*sizeof(double));
  for(i=0; i<VSIZE; i++)
  {
    if (rank == 0) 
    {
       for (proc=1; proc<numProcesses; proc++) //for each process except the 0
       {
         portion = (NUM_ROWS_A/(numProcesses-1)); // calculate submatrix-size
         low_bound = (proc-1)*portion;
         if (((proc+1)==numProcesses) && ((VSIZE%(numProcesses-1))!=0)) //if division not possible equally
         {
           upper_bound=VSIZE; //last process gets remaining rows
         }
         else
         {
           upper_bound=low_bound+portion; //every process with same data
         }

//MPI non-blocking send function calls
         MPI_Isend(&low_bound, 1, MPI_INT, proc, tag, MPI_COMM_WORLD, &request);
         MPI_Isend(&upper_bound, 1, MPI_INT, proc, tag+1, MPI_COMM_WORLD, &request);
         MPI_Isend(&mat_a[low_bound][0], (upper_bound-low_bound)*VSIZE, MPI_DOUBLE, proc, tag+2, MPI_COMM_WORLD, &request);
         }
       }
//broadcast b to all the slaves
       MPI_Bcast(b, VSIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    if (rank > 0)
    {
	  MPI_Recv(&low_bound, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(&upper_bound, 1, MPI_INT, 0, tag+1, MPI_COMM_WORLD, &status);
	  MPI_Recv(&mat_a[low_bound][0], (upper_bound-low_bound)*VSIZE, MPI_DOUBLE, 0, tag+2, MPI_COMM_WORLD, &status);

	  for (i = low_bound; i < upper_bound; i++)
      {
	    for (j = 0; j < NUM_COLUMNS_B; j++)
        {
       	  for (k = 0; k < NUM_ROWS_B; k++) 
          {
	        mat_result[i][j] += (mat_a[i][k] * b[j]);
	      }
	    }
	  }
    }
  }

// gamma*b
  gb = (double*)malloc(VSIZE*sizeof(double));
  gb = multiply(gammaVector, b, VSIZE);

/*---DECOMMENT FOR DEBUGGING PURPOSES---
for(i=0; i<VSIZE; i++)
printf("gb[%d] = %f", i, gb[i]);
*/

/*---DECOMMENT FOR OPENMP---
#pragma omp parallel sections
*/
{
// y = a + Ab
/*---DECOMMENT FOR OPENMP---
#pragma omp parallel section
*/
{
  y = (double*)malloc(VSIZE*sizeof(double));
  y = add(Ab, a, 3);
}

// x = a + gb
/*---DECOMMENT FOR OPENMP---
#pragma omp parallel section
*/
{
  x = (double*)malloc(VSIZE*sizeof(double));
  x = add(gb, a, VSIZE);
}
}

// alpha = xT*y
  alpha = (double*)malloc(sizeof(double));
  alpha = multiply(x, y, VSIZE);


// Outputs
  fprintf(stdout, "-----------------Outputs-----------------\n");
  fprintf(stdout, " x = [%f %f %f]\n", x[0], x[1], x[2]);
  fprintf(stdout, " y = [%f %f %f]\n", y[0], y[1], y[2]);
  fprintf(stdout, " alpha = %f\n", *alpha);
  fprintf(stdout, "-------------------End-------------------\n");

// freeing memory
  free(y);
  free(x);
  free(Ab);
  free(gb);
  free(gammaVector);
  free(alpha);
 
  MPI_Finalize();
  return EXIT_SUCCESS;
}

/*----------------------------------------------------------------------------*/
/* Function add */
/* Parameters: */
/* const double* vector1......first vector to be added */
/* const double* vector2......second vector to be added */
/* int size...................vector size */
/* Return value: */
/* double*....................result of the addition */
/*----------------------------------------------------------------------------*/
//#pragma Omp parallel schedule(guided,1)
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
/* Function multiply */
/* Parameters: */
/* const double* vector1......first vector to be added */
/* const double* vector2......second vector to be added */
/* int size...................vector size */
/* Return value: */
/* double*....................result of the multiplication */
/*----------------------------------------------------------------------------*/
//#pragma Omp parallel schedule(guided,1)
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
