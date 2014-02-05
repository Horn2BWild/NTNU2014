#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define VECTORSIZE 32

int main(int argc, char** argv)
{
  int rank, size, i, tag;
  MPI_Status status;
  char message[20];
  int j=0;
  double sum=0;
  double* sendvec=(double*)malloc(VECTORSIZE*sizeof(double));
  //double* receivevec=(double*)malloc(VECTORSIZE*sizeof(double));
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double* receivevec=(double*)malloc(sizeof(double)*((VECTORSIZE/size)+1));
  tag = 100;

  if (rank == 0) 
  {
    for(i=0;i<VECTORSIZE;i++)
    {
      sendvec[i]=i;
    }

    MPI_Scatter(sendvec, sizeof(sendvec)/(size*sizeof(double)), MPI_DOUBLE, receivevec, sizeof(sendvec)/(size*sizeof(double)), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //  for (i=1; i < size; ++i)
  //    MPI_Send(vec, sizeof(vec), MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
  } 
  else
  {

   // MPI_Recv(receivevec, sizeof(vec), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    for(j=0;j<VECTORSIZE/size;j++)
    {
      sum+=receivevec[j];
    }
  }

  for(i=0;i<size;i++)
  {
    printf("process %d: %f\n", i, sum);
  }

  MPI_Finalize();
  return 0;
}
