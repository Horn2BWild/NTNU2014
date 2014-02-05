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
  double* vec=(double*)malloc(VECTORSIZE*sizeof(double));

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  tag = 100;

  if (rank == 0) 
  {
    for(i=0;i<VECTORSIZE;i++)
    {
      vec[i]=1.0/pow((i+1),2);
    }
    for (i=1; i < size; ++i)
      MPI_Send(vec, sizeof(vec), MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
  } 
  else
  {
    MPI_Recv(vec, sizeof(vec), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    for(j=0;j<VECTORSIZE;j++)
    {
      vec[j]+=rank;
    }
  }

  for(i=0;i<VECTORSIZE;i++)
  {
    printf("process %d: %s\n", rank, message);
  }

  MPI_Finalize();
  return 0;
}
