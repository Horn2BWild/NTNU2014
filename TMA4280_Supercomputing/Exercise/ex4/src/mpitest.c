#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "common.h"

#define VECTORSIZE 8

int main(int argc, char** argv)
{
  int rank, size, i, tag;
  MPI_Status status;
  char message[20];
  int j=0;
  double sum=0;
  int *displ, *sublength;
  double* sendvec=(double*)malloc(VECTORSIZE*sizeof(double));
  //double* receivevec=(double*)malloc(VECTORSIZE*sizeof(double));
  //MPI_Init(&argc, &argv);
  //MPI_Comm_size(MPI_COMM_WORLD, &size);
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double* receivevec=(double*)malloc(sizeof(double)*((VECTORSIZE/size)+1));
  tag = 100;

  init_app(argc, argv, &rank, &size);
  splitVector(VECTORSIZE, size, &sublength, &displ);

  if (rank == 0) 
  {
    for(i=0;i<VECTORSIZE;i++)
    {
      sendvec[i]=i;
    }
  /*---DECOMMENT FOR DEBUGGING PURPOSES---*/
    for(i=0; i<VECTORSIZE;i++)
    {
      fprintf(stdout, "process %d\n  sendvec[%d]=%f\n", rank, i, sendvec[i]);
    }
  /**/
 //   MPI_Scatter(sendvec, VECTORSIZE/size, MPI_DOUBLE, receivevec, VECTORSIZE/size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(i=0; i<size; i++)
    {
      fprintf(stdout, "size: %d, displacement %d\n", sublength[i], displ[i]);
    }
    for (i=0; i < size; ++i)
    {
      fprintf(stdout, "---send for proc %d\n", i);
      double* vsend=&(sendvec[displ[i]]);
      for(j=0; j<sublength[i]; j++)
      {
        fprintf(stdout, "----process %d vsend[j]=%f", i, vsend[j]);
      }
      MPI_Send(vsend, sublength[i], MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
    }
  }
 // else
 // {

    MPI_Recv(receivevec, sublength[i], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    fprintf(stdout, "process %d: data received", rank);
    for(j=0;j<sublength[i];j++)
    {
      fprintf(stdout, "process %d\n  element %d: %f\n  sum: %f\n\n", rank, j, receivevec[j], sum);
      sum+=receivevec[j];
    }
  //}

  for(i=0;i<size;i++)
  {
    printf("process %d: %f\n\n", i, sum);
  }

  close_app();
  return 0;
}
