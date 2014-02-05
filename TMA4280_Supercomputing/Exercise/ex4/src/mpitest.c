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
  double* sum=(double*)malloc(sizeof(double));
  int *displ, *sublength;
  double* globalsum=(double*)malloc(sizeof(double));
  double* sendvec=(double*)malloc(VECTORSIZE*sizeof(double));

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
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    for(i=0; i<VECTORSIZE;i++)
    {
      fprintf(stdout, "process %d\n  sendvec[%d]=%f\n", rank, i, sendvec[i]);
    }
  */
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    for(i=0; i<size; i++)
    {
      fprintf(stdout, "size: %d, displacement %d\n", sublength[i], displ[i]);
      fprintf(stdout, "address %x\n", &(sendvec[displ[i]]));
    }
  */
    for (i=0; i < size; ++i)
    {
      //fprintf(stdout, "---send for proc %d\n", i);
      double* vsend=&(sendvec[displ[i]]);
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
      for(j=0; j<sublength[i]; j++)
      {
        fprintf(stdout, "----process %d vsend[j]=%f\n", i, vsend[j]);
      }
  */
      MPI_Send(vsend, sublength[i], MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
    }
  }

  receivevec=(double*)malloc(sizeof(double)*sublength[rank]);
  MPI_Recv(receivevec, sublength[rank], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
  fprintf(stdout, "process %d: data received\n", rank);
  for(j=0;j<sublength[rank];j++)
  {
    fprintf(stdout, "process %d\n  element %d: %f\n  sum: %f\n", rank, j, receivevec[j], (*sum));
    (*sum)+=receivevec[j];
  }
  //}
  //MPI_Reduce(sum, globalsum, sizeof(double), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
  MPI_Allreduce(sum, globalsum, sizeof(double), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  printf("sum %f\n", *globalsum);

  close_app();
  return 0;
}
