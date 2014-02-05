/**
* Author: Andreas J. Hoermer
* Email: andreas (at) hoermer.at
*
* Filename: ex4.c
* LastChange: 04.02.2014
**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

#include <omp.h>

//function prototypes
double mathPi();

double sum(double* vec, int length)
{
  int i=0;
  double partsum=0.0;
#pragma omp parallel for schedule(guided,1) reduction(+:partsum)
  for(i=0; i<length; i++)
  {
    partsum+=vec[i];
  }
  return partsum;
}

//! \usage: ex4 <lowerbound> <upperbound>
//e.g. 3..14
int main(int argc, char** argv)
{
  double startTime=WallTime();
  MPI_Status status;
  double endTime=0.0;
  double S=pow(mathPi(),2)/6; //reference value S
  double Sn=0.0; //approximated Sn
  int i=0;
  int j=0;
  double diff=0.0; //difference S-Sn
  int rank=0;
  int size=0;
  int klower=0;
  int P=0;
  int tag=100;
  int kupper=0;
  int dbgloop=0;
  double Snpartial=0.0;
  double* globalsum=(double*)malloc(sizeof(double));
  int *displ, *sublength;
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
fprintf(stdout, "------precalculated values------\n");
fprintf(stdout, "-- PI: %f\n", mathPi());
fprintf(stdout, "-- S: %f\n", S);
fprintf(stdout, "--------------------------------\n");
*/

  if(argc<3 || argc>3)
  {
    fprintf(stdout, "\nusage: ex4 <lowerbound> <upperbound>\nexiting...\n\n");
    return EXIT_FAILURE;
  }

  klower=atoi(argv[1]);
  kupper=atoi(argv[2]);
 // P=atoi(argv[3]);
  
  int ktemp=klower;
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
fprintf(stdout, "------command line arguments------\n");
fprintf(stdout, "-- lower bound: %d\n", klower);
fprintf(stdout, "-- upper bound: %d\n", kupper);
fprintf(stdout, "----------------------------------\n");
*/

  int vectorlength=pow(2,kupper); //maximum vector length
  init_app(argc, argv, &rank, &size);

  double* sendvec=(double*)malloc(vectorlength*sizeof(double));
  double* receivevec;
  double* localsum=(double*)malloc(sizeof(double));
  //Vector v=createVector(vectorlength); //storage vector for sum elements
//calculate vector elements

  if(rank==0)
  {
    for(i=1;i<=kupper; i++)
    {
		for(j=pow(2,i-1); j<pow(2,i); j++) //starting from element 1
		{
		    sendvec[j-1]=1.0/pow(j,2);
		}
    }
    fprintf(stdout,"memory allocated and calculated values\n");
  

    for(i=1;i<=kupper; i++)
    {
	  if(rank==0)
	  {
		for(dbgloop=0;dbgloop<pow(2,i);dbgloop++)
		{
		  fprintf(stdout, "sendvec[%d]: %f\n", dbgloop, sendvec[dbgloop]);
		}
		  splitVector(pow(2,i), size, &sublength, &displ);
		  fprintf(stdout, "k: %d size: %d\n", i, size);
		for (j=0; j < size; j++)
		{

		    fprintf(stdout, "sublength[%d]: %d, displacement[%d]: %d\n", j, sublength[j], j, displ[j]);
		}
		for (j=0; j < size; j++)
		{
			  fprintf(stdout, "---send for proc %d\n", j);

			  double* vsend=&(sendvec[displ[j]]);
			  for(dbgloop=0; dbgloop<sublength[j]; dbgloop++)
			  {
				fprintf(stdout, "----process %d\nvsend[j]=%f\n", j, vsend[dbgloop]);
			  }
		  	  MPI_Send(vsend, sublength[j], MPI_DOUBLE, j, tag, MPI_COMM_WORLD);

		} 
	  }
}

	  receivevec=(double*)malloc(sizeof(double)*sublength[rank]);
	  MPI_Recv(receivevec, sublength[rank], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
	  fprintf(stdout, "process %d: data received\n", rank);
	  (*localsum)=sum(receivevec, sublength[rank]);

	  //MPI_Reduce(sum, globalsum, sizeof(double), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	  MPI_Allreduce(localsum, globalsum, sizeof(double), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	  diff=S-(*globalsum);
	  fprintf(stdout, "k=%d\n elements:%d\n S-Sn:%lf\n--------------------\n", i, j, diff);
    MPI_Barrier(MPI_COMM_WORLD);
    }
 


  endTime=WallTime();

  fprintf(stdout, "total run time: %lf\n\n", endTime-startTime);

  free(sendvec);
  close_app();
  return EXIT_SUCCESS;
}

//! \brief calculating the number pi
//! \return pi
double mathPi()
{
  return 4*atan(1);
}
