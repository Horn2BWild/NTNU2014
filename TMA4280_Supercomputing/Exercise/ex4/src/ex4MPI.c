/** 
  * Author: Andreas J. Hoermer
  * Email: andreas (at) hoermer.at 
  * 
  * Filename: ex4MPI.c
  * LastChange: 04.02.2014
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

#include <omp.h>

//function prototypes
double mathPi();

double dosum(Vector v)
{
  int i;
  double tempsum=0.0;
  for (i=0;i<v->len;i++) {
    tempsum+=v->data[i];
  }

#ifdef HAVE_MPI
  for (i=0;i<v->comm_size;++i) {
    MPI_Reduce(temp->data+v->displ[i], u->data, v->sizes[i],
               MPI_DOUBLE, MPI_SUM, i, *v->comm);
  }
  freeVector(temp);
#endif
  return alpha;
}


//! \usage: ex4 <lowerbound> <upperbound>
//e.g. 3..14
int main(int argc, char** argv)
{ 
  double startTime=WallTime();
  double endTime=0.0;
  double S=pow(mathPi(),2)/6; //reference value S
  double Sn=0.0; //approximated Sn
  int i=0;
  int j=0;
  double diff=0.0; //difference S-Sn
  int rank=0;
  int size=0;
  int k=0;
  int P=0;
  int *sublength;
  int *displacement;
  double Snpartial=0.0;
  int length=0;
  int tag=100;
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "------precalculated values------\n");
    fprintf(stdout, "-- PI: %f\n", mathPi());
    fprintf(stdout, "-- S: %f\n", S);
    fprintf(stdout, "--------------------------------\n");
  */

  if(argc<3 || argc>3)
  {
    fprintf(stdout, "\nusage: ex4 <k> <P>\nexiting...\n\n");
    return EXIT_FAILURE;
  }

  k=atoi(argv[1]);
  P=atoi(argv[2]);
  
  int ktemp=k;
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "------command line arguments------\n");
    fprintf(stdout, "-- k %d\n", k);
    fprintf(stdout, "----------------------------------\n");
  */

  int vectorlength=pow(2,k); //maximum vector length
  init_app(argc, argv, &rank, &size);


  Vector v=createVectorMPI(vectorlength, &WorldComm, 0);

  if(rank==0)
  {
  //generating vector elements
  for(i=1; i<=vectorlength; i++)
  {
      v->data[i-1]=1.0/pow(i,2); 
  }

  sum = dosum(v);


  //partitioning vector elements
  for(i=0; i<P; i++)
  {
  //distribution
    //distributing length
    MPI_send(sublength, 1, MPI_DOUBLE, i, tag, &WorldComm);
    //distributing subvector
    MPI_Send((v->data)), sublength, MPI_DOUBLE, i, tag, &WorldComm);
  }
  //collecting and sum up
  MPI_Reduce (&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, &WorldComm);
  //printout
  }
  else
  {
  //receiving
    //receiving length
    MPI_Recv(&length, 1, MPI_CHAR, 0 ,tag, &WorldComm,&status);
    //receiving subvector
    double* partialVector=(double*)malloc(length*sizeof(double));
    MPI_Recv(partialVector, length, MPI_DOUBLE, 0 ,tag, &WorldComm,&status);
  //calculate partial sums
    Sn=0;
    for(i=0; i<=length; i++)
  {
      Sn+=partialVector[i]

  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "---------sum calculation--------\n");
    fprintf(stdout, "-- v->data[%d]: %f\n", j, v->data[j-1]);
    fprintf(stdout, "-- Sn: %f\n", Sn);
    fprintf(stdout, "--------------------------------\n");
    if(i%10==0)
    {
      getchar();
    }
  */
  }
//sending back partial sum
  }

  /*---DECOMMENT FOR DEBUGGING PURPOSES--- 
    fprintf(stdout, "-- Vector length %d\n", vectorlength);
    fprintf(stdout, "----------------------------------\n");
  */

//calculate vector elements
 // #pragma omp parallel for schedule(guided,1) reduction(+:Sn)

  diff=S-Sn;
  fprintf(stdout, "k=%d\n  elements:%d\n  %lf\n--------------------\n", k, i-1, diff);
 
  endTime=WallTime();
  fprintf(stdout, "total run time: %lf\n\n", endTime-startTime);

  freeVector(v);
  close_app();
  return EXIT_SUCCESS;
}

//! \brief calculating the number pi
//! \return pi
double mathPi()
{
  return 4*atan(1);
}
