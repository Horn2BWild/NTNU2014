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
  int **sublength;
  int **displacement;
  double Snpartial=0.0;
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

#ifdef HAVE_MPI
  Vector v=createVectorMPI(vectorlength, &WorldComm, 0);
#else
  Vector v=createVector(vectorlength); //storage vector for sum elements
#endif

  if(rank==0)
  {
  //generating vector elements
  for(i=1; i<=vectorlength; i++)
  {
      //calculating j, storing in j-1
      //e.g. calculating 1st element, storing in data[0]
      //otherwise buffer overflow at last element
      v->data[i-1]=1.0/pow(i,2); 
  }
  //partitioning vector elements
  splitVector(vectorlength, P, sublength, displacement);
  for(i=0; i<P; i++)
  {
  //distributing vector elements
  MPI_Send(
  }
  //collecting and sum up

  //printout
  }
  else
  {
  //calculate partial sums
  
  }

  /*---DECOMMENT FOR DEBUGGING PURPOSES--- 
    fprintf(stdout, "-- Vector length %d\n", vectorlength);
    fprintf(stdout, "----------------------------------\n");
  */

//calculate vector elements
 // #pragma omp parallel for schedule(guided,1) reduction(+:Sn)
  for(i=1; i<=vectorlength; i++)
  {
      Sn+=v->data[i-1];

  /*---DECOMMENT FOR DEBUGGING PURPOSES---*/
    fprintf(stdout, "---------sum calculation--------\n");
    fprintf(stdout, "-- v->data[%d]: %f\n", j, v->data[j-1]);
    fprintf(stdout, "-- Sn: %f\n", Sn);
    fprintf(stdout, "--------------------------------\n");
    if(i%10==0)
    {
      getchar();
    }
 /* */
  }
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
