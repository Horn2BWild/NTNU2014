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
  int klower=0;
  int kuppper=0;

  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "------precalculated values------\n");
    fprintf(stdout, "-- PI: %f\n", mathPi());
    fprintf(stdout, "-- S: %f\n", S);
    fprintf(stdout, "--------------------------------\n");
  */

  if(argc<3)
  {
    fprintf("\nusage: ex4 <lowerbound> <upperbound>\nexiting...\n\n");
    return EXIT_FAILURE;
  }

  klower=atoi(argv[1]);
  kupper=atoi(argv[2]);

  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "------command line arguments------\n");
    fprintf(stdout, "-- lower bound: %d\n", klower);
    fprintf(stdout, "-- upper bound: %d\n", kupper);
    fprintf(stdout, "----------------------------------\n");
  */

  int vectorlength=pow(2,kupper); //maximum vector length

  init_app(argc, argv, &rank, &size);

  Vector v=createVector(vectorlength); //storage vector for sum elements


//calculate vector elements
  #pragma omp parallel for schedule(guided,1) reduction(+:Sn)
  for(i=1; i<vectorlength; i++)
  {
    v->data[i]=1/pow(i,2);
    Sn += v->data[i];
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "---------sum calculation--------\n");
    fprintf(stdout, "-- v->data[%d]: %f\n", i, v->data[i]);
    fprintf(stdout, "-- Sn: %f\n", Sn);
    fprintf(stdout, "--------------------------------\n");
    if(i%10==0)
    {
      getchar();
    }
  */

    for(j=klower; j<=kupper; j++)
    {
      if(i==pow(2,j)-1)
      {
        diff=S-Sn;
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "-----------calculation----------\n");
    fprintf(stdout, "-- Sn: %f\n", Sn);
    fprintf(stdout, "-- S: %f\n", S);
    fprintf(stdout, "-- diff: %f\n", diff);
    fprintf(stdout, "--------------------------------\n");
  */

        fprintf(stdout, "k=%d\n  elements:%d\n  %lf\n--------------------\n", j, i+1, diff);
      }
    }
  }
  endTime=WallTime();

  fprintf(stdout, "total run time: %f\n\n", endTime-startTime);

  close_app();
  return EXIT_SUCCESS;
}

double mathPi()
{
  return 4*atan(1);
}
