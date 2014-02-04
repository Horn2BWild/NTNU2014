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
  int kupper=0;
  double Snpartial=0.0;
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
  
  int ktemp=klower;
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
  for(i=0; i<=kupper; i++)
  {
    Snpartial=0;
    for(j=pow(2,i); j<pow(2,i+1); j++)
    {
      v->data[j]=1/pow(j,2);
      Snpartial+=v->data[j];
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "---------sum calculation--------\n");
    fprintf(stdout, "-- v->data[%d]: %f\n", i, v->data[i]);
    fprintf(stdout, "-- Snpartial: %f\n", Snpartial);
    fprintf(stdout, "--------------------------------\n");
    if(i%10==0)
    {
      getchar();
    }
  */
    }
    Sn += Snpartial;
    if(i>=klower && i<=kupper)
    {
      diff=S-Sn;
      fprintf(stdout, "k=%d\n  elements:%d\n  %lf\n--------------------\n", i, j, diff);
    }
  }
  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "-----------calculation----------\n");
    fprintf(stdout, "-- Sn: %f\n", Sn);
    fprintf(stdout, "-- S: %f\n", S);
    fprintf(stdout, "-- diff: %f\n", diff);
    fprintf(stdout, "--------------------------------\n");
  */

  endTime=WallTime();

  fprintf(stdout, "total run time: %lf\n\n", endTime-startTime);

  close_app();
  return EXIT_SUCCESS;
}

double mathPi()
{
  return 4*atan(1);
}
