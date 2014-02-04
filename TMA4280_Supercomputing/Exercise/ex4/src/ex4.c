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
  #pragma omp parallel for schedule(guided,1) reduction(+:Sn)
  for(i=1; i<=kupper; i++)
  {
    Snpartial=0.0;
 //   #pragma omp parallel for schedule(guided,1) reduction(+:Snpartial)
    for(j=pow(2,i-1); j<pow(2,i); j++) //starting from element 1
    {
      //calculating j, storing in j-1
      //e.g. calculating 1st element, storing in data[0]
      //otherwise buffer overflow at last element
      v->data[j-1]=1/pow(j,2); 
      Snpartial+=v->data[j-1];

  /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "---------sum calculation--------\n");
    fprintf(stdout, "-- v->data[%d]: %f\n", j, v->data[j-1]);
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
