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

double mathPi();

int main(int argc, char** argv)
{
  int vectorlength=pow(2,14);
  double S=pow(mathPi(),2)/6;
  double Sn=0.0;
  Vector v=createVector(vectorlength);
  int i=0;
  int j=0;
  double diff=0.0;

  /*---DECOMMENT FOR DEBUGGING PURPOSES---
  fprintf(stdout, "------precalculated values------\n");
  fprintf(stdout, "-- PI: %f\n", mathPi());
  fprintf(stdout, "-- S: %f\n", S);
  fprintf(stdout, "--------------------------------\n");
  */

//calculate vector elements
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

    for(j=3; j<=14; j++)
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


  return EXIT_SUCCESS;
}

double mathPi()
{
  return 4*atan(1);
}
