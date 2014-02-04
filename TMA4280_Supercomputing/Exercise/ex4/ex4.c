#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

int main(int argc, char** argv)
{
  int vectorlength=pow(2,14);
  double S=pow(pi(),2)/6;
  double Sn=0.0;
  Vector v=createVector(vectorlength);
  int i=0;
  int j=0;
  double diff=0.0;

//calculate vector elements
  for(i=0; i<vectorlength; i++)
  {
    v->data[i]=1/pow(i,2);
    Sn += v->data[i];

    for(j=3; j<=14; j++)
    {
      if(i==pow(2,j))
      {
        diff=S-Sn;
        fprintf(stdout, "k=%d:\t%d", j, diff);
      }
    }
  }


  return EXIT_SUCCESS;
}

double pi()
{
  return 4*atan(1);
}
