#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

int main(int argc, char** argv)
{
  int vectorlength=pow(2,14);
  double S=pow(pi(),2)/6;
  Vector v=createVector(vectorlength);
  int i=0;

//calculate vector elements
  for(i=0; i<vectorlength; i++)
  {
    v->data[i]=1/pow(i,2);
  }

  return EXIT_SUCCESS;
}

double pi()
{
  return 4*atan(1);
}
