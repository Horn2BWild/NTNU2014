#define VECTORSIZE 20

double vector[VECTORSIZE]={0}
int i=0;

for(i=0; i<VECTORSIZE/3; i++)
{
//access each 3rd element
  printf("a[%d]: %f", i, 3*i);
  if(VECTORSIZE%3>0)
  {
    printf("b[%d]: %f", i, 3*i+1);
    if(VECTORSIZE%3>1)
    {
      printf("c[%d]: %f", i, 3*i+2);
    }
  }
}

