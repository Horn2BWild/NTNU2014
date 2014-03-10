/**************************************************************************** *
* Author: Andreas J. Hoermer                                                  *
* Email: andreas (at) hoermer.at                                              *
*                                                                             *
* Filename: poissonparallel.c                                                 *
* LastChange: 04.03.2014                                                      *
*                                                                             *
* --------------------------------------------------------------------------- *
* PROGRAM DESCRIPTION                                                         *
* --------------------------------------------------------------------------- *
*                                                                             *
*                                                                             *
*                                                                             *
*                                                                             *
*                                                                             *
* --------------------------------------------------------------------------- *
* USAGE                                                                       *
* --------------------------------------------------------------------------- *
* ./ex6parallel ...                                                           *
*                                                                             *
*                                                                             *
*                                                                             *
*                                                                             *
* *************************************************************************** */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "common.h"

#define DEBUG_2KCHECK 0
#define DEBUG_TESTMATRIX 1
#define DEBUG_ALL 0


typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

void trans (double *a, int n);


int main(int argc, char **argv)
{
    int rank=0;                                     //current proc id
    int size=0;                                     //number of proc
    int tag=100;                                    //tag used for MPI comm
    MPI_Status status;
    Real *diag, **b, **bt, *z;
    Real pi, h, umax;
    int i, j, n, m, nn;
    int *displ;           // displacements
    int *scnt;            // send/receive count

    /* the total number of grid points in each spatial direction is (n+1) */
    /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
    /* this version requires n to be a power of 2 */

    if(argc < 2)
    {
        printf("need a problem size\n");
        return EXIT_FAILURE;
    }


    n  = atoi(argv[1]);
    m  = n-1;
    nn = 4*n;

    //check for 2^k
    unsigned int tempn=(unsigned int)n;
    int onecount=0;
    while(tempn>0)
    {
#if DEBUG_2KCHECK
        fprintf(stdout, "-inside while\n");
        fflush(stdout);
#endif
        if((tempn&1)!=0)
        {
#if DEBUG_2KCHECK
            fprintf(stdout, "--inside if\n");
            fflush(stdout);
#endif
            onecount++;
            if(onecount>1)
            {
                fprintf(stdout, "n has to be 2^k\n");
                fflush(stdout);
                return EXIT_FAILURE;
            }
        }
#if DEBUG_2KCHECK
        fprintf(stdout, "--outside if\n");
        fflush(stdout);
#endif
        tempn=tempn>>1;
#if DEBUG_2KCHECK
        fprintf(stdout, "tempn: %d\n",tempn);
        fflush(stdout);
#endif
    }
#if DEBUG_2KCHECK
    fprintf(stdout, "-outside while\n");
    fflush(stdout);
#endif
  //  fprintf(stdout, "step 1\n");
    init_app(argc, argv, &rank, &size);

    //calculate partial vector sizes to send
//    scnt=(int*)malloc(size*sizeof(int));
//    displ=(int*)malloc(size*sizeof(int));

//    int elementcount=m/size;
//    for(i=1; i<size; i++)
//    {
//        scnt[i]=elementcount;
//    }
//    scnt[0]=m%size;
 //   displ[0]=0;
splitVector(m, size, &scnt, &displ);
//    fprintf(stdout, "step 2\n");

    //calculate partial displacements to send
 //   int elementdisplacement=0;
 //   displ[0]=0;
 //   for(i=1; i<size; i++)
 //   {
 //       elementdisplacement+=scnt[i-1];
  //      displ[i]=elementdisplacement;
 //   }

if(rank==0){
fprintf(stdout, "proc\tsize\tdispl\n");
for(i=0; i<size; i++)
{
  fprintf(stdout, "%d\t%d\t%d\n", i, scnt[i], displ[i]);
}
fprintf(stdout, "-------------------\n");
}
    diag = createRealArray (m);
    b    = createReal2DArray (m,m);
    bt   = createReal2DArray (m,m);
    z    = createRealArray (nn);

    h    = 1./(Real)n;
    pi   = 4.*atan(1.);
 //   fprintf(stdout, "step 3\n");
#pragma omp parallel for schedule(guided,1)
    for (i=0; i < m; i++)
    {
        diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
    }

 //   fprintf(stdout, "step 4\n");
#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        for (i=0; i < m; i++)
        {
            b[j][i] = h*h;
        }
    }
   //     fprintf(stdout, "step 5\n");
#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        fst_(b[j], &n, z, &nn);
    }
  //  fprintf(stdout, "step 6\n");
#if DEBUG_TESTMATRIX
int value=0;
for(i=0; i<m; i++)
{
  for(j=0; j<m; j++)
  {
    value++;
    b[i][j]=value;
  }
}
  //  fprintf(stdout, "step 7\n");
if(rank==0){
fprintf(stdout, "\n-----------------------\n");
for(i=0; i<m; i++)
{
  for(j=0; j<m; j++)
  {
fprintf(stdout, "%.1f\t", b[i][j]);
  }
  fprintf(stdout, "\n");
}
fprintf(stdout, "-----------------------\n");
}
#endif
  //  fprintf(stdout, "step 8\n");


int proccnt=0;
int rowcnt=0;
int elementcnt=0;
int columncnt=0;
int vectorposition=0;
Real *sendvector=malloc(m*scnt[rank]*sizeof(Real));
Real *receivevector=malloc(m*scnt[rank]*sizeof(Real));

fprintf(stdout, "size of vector on process %d: %d\n", rank, m*scnt[rank]);

  for(rowcnt=0; rowcnt<m; rowcnt++)
  {
//for(proccnt=0; proccnt<size; proccnt++)
//{
    for(elementcnt=displ[rank]; elementcnt<displ[rank]+scnt[rank]; elementcnt++)
    {

//fprintf(stdout, "row: %d element %d rank %d: %f\n", rowcnt, elementcnt, rank, b[rowcnt][elementcnt]);
      sendvector[vectorposition]=b[rowcnt][elementcnt];
      vectorposition++;
    }
  }
//}
/*vectorposition=0;
for(proccnt=0; proccnt<size; proccnt++)
{
for(columncnt=0; columncnt<scnt[rank]; columncnt++)
{
for(rowcnt=0; rowcnt<scnt[proccnt]; rowcnt++)
{
fprintf(stdout, "%.1f\t", sendvector[vectorposition]);
vectorposition++;
}
fprintf(stdout, "\n");
}
}
*/


fprintf(stdout, "---SENDVECTOR process %d---\n", rank);
for(i=0; i<vectorposition; i++)
{
  fprintf(stdout, "%d.%d: %.1f\t",rank, i, sendvector[i]);
}
fprintf(stdout, "\n");

    int *MPIdispl=malloc(size*sizeof(int));           // displacements
    int *MPIscnt=malloc(size*sizeof(int));            // send/receive count

for(i=0; i<size; i++)
{
MPIdispl[i]=0;
MPIscnt[i]=0;
}

for(i=1; i<size; i++)
{
//fprintf(stdout, "scnt[%d]: %d, scnt[%d]: %d\n", rank, scnt[rank], i, scnt[i]);
MPIdispl[i]=MPIdispl[i-1]+scnt[rank]*scnt[i-1];
}

for(i=0; i<size; i++)
{
fprintf(stdout, "rank %d: MPIdispl[%d]: %d\n", rank, i, MPIdispl[i]);
}
for(i=0; i<size; i++)
{
MPIscnt[i]=scnt[i]*scnt[rank];
}
for(i=0; i<size; i++)
{
//fprintf(stdout, "MPI p%d to %d: scnt: %d displ %d\n", i, rank, MPIscnt[i], MPIdispl[i]);
}

  MPI_Alltoallv (
    &sendvector,	/* address of data to send  */
		MPIscnt,	/* number of items to send to processes  */
		MPIdispl, /* displacements for each process */
		MPI_DOUBLE,	/* type of data  */
		&receivevector,	/* address for receiving the data  */
		/* NOTE: send data and receive data may NOT overlap */
		MPIscnt,	/* number of items to receive
				   from any process  */
    MPIdispl,
		MPI_DOUBLE,	/* type of receive data  */
		MPI_COMM_WORLD);

 //   transpose (bt,b,m);

   #if DEBUG_TESTMATRIX
vectorposition=m*scnt[rank];
fprintf(stdout, "receive vector size proc %d\n", vectorposition);
fprintf(stdout, "---RECEIVEVECTOR process %d---\n", rank);
for(i=0; i<vectorposition; i++)
{
  fprintf(stdout, "%d.%d: %.1f\t",rank, i, receivevector[i]);
fflush(stdout);
}

#endif
 //   fprintf(stdout, "step 9\n");
   for (i = 0; i <size; i++){
      trans (&bt[displ[i]][0], scnt[i]);
   }

   #if DEBUG_TESTMATRIX
   if(rank==0){
   fprintf(stdout, "\n-----------------------\n");
for(i=0; i<m; i++)
{
  for(j=0; j<m; j++)
  {
fprintf(stdout, "%.1f\t", bt[i][j]);
  }
  fprintf(stdout, "\n");
}
fprintf(stdout, "-----------------------\n");
   }
#endif
//s    fprintf(stdout, "step 10\n");
#pragma omp parallel for schedule(guided,1)
    for (i=0; i < m; i++)
    {
        fstinv_(bt[i], &n, z, &nn);
    }
 //   fprintf(stdout, "step 11\n");
#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        for (i=0; i < m; i++)
        {
            bt[j][i] = bt[j][i]/(diag[i]+diag[j]);
        }
    }
   //     fprintf(stdout, "step 12\n");
#pragma omp parallel for schedule(guided,1)
    for (i=0; i < m; i++)
    {
        fst_(bt[i], &n, z, &nn);
    }

 //   fprintf(stdout, "step 13\n");

    transpose (b,bt,m);

#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        fstinv_(b[j], &n, z, &nn);
    }
  //  fprintf(stdout, "step 14\n");
    umax = 0.0;
#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        for (i=0; i < m; i++)
        {
            if (b[j][i] > umax) umax = b[j][i];
        }
    }
 //       fprintf(stdout, "step 15\n");
    printf (" umax = %e \n",umax);
    close_app();
    return EXIT_SUCCESS;
}

void transpose (Real **bt, Real **b, int m)
{
    int i, j;
    for (j=0; j < m; j++)
    {
        for (i=0; i < m; i++)
        {
            bt[j][i] = b[i][j];
        }
    }

}



Real *createRealArray (int n)
{
    Real *a;
    int i;
    a = (Real *)malloc(n*sizeof(Real));
    for (i=0; i < n; i++)
    {
        a[i] = 0.0;
    }
    return (a);
}

Real **createReal2DArray (int n1, int n2)
{
    int i, n;
    Real **a;
    a    = (Real **)malloc(n1   *sizeof(Real *));
    a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
    for (i=1; i < n1; i++)
    {
        a[i] = a[i-1] + n2;
    }
    n = n1*n2;
    memset(a[0],0,n*sizeof(Real));
    return (a);
}

void trans (double *a, int n)
/* transpose square matrix a, dimension nxn */

{
  int i, j;
  int ij, ji, l;
  double tmp;
  ij = 0;
  l = -1;
  for (i = 0; i < n; i++)
    {
      l += n + 1;
      ji = l;
      ij += i + 1;
      for (j = i+1; j < n; j++)
	{
	  tmp = a[ij];
	  a[ij] = a[ji];
	  a[ji] = tmp;
	  ij++;
	  ji += n;
	}
    }
}
