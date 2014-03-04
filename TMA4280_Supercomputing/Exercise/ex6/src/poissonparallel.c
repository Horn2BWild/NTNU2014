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


typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);


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

    init_app(argc, argv, &rank, &size);

    //calculate partial vector sizes to send
    scnt=(int*)malloc(size*sizeof(int));
    displ=(int*)malloc(size*sizeof(int));

    int elementcount=m/size;
    for(i=1; i<size; i++)
    {
        scnt[i]=elementcount;
    }
    scnt[0]=m%size;

    //calculate partial displacements to send
    int elementdisplacement=0;
    displ[0]=0;
    for(i=1; i<size; i++)
    {
        elementdisplacement+=scnt[i-1];
        displ[i]=elementdisplacement;
    }


    diag = createRealArray (m);
    b    = createReal2DArray (m,m);
    bt   = createReal2DArray (m,m);
    z    = createRealArray (nn);

    h    = 1./(Real)n;
    pi   = 4.*atan(1.);

#pragma omp parallel for schedule(guided,1)
    for (i=0; i < m; i++)
    {
        diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
    }
#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        for (i=0; i < m; i++)
        {
            b[j][i] = h*h;
        }
    }
#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        fst_(b[j], &n, z, &nn);
    }


  MPI_Alltoallv (&b[0][0],	/* address of data to send  */
		scnt,	/* number of items to send to processes  */
		displ, /* displacements for each process */
		MPI_DOUBLE,	/* type of data  */
		&bt[0][0],	/* address for receiving the data  */
		/* NOTE: send data and receive data may NOT overlap */
		scnt,	/* number of items to receive
				   from any process  */
    displ,
		MPI_DOUBLE,	/* type of receive data  */
		MPI_COMM_WORLD);

    transpose (bt,b,m);
#pragma omp parallel for schedule(guided,1)
    for (i=0; i < m; i++)
    {
        fstinv_(bt[i], &n, z, &nn);
    }

#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        for (i=0; i < m; i++)
        {
            bt[j][i] = bt[j][i]/(diag[i]+diag[j]);
        }
    }
#pragma omp parallel for schedule(guided,1)
    for (i=0; i < m; i++)
    {
        fst_(bt[i], &n, z, &nn);
    }



    transpose (b,bt,m);

#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        fstinv_(b[j], &n, z, &nn);
    }

    umax = 0.0;
#pragma omp parallel for schedule(guided,1)
    for (j=0; j < m; j++)
    {
        for (i=0; i < m; i++)
        {
            if (b[j][i] > umax) umax = b[j][i];
        }
    }
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
