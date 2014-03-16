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
#define DEBUG_TESTMATRIX 0
#define DEBUG_MPISUMMARY 0
#define DEBUG_ALL 0
#define DEBUG_MPIINFO 0
#define DEBUG_TESTTRANSPOSE 0


typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

void trans (double *a, int n);


void mergeAndPrintMpiMat(Real **b, int size, int rank, int m, int* scnt, int* displ);
void printMat(Real **b, int m);

void mergeAndPrintMpiMat(Real **b, int size, int rank, int m, int* scnt, int* displ)
{
    MPI_Status status;
    int tag = 10;

    MPI_Send(b[0], scnt[rank]*m, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

    if (rank == 0)
    {
        Real** test = createReal2DArray(m,m);

        int i = 0;
        for (i = 0; i < size; ++i)
        {
            MPI_Recv(test[displ[i]], scnt[i]*m, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
        }
        printMat(test, m);
    }
}

void printMat(Real **b, int m)
{
    printf("\n-----------------------\n");
    int i,j;
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            printf("%f ", b[j][i]);
        }
        printf("\n");
    }
    printf("-----------------------\n");

}

void transposeMPI(Real **bt, Real**b, int m, int rank, int size, int* scnt, int* displ)
{
    int i=0;
    int proccnt=0;
    int rowcnt=0;
    int elementcnt=0;
    int columncnt=0;
    int vectorposition=0;
    Real *sendvector=malloc(m*scnt[rank]*sizeof(Real));
    Real *receivevector=malloc(m*scnt[rank]*sizeof(Real));

    for(rowcnt=0; rowcnt<m; rowcnt++)
    {
        for(columncnt=0; columncnt<scnt[rank]; columncnt++)
        {
            sendvector[vectorposition]=b[columncnt][rowcnt];
            vectorposition++;
        }
    }


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
        MPIscnt[i]=scnt[i]*scnt[rank];
    }
#if DEBUG_MPIINFO
    for(i=0; i<size; i++)
    {
        fprintf(stdout, "rank %d: MPIdispl[%d]: %d\n", rank, i, MPIdispl[i]);
    }
        for(i=0; i<size; i++)
    {
                fprintf(stdout, "rank %d: MPIscnt[%d]: %d\n", rank, i, MPIscnt[i]);
    }
    #endif



#if DEBUG_MPISUMMARY
fprintf(stdout, "--------------------------------------\n");
fprintf(stdout, "sendvector:");
for(i=0; i<vectorposition; i++)
{
  fprintf(stdout, "%f\t", sendvector[i]);
}
fprintf(stdout, "\n");
fprintf(stdout, "MPIscnt:");
for(i=0; i<size; i++)
{
  fprintf(stdout, "%d\t", MPIscnt[i]);
}
fprintf(stdout, "\n");
fprintf(stdout, "MPIdispl:");
for(i=0; i<size; i++)
{
  fprintf(stdout, "%d\t", MPIdispl[i]);
}
fprintf(stdout, "\n");

fprintf(stdout, "rank: %d before alltoallv", rank);
fflush(stdout);
#endif


    MPI_Alltoallv (
        sendvector,	/* address of data to send  */
        MPIscnt,	/* number of items to send to processes  */
        MPIdispl, /* displacements for each process */
        MPI_DOUBLE,	/* type of data  */
        receivevector,	/* address for receiving the data  */
        /* NOTE: send data and receive data may NOT overlap */
     MPIscnt,	/* number of items to receive
          from any process  */
        MPIdispl,
        MPI_DOUBLE,	/* type of receive data  */
        MPI_COMM_WORLD);

#if DEBUG_TESTMATRIX
  //  vectorposition=m*scnt[rank];
    fprintf(stdout, "receive vector size proc %d: %d\n", rank, vectorposition);
    fprintf(stdout, "---RECEIVEVECTOR process %d---\n", rank);
    for(i=0; i<vectorposition; i++)
    {
        fprintf(stdout, "%d.%d: %.1f\t",rank, i, receivevector[i]);
        fflush(stdout);
    }

#endif
//   fprintf(stdout, "proc %d step 9\n", rank);

    vectorposition = 0;
    for (proccnt = 0; proccnt < size; ++proccnt)
    {
        for(columncnt=0; columncnt<scnt[rank];columncnt++)
        {
            for (rowcnt = displ[proccnt]; rowcnt < displ[proccnt] + scnt[proccnt]; ++rowcnt)
            {

                bt[columncnt][rowcnt]=receivevector[vectorposition];
                vectorposition++;
               // fprintf(stdout, "%f\t", bt[columncnt][rowcnt]);
            }
        }
    }
}


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


    splitVector(m, size, &scnt, &displ);


    if(rank==0)
    {
        fprintf(stdout, "proc\tsize\tdispl\n");
        for(i=0; i<size; i++)
        {
            fprintf(stdout, "%d\t%d\t%d\n", i, scnt[i], displ[i]);
        }
        fprintf(stdout, "-------------------\n");
    }
    diag = createRealArray (m);
    b    = createReal2DArray (scnt[rank],m);
    bt   = createReal2DArray (scnt[rank],m);
    z    = createRealArray (nn);

    h    = 1./(Real)n;
    pi   = 4.*atan(1.);
//   fprintf(stdout, "step 3\n");
 #pragma omp parallel for schedule(guided,1)
    for (i=0; i < m; i++)
    {
        diag[i] = 2.*(1.-cos((i+1) * pi/(Real)n));
    }

//   fprintf(stdout, "step 4\n");
 #pragma omp parallel for schedule(guided,1)
    for (j=0; j < scnt[rank]; j++)
    {
        for (i=0; i < m; i++)
        {
            b[j][i] = h*h*5*M_PI*M_PI*(
                sin(M_PI*(j+1+displ[rank])*h)*
                sin(2*M_PI*(i+1)*h));
#if DEBUG_TESTTRANSPOSE
            b[j][i] = i*m + j + displ[rank];
#endif
        }
    }

#if DEBUG_TESTTRANSPOSE
//********** testing stuff
    mergeAndPrintMpiMat(b, size, rank, m, scnt, displ);
    transposeMPI(bt,b,m,rank,size, scnt, displ);
    mergeAndPrintMpiMat(bt, size, rank, m, scnt, displ);

    close_app();
    return EXIT_SUCCESS;
//********** end testing stuff
#endif



    //     fprintf(stdout, "step 5\n");
// #pragma omp parallel for schedule(guided,1)
    for (j=0; j < scnt[rank]; j++)
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
    if(rank==0)
    {
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

    // mergeAndPrintMpiMat(b, size, rank, m, scnt, displ);
transposeMPI(bt,b,m,rank,size, scnt, displ);
    // mergeAndPrintMpiMat(bt, size, rank, m, scnt, displ);

// #pragma omp parallel for schedule(guided,1)
    for (i=0; i < scnt[rank]; i++)
    {
        fstinv_(bt[i], &n, z, &nn);
    }
 //  fprintf(stdout, "proc %d step 11\n", rank);
// #pragma omp parallel for schedule(guided,1)
    for (j=0; j < scnt[rank]; j++)
    {
        for (i=0; i < m; i++)
        {
            bt[j][i] = bt[j][i]/(diag[i]+diag[j+displ[rank]]);
        }
    }
   //      fprintf(stdout, "proc %d step 12\n", rank);
// #pragma omp parallel for schedule(guided,1)
    for (i=0; i < scnt[rank]; i++)
    {
        fst_(bt[i], &n, z, &nn);
    }

 //  fprintf(stdout, "proc %d step 13\n", rank);

transposeMPI(b,bt,m,rank,size, scnt, displ);

// #pragma omp parallel for schedule(guided,1)
    for (j=0; j < scnt[rank]; j++)
    {
        fstinv_(b[j], &n, z, &nn);
    }
  //    fprintf(stdout, "proc %d step 14\n", rank);
    umax = 0.0;
 #pragma omp parallel for schedule(guided,1)
    for (j=0; j < scnt[rank]; j++) {
        for (i=0; i < m; i++) {
            double value = fabs(b[j][i] -
                sin(M_PI * (j+1+displ[rank])*h)*
                sin(2*M_PI*(i+1)*h));
            if (value > umax) umax = value;
        }
    }
    // for (j=0; j < m; j++)
    // {
    //     for (i=0; i < m; i++)
    //     {
    //         if (b[j][i] > umax) umax = b[j][i];
    //     }
    // }
    double globalMax = 0;
    MPI_Reduce(&umax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        fprintf(stdout, "proc %d done...\n", rank);
        printf (" umax = %f \n",globalMax);
    }
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
