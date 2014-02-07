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
#include <mpi.h>

//function prototypes
double mathPi();
/*NOT USED AT THE MOMENT
int cleanup(void** memory, int length, char* message, int exitcode)
{
    int i=0;
    for(i=0; i<length; i++)
    {
        free(memory[i]);
    }
    close_app();
    if(exitcode==EXIT_FAILURE)
    {
        fprintf(stderr, "%s\naborting...\n", message);
        return EXIT_FAILURE;
    }
    else
    {
        fprintf(stdout, "%s\n", message);
        return EXIT_SUCCESS;
    }
}
*/

double sum(double* vec, int length)
{
    int i=0;
    double partsum=0.0;
#pragma omp parallel for schedule(guided,1) reduction(+:partsum)
    for(i=0; i<length; i++)
    {
        partsum+=vec[i];
    }
    return partsum;
}

//! \usage: ex4 <lowerbound> <upperbound>
//e.g. 3..14
int main(int argc, char** argv)
{
    double startTime=WallTime();
    MPI_Status status;
    double endTime=0.0;
    double S=pow(mathPi(),2)/6; //reference value S
    double Sn=0.0; //approximated Sn
    int i=0;
    int j=0;
    double diff=0.0; //difference S-Sn
    int rank=0;
    int size=0; //no proc
    int klower=0; //lower bound for 2^k
    int tag=100;
    int kupper=0; //upper bound for 2^k calculation
    int dbgloop=0;
    double Snpartial=0.0;
    double* globalsum=(double*)malloc(sizeof(double));
    int *displ, *sublength;
    /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "------precalculated values------\n");
    fprintf(stdout, "-- PI: %f\n", mathPi());
    fprintf(stdout, "-- S: %f\n", S);
    fprintf(stdout, "--------------------------------\n");
    */

    if(argc<3 || argc>3)
    {
        fprintf(stdout, "\nusage: ex4 <lowerbound> <upperbound>\nexiting...\n\n");
        free(globalsum);
        return EXIT_FAILURE;
    }
    init_app(argc, argv, &rank, &size);
    klower=atoi(argv[1]);
    kupper=atoi(argv[2]);
    int vectorlength=pow(2,kupper); //maximum vector length 2^k
    fprintf(stdout, "vectorlength: %d\n", vectorlength);

    double* sendvec=(double*)malloc(vectorlength*sizeof(double));
    double* receivevec=(double*)malloc(vectorlength*sizeof(double));
    double* localsum=(double*)malloc(sizeof(double));

// P=atoi(argv[3]);
    /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "------command line arguments------\n");
    fprintf(stdout, "-- lower bound: %d\n", klower);
    fprintf(stdout, "-- upper bound: %d\n", kupper);
    fprintf(stdout, "----------------------------------\n");
    */



    if(sendvec==NULL)
    {
        fprintf(stdout, "could not allocate memory: sendvec\n aborting...");
        return EXIT_FAILURE;
    }
    if(localsum==NULL)
    {
        fprintf(stdout, "could not allocate memory: localsum\n aborting...");
        return EXIT_FAILURE;
    }
//calculate vector elements

    if(rank==0)
    {
#pragma omp parallel for schedule(guided,1)
        for(i=1; i<=vectorlength; i++)
        {
            sendvec[i-1]=1.0/pow(i,2);
        }
    }

    for(i=klower; i<=kupper; i++)
    {
        splitVector(pow(2,i), size, &sublength, &displ); //split vector for every k

        MPI_Scatterv(sendvec, sublength, displ, MPI_DOUBLE, receivevec, sizeof(double)*vectorlength, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        *localsum=0;
        if(rank>0)
        {
            for(j=0;j<sublength[rank];j++)
            {
                fprintf(stdout, "receivevec[%d]: %f\t", j, receivevec[j]);
                if(j%5==0)
                {
                    fprintf(stdout,"\n");
                }
            }
            fprintf(stdout, "\n");
    //        (*localsum)=sum(receivevec, sublength[rank]);
        }


        MPI_Barrier(MPI_COMM_WORLD);
        //summing up all local sums
        MPI_Allreduce(localsum, globalsum, sizeof(double), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if(rank==0 && i>=klower && i<=kupper)
        {
            diff=S-(*globalsum);
            fprintf(stdout,"----------Result----------\n");
            fprintf(stdout,"--k=%d elements=%d\n", i, (int)pow(2,i));
            fprintf(stdout,"--Sn=%f\n", *globalsum);
            fprintf(stdout,"--S=%f\n", S);
            fprintf(stdout,"----diff=%f\n", diff);


            endTime=WallTime();


            fprintf(stdout, "total run time: %lf\n\n", endTime-startTime);
            fprintf(stdout,"--------------------------\n");
        }
    }
/*
    free(receivevec);
    free(sendvec);
    free(localsum);
    free(globalsum);
    */
    close_app();


    return EXIT_SUCCESS;
}

//! \brief calculating the number pi
//! \return pi
double mathPi()
{
    return 4*atan(1);
}
