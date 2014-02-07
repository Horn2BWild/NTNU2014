/**
* Author: Andreas J. Hoermer
* Email: andreas (at) hoermer.at
*
* Filename: ex4.c
* LastChange: 07.02.2014
**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

#include <omp.h>
#include <mpi.h>

//function prototypes
double mathPi();

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
    int i=0;
    int j=0;
    int rank=0;                                     //current proc id
    int size=0;                                     //number of proc
    int klower=0;                                   //lower bound for 2^k
    int tag=100;                                    //tag used for MPI comm
    int kupper=0;                                   //upper bound for 2^k calc
    int *displ;                                     //displacements for proc
    int *sublength;                                 //vectorlength for proc
    int vectorlength=0;                             //maximum vector size

    double startTime=WallTime();                    //timestamp of program start
    double endTime=0.0;                             //timestamp of program end
    double S=pow(mathPi(),2)/6;                     //reference value S
    double Sn=0.0;                                  //approximated Sn
    double diff=0.0;                                //difference S-Sn
    double* globalsum=(double*)malloc(sizeof(double));
    double* sendvec=NULL;                           //MPI send vector
    double* receivevec=NULL;                        //MPI receive vector
    double* localsum=NULL;                         //sum of all partial sums

    MPI_Status status;

    /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "------precalculated values------\n");
    fprintf(stdout, "-- PI: %f\n", mathPi());
    fprintf(stdout, "-- S: %f\n", S);
    fprintf(stdout, "--------------------------------\n");
    */

    //check for correct number of input arguments
    if(argc<3 || argc>3)
    {
        fprintf(stdout, "\nusage: ex4 <lowerbound> <upperbound>\nexiting...\n");
        return EXIT_FAILURE;
    }

    init_app(argc, argv, &rank, &size);

    klower=atoi(argv[1]);
    kupper=atoi(argv[2]);
    vectorlength=pow(2,kupper); //maximum vector length 2^k

    //fprintf(stdout, "vectorlength: %d\n", vectorlength);

    sendvec=(double*)malloc(vectorlength*sizeof(double));
    receivevec=(double*)malloc(vectorlength*sizeof(double));
    localsum=(double*)malloc(sizeof(double));

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
        //split vector for every k
        splitVector(pow(2,i), size, &sublength, &displ);
        //send partial vectors to every proc
        MPI_Scatterv(sendvec, sublength, displ, MPI_DOUBLE, receivevec,
                     sizeof(double)*vectorlength, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //calculate local sum on every proc
        *localsum=0;
        (*localsum)=sum(receivevec, sublength[rank]);
        //fprintf(stdout, "localsum proc %d: %f\n", rank, *localsum);

        //wait for every proc to have the partial sum calculated
        MPI_Barrier(MPI_COMM_WORLD);

        //summing up all local sums
        MPI_Allreduce(localsum, globalsum, sizeof(double),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //output of calculated data
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
