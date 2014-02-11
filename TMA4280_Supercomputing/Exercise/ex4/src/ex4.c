/**************************************************************************** *
* Author: Andreas J. Hoermer                                                  *
* Email: andreas (at) hoermer.at                                              *
*                                                                             *
* Filename: ex4.c                                                             *
* LastChange: 07.02.2014                                                      *
*                                                                             *
* --------------------------------------------------------------------------- *
* PROGRAM DESCRIPTION                                                         *
* --------------------------------------------------------------------------- *
* This program is for summing up vectors with size 2^k containing values      *
* 1/2^k. For calculating the sum there are four options possible              *
*   --single processor, single thread                                         *
*   --single processor, multithreaded (using openMP)                          *
*   --multi processor, single thread (using MPI)                              *
*   --multi processor, multithreaded (using a openMP/MPI combination          *
*                                                                             *
* This program is intended to show different methods of distributed computing *
* --------------------------------------------------------------------------- *
* USAGE                                                                       *
* --------------------------------------------------------------------------- *
* ./ex4 <lower k> <upper k>                                                   *
*     -lower k ... lower bound of 2^k with console output                     *
*     -upper k ... upper bound of console output                              *
* mpirun -np <p> ./ex4 <lower k> <upper k>                                    *
*     -p ... number of proc                                                   *
* *************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

#include <omp.h>
#include <mpi.h>

#define BUFFERSIZE 128
#define TIMING 1             //set to 1 for enabling timing outputs
#define VERBOSE 0            //set to 1 for verbose output

//function prototypes
double mathPi();
double sum(double* vec, int length);
int errorHandlingMPI(int returnvalue, char* functionname);

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
    int returnvalue=0;                              //return value MPI functions

    double startTime=WallTime();                    //timestamp of program start
    double endTime=0.0;                             //timestamp of program end
    double currentTime=startTime;                   //timestamp of current step
    double S=pow(mathPi(),2)/6;                     //reference value S
    double Sn=0.0;                                  //approximated Sn
    double diff=0.0;                                //difference S-Sn
    double* globalsum=(double*)malloc(sizeof(double));
    double* sendvec=NULL;                           //MPI send vector
    double* receivevec=NULL;                        //MPI receive vector
    double* localsum=NULL;                         //sum of all partial sums
    MPI_Comm k_comm=MPI_COMM_WORLD;

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

    //getting lower and upper limit as numbers
    klower=atoi(argv[1]);
    kupper=atoi(argv[2]);
    vectorlength=pow(2,kupper); //maximum vector length 2^k

    //fprintf(stdout, "vectorlength: %d\n", vectorlength);

    sendvec=(double*)malloc(vectorlength*sizeof(double));
    receivevec=(double*)malloc(vectorlength*sizeof(double));
    localsum=(double*)malloc(sizeof(double));

    /*---DECOMMENT FOR DEBUGGING PURPOSES---
    fprintf(stdout, "------command line arguments------\n");
    fprintf(stdout, "-- lower bound: %d\n", klower);
    fprintf(stdout, "-- upper bound: %d\n", kupper);
    fprintf(stdout, "----------------------------------\n");
    */

    if(sendvec==NULL)
    {
        fprintf(stdout, "could not allocate memory: sendvec\n aborting...\n");
        return EXIT_FAILURE;
    }
    if(localsum==NULL)
    {
        fprintf(stdout, "could not allocate memory: localsum\n aborting...\n");
        return EXIT_FAILURE;
    }

#if TIMING
    fprintf(stdout, "Initialization time: %f\n", WallTime()-startTime);
    currentTime=WallTime();
#endif

    //calculating vector elements on root process
    if(rank==0)
    {
        //calculating multithreaded with openMP
#pragma omp parallel for schedule(guided,1)
        for(i=1; i<=vectorlength; i++)
        {
            sendvec[i-1]=1.0/pow(i,2);
        }
#if TIMING
        fprintf(stdout, "vector calculation time: %f\n", WallTime()-currentTime);
        currentTime=WallTime();
#endif
    }
    //calculate sum for every 2^i
    for(i=klower; i<=kupper; i++)
    {
      if(rank<pow(2,i)-1{
        printf("proc %d loop %d\n", rank, i);
        //split vector for every k
        splitVector(pow(2,i), size, &sublength, &displ);

        //        MPI_Comm_split(MPI_COMM_WORLD, i, rank, &k_comm);
        printf("proc %d after split vec\n", rank);

            printf("proc %d in if\n", rank);
            //send partial vectors to every proc

            returnvalue=MPI_Scatterv(sendvec, sublength, displ, MPI_DOUBLE,
                                     receivevec, sizeof(double)*vectorlength,
                                     MPI_DOUBLE, 0, MPI_COMM_WORLD);
            //MPI error handling

            printf("proc %d after scatter\n", rank);

            if (returnvalue != MPI_SUCCESS)
            {
                char error_string[BUFFERSIZE];
                int length_of_error_string, error_class;

                MPI_Error_class(returnvalue, &error_class);
                MPI_Error_string(error_class, error_string, &length_of_error_string);
                fprintf(stderr, "%3d: %s\n", rank, error_string);
                MPI_Error_string(returnvalue, error_string, &length_of_error_string);
                fprintf(stderr, "%3d: %s\n", rank, error_string);
                MPI_Abort(MPI_COMM_WORLD, returnvalue);
            }


                //calculate local sum on every proc
                *localsum=0;
                (*localsum)=sum(receivevec, sublength[rank]);
                printf("proc %d after calculation\n", rank);

            //wait for every proc to have the partial sum calculated

            //   returnvalue=MPI_Barrier(MPI_COMM_WORLD);
            printf("proc %d after barrier\n", rank);
            //MPI error handling
            if (returnvalue != MPI_SUCCESS)
            {
                char error_string[BUFFERSIZE];
                int length_of_error_string, error_class;

                MPI_Error_class(returnvalue, &error_class);
                MPI_Error_string(error_class, error_string, &length_of_error_string);
                fprintf(stderr, "%3d: %s\n", rank, error_string);
                MPI_Error_string(returnvalue, error_string, &length_of_error_string);
                fprintf(stderr, "%3d: %s\n", rank, error_string);
                MPI_Abort(MPI_COMM_WORLD, returnvalue);
            }

            //summing up all local sums
            returnvalue=MPI_Allreduce(localsum, globalsum, sizeof(double),
                                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            printf("proc %d after allreduce\n", rank);
            //MPI error handling

        if (returnvalue != MPI_SUCCESS)
        {
            char error_string[BUFFERSIZE];
            int length_of_error_string, error_class;

            MPI_Error_class(returnvalue, &error_class);
            MPI_Error_string(error_class, error_string, &length_of_error_string);
            fprintf(stderr, "%3d: %s\n", rank, error_string);
            MPI_Error_string(returnvalue, error_string, &length_of_error_string);
            fprintf(stderr, "%3d: %s\n", rank, error_string);
            MPI_Abort(MPI_COMM_WORLD, returnvalue);
        }
}
        //output of calculated data
        if(rank==0 && i>=klower && i<=kupper)
        {
            diff=S-(*globalsum);

#if VERBOSE
            fprintf(stdout,"----------Result----------\n");
            fprintf(stdout,"--k=%d elements=%d\n", i, (int)pow(2,i));
            fprintf(stdout,"--Sn=%f\n", *globalsum);
            fprintf(stdout,"--S=%f\n", S);
            fprintf(stdout,"----diff=%.10g\n", diff);
#else
            fprintf(stdout,"k=%d, diff=%.10g\n", i, diff);
#endif


#if TIMING
            fprintf(stdout, "sum up time k=%d: %f\n", i, WallTime()-currentTime);
            currentTime=WallTime();
#endif

        }


    }

    endTime=WallTime();
    fprintf(stdout, "total run time: %lf\n\n", endTime-startTime);

    close_app();
    return EXIT_SUCCESS;
}

//! \brief calculating the number pi
//! \return pi
double mathPi()
{
    return 4*atan(1);
}

//! \brief calculating sum over vector elements
//! \param vec input vector (double)
//! \param length length of the input vector
//! \return sum
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
