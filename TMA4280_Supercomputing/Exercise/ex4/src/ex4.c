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
    double* receivevec;
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
            //      fprintf(stdout, "sendvec[%d]: %f", i, sendvec[i-1]);
        }
        /*---DECOMMENT FOR DEBUGGING PURPOSES---
          fprintf(stdout,"memory allocated and calculated values\n");
        */
    }

    for(i=klower; i<=kupper; i++)
    {
        splitVector(pow(2,i), size, &sublength, &displ); //split vector for every k
        if(rank==0) //master process
        {
            /*---DECOMMENT FOR DEBUGGING PURPOSES---
            	for(dbgloop=0;dbgloop<pow(2,i);dbgloop++)
            	{
            	  fprintf(stdout, "sendvec[%d]: %f\n", dbgloop, sendvec[dbgloop]);
            	}
            	  fprintf(stdout, "k: %d size: %d\n", i, size);
            */
            /*---DECOMMENT FOR DEBUGGING PURPOSES---*/
                        	for (j=0; j < size; j++)
                        	{

                        	    fprintf(stdout, "sublength[%d]: %d, displacement[%d]: %d\n", j, sublength[j], j, displ[j]);
                        	}
                      /* */
            for (j=0; j < size; j++) //send for each slave process
            {
                /*---DECOMMENT FOR DEBUGGING PURPOSES---*/
                fprintf(stdout, "---send for proc %d\n", j);
                /* */

                double* vsend=&(sendvec[displ[j]]);
                /*---DECOMMENT FOR DEBUGGING PURPOSES---
                		  for(dbgloop=0; dbgloop<sublength[j]; dbgloop++)
                		  {
                			fprintf(stdout, "----process %d\nvsend[%d]=%f\n", j, dbgloop, vsend[dbgloop]);
                		  }
                */
          //      fprintf(stdout, "---data for proc %d just NOT sent\n", j);
          //      fprintf(stdout, "------proc %d sublength %d displ %d address %x\n", j, sublength[j], displ[i], vsend);
                MPI_Send(vsend, sublength[j], MPI_DOUBLE, j, tag, MPI_COMM_WORLD);
          //      fprintf(stdout, "---data for proc %d sent\n", j);

            }
        }


        receivevec=(double*)malloc(sizeof(double)*sublength[rank]);
        fprintf(stdout, "proc %d receivevec created\n", rank);
        if(receivevec==NULL)
        {
            fprintf(stderr, "could not allocate memory: receivevec\n aborting...");
            return EXIT_FAILURE;
        }
        MPI_Recv(receivevec, sublength[rank], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        /*---DECOMMENT FOR DEBUGGING PURPOSES---*/
        fprintf(stdout, "process %d: data received\n", rank);
        /**/
        (*localsum)=sum(receivevec, sublength[rank]);


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

    close_app();


    return EXIT_SUCCESS;
}

//! \brief calculating the number pi
//! \return pi
double mathPi()
{
    return 4*atan(1);
}
