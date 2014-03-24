#!/bin/bash

MAXPROBLEMSIZE_P1=4096
MAXPROBLEMSIZE_P2=16384
MAXTHREADS=12
MAXPROCESSES=36

# args: mpiProcesses, ompThreads, problemSize
function runParallel { 
	OMP_NUM_THREADS=$2  mpiExec -np $1 timingTest/ex6parallel $3
}

function runTest1 {
	for (( size = 4; size <= $MAXPROBLEMSIZE_P1 ; size=size*2 )); do
		for (( processes = 1; processes <= $MAXPROCESSES; processes++ )); do
			for (( threads = 1; threads <= $MAXTHREADS; threads++)); do
				runParallel $processes $threads $size
			done
		done
	done
}

function runTest2 {
	for (( size = 4; size <= $MAXPROBLEMSIZE_P2 ; size=size*2 )); do
		runParallel $processes $threads $size
	done
}

function runTest3 {
	runParallel 1 $MAXTHREADS $MAXPROBLEMSIZE_P2
	runParallel $MAXPROCESSES 1 $MAXPROBLEMSIZE_P2
	runParallel $MAXPROCESSES $MAXTHREADS $MAXPROBLEMSIZE_P2
}



echo "# Timing test output is: numMpiProcess numOpenmpThreads problemSize time maxError"
runTest1
echo "# Constant number of MPI processes and OpenMP threads. 36 processes, 12 threads"
runTest2
echo  "# Mixed vs single speedup test"
runTest3


# echo "Done!"