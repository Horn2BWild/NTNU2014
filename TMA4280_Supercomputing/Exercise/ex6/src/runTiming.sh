#!/bin/bash

MAXPROBLEMSIZE=1024
MAXTHREADS=4
MAXPROCESSES=4

function buildRelease {
	if [ ! -d "release" ]; then
		mkdir release
	fi
	cd release
	cmake -DCMAKE_BUILD_TYPE=Release .. -DTIMING_TEST=ON
	TIMING_TEST=1 make
	cd ..
}

# args: mpiProcesses, ompThreads, problemSize
function runParallel { 
	OMP_NUM_THREADS=$2  mpiExec -np $1 release/ex6parallel $3
}

function runTests {
	for (( size = 2; size <= $MAXPROBLEMSIZE; size=size*2 )); do
		for (( processes = 1; processes <= $MAXPROCESSES; processes++ )); do
			for (( threads = 1; threads <= $MAXTHREADS; threads++ )); do
				runParallel $processes $threads $size
			done
		done
	done
}

echo "Timing test output is: numMpiProcess numOpenmpThreads problemSize time"

buildRelease
runTests

echo "Done!"