#!/bin/bash

function buildTiming {
	if [ ! -d "timingTest" ]; then
		mkdir timingTest
	fi
	cd timingTest
	CXX=icpc CC=icc FC=ifort cmake -DCMAKE_BUILD_TYPE=Release .. -DTIMING_TEST=ON
	TIMING_TEST=1 make
	cd ..
}

module load cmake
module load intelcomp/13.0.1
module load openmpi/1.4.3-intel
module load lapack

buildTiming