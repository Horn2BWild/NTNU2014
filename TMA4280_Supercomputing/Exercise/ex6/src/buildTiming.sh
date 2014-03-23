#!/bin/bash

function buildTiming {
	if [ ! -d "timingTest" ]; then
		mkdir timingTest
	fi
	cd timingTest
	cmake -DCMAKE_BUILD_TYPE=Release .. -DTIMING_TEST=ON
	TIMING_TEST=1 make
	cd ..
}

buildTiming