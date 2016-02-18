#!/bin/bash/

.PHONY: clean all

poisson_omp: poisson_omp.f90
	gfortran -o poisson_omp -fopenmp poisson_omp.f90

poisson_serial: poisson_serial.f90
	gfortran -o poisson_serial poisson_serial.f90

poisson_serial_N10: poisson_serial_N10.f90
	gfortran -o poisson_serial_N10 poisson_serial_N10.f90

all: poisson_omp poisson_serial poisson_serial_N10

clean:
	rm -f *.o *.mod
