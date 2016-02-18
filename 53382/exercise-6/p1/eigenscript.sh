#!/bin/bash

COUNTER=0
n=50
file1='time_eigenO3.data'
file2='time_eigenO0.data'
OPTIMIZATION=-O3

rm -f $file1 $file2

cd ~/Helsinki/THPC/crepaldi_caike_ex06/ex6p1/lapackf90/src
make cleanall
make FFLAGS=-O3
echo "lapack library: Compilation ($OPTIMIZATION)---> DONE"
cd ..
cd ..
gfortran eigendouble.f90 -llapack -Llapackf90/src
echo "eigendouble.f90: Compilation ---> DONE"

echo "File: $file1"
#echo "Optimization = $OPTIMIZATION"

while [ $n -le 2000 ]; do
    ./a.out $n 8540585 >> $file1
    echo "Iteration $COUNTER ---> DONE (n=$n)"
    let COUNTER=COUNTER+1
    let n=n+50
done

COUNTER=0
n=50
OPTIMIZATION=-O0

cd ~/Helsinki/THPC/crepaldi_caike_ex06/ex6p1/lapackf90/src
make cleanall
make FFLAGS=-O0
echo "lapack library: Compilation ($OPTIMIZATION)---> DONE"
cd ..
cd ..
gfortran eigendouble.f90 -llapack -Llapackf90/src
echo "eigendouble.f90: Compilation ---> DONE"

echo "File: $file2"
#echo "Optimization = $OPTIMIZATION"

while [ $n -le 2000 ]; do
    ./a.out $n 8540585 >> $file2
    echo "Iteration $COUNTER ---> DONE (n=$n)"
    let COUNTER=COUNTER+1
    let n=n+50
done
