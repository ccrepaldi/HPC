#!/bin/bash

echo "Runing with $1 threads."

for i in $( seq 1 25  ); do
    ../src/poisson_omp $1 > output.dat
    tail -1 output.dat | cut -d':' -f2 >> times${1}.txt
    echo "Iter: $i done!"
done

rm -f output.dat
    
