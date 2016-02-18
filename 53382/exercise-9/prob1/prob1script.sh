#!/bin/bash
for i in {2..24}
do
    mpirun -np $i prog >> prob1.out
    echo "-n = $i Done!"
done
