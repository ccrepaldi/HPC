#!/bin/bash

FILE=$1
COUNTER=0

rm -f times.data

while read line; do
    let COUNTER=COUNTER+1
    done < $FILE
WL=0

while read line; do
    let WL=WL+1
    FFLAGS=$line
    gfortran $FFLAGS -o prog maximum.f90
    ./prog >> times.data
    echo "command $FFLAGS done: ($WL/$COUNTER)"
done < $FILE
