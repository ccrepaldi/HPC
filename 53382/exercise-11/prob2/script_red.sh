#!/usr/local/bin/bash

rm -f prob2.output

for i in {1..20}
do
    ./a.out $i 10000000 >> prob2.output
    echo "Threads: $i -> Done!"
done
