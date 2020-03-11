#!/bin/bash

N=(500, 1000, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000)

Nfile=(00500 01000 02000 04000 06000 08000 10000 12000 14000 16000)

echo "Testing for 2 cores"

for i in {0..10}
do
    time ./galsim ${N[$i]} input_data/ellipse_N_${Nfile[$i]}.gal 200 1e-5 0.255 0
done
