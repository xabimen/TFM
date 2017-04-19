#!/bin/bash

rm -r low_interaction/g_0.1.dat

for i in $(seq 1 2 9)
do
    echo $i
    sed -i a.bac "7s/[0-9]/$i/" input_parameters.txt
	cat input_parameters.txt
    ener=$(make | tail -n 1)
	echo $i $ener >> low_interaction/g_0.1.dat 
done
