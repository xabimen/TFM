#!/bin/bash

tail -n 11000 energy.dat > temp.dat
sed 's/\./,/g' temp.dat > temp2.dat
mean=$(datamash -W mean 2 < temp2.dat)
stdev=$(datamash -W sstdev 2 < temp2.dat)
temp=$(grep 'Temp' input_parameters.txt | cut -d ' ' -f 2)
rm temp.dat temp2.dat
echo $temp $mean $stdev
