#!/bin/sh
cd ../
for i in {0..41}; do
echo $i
./QS3.exe < input/input_$i.dat > output/output_$i.dat
done
