#!/bin/sh
cd ../source
make
cd ../source_DSFLan
gfortran -o ../output/DSFLan.exe -Ofast main.f90 
cd ..
./QS3.exe < ./input/input_0.dat
cd output
./DSFLan.exe
