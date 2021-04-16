#!/bin/sh
cd ../source
make
cd ..
#**************
#  ex1) FM HB model on the 6x6 square lattice 
#       with u(1) and translational symmetry
#       Solver: Thick-restart Lanczos
#**************
mkdir output work
./QS3.exe < ./input/input.dat
#cd output
#./DSFLan.exe
#cd ../
#**************
#  ex2) AFM HB model on the 10x10x10 cubic lattice 
#       with u(1) and translational symmetry
#       Solver: Lanczos
#**************
mkdir output_1
./QS3.exe < ./input_1/input.dat
#cd output_1
#./DSFLan.exe
#cd ../
#**************

#cd ../source_DSFLan
#gfortran -o ../output/DSFLan.exe -Ofast main.f90 
