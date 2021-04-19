#!/bin/sh
echo "cd ../source"
cd ../source
echo "make"
make
echo "cd ../source_only_u1"
cd ../source_only_u1
echo "make"
make
echo "cd ../"
cd ../

echo "**************"
echo  "  ex ) FM HB model on the 6x6 square lattice"
echo  "       with u(1) and translational symmetry"
echo  "       Solver: Thick-restart Lanczos"
echo "**************"
echo "mkdir output work"
mkdir output work
echo "./QS3.exe < ./input/input.dat > output/output.dat"
./QS3.exe < ./input/input.dat > output/output.dat

echo "**************"
echo  "  ex0) run ex )"
echo  "       Solver: Full diagonalization"
echo "**************"
echo "mkdir output_0"
mkdir output_0
echo "./QS3.exe < ./input_0/input.dat > output_0/output.dat"
./QS3.exe < ./input_0/input.dat > output_0/output.dat

echo "**************"
echo "  ex1) AFM HB model on the 10x10x10 cubic lattice" 
echo "       with u(1) and translational symmetry"
echo "       Solver: Lanczos"
echo "**************"
echo "mkdir output_1"
mkdir output_1
echo "./QS3.exe < ./input_1/input.dat > output_1/output.dat"
./QS3.exe < ./input_1/input.dat > output_1/output.dat

echo "**************"
echo "  ex2) run ex1) by use of the only U(1) symmetry"
echo "       Solver: Lanczos"
echo "**************"
echo "mkdir output_2"
mkdir output_2
echo "./QS3_only_u1.exe < ./input_2/input.dat > output_2/output.dat"
./QS3_only_u1.exe < ./input_2/input.dat > output_2/output.dat

echo "**************"
#cd ../source_DSFLan
#gfortran -o ../output/DSFLan.exe -Ofast main.f90 
