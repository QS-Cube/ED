#!/bin/sh
#################################
# ARG1 = ifort / gfortran
# ARG2 = lapack / mkl
#################################
ARG1=gfortran 
ARG2=mkl
#################################
echo "cd ../"
cd ../
echo "./configure FC=${ARG1} --with-lapak=${ARG2}"
./configure FC=${ARG1} --with-lapack=${ARG2}
echo "make"
make
echo "mv source/QS3 QS3.exe"
mv source/QS3 ./QS3.exe 1> /dev/null 2>&1
echo "mv source_only_u1/QS3_only_u1 ./QS3_only_u1.exe"
mv source_only_u1/QS3_only_u1 ./QS3_only_u1.exe 1> /dev/null 2>&1
echo "mv source_DSFLan/DSFLan ./DSFLan.exe"
mv source_DSFLan/DSFLan ./DSFLan.exe 1> /dev/null 2>&1

echo "**************"
echo  "  ex 1) FM HB model on the 6x6 square lattice"
echo  "       with u(1) and translational symmetry"
echo  "       Solver: Thick-restart Lanczos"
echo "**************"
echo "mkdir output_ex1 work"
mkdir output_ex1 work
echo "./QS3.exe < ./input_ex1/input.dat > output_ex1/output.dat"
./QS3.exe < input_ex1/input.dat > output_ex1/output.dat
echo "cd output_ex1"
cd output_ex1
echo "../DSFLan.exe"
../DSFLan.exe
echo "cd ../"
cd ../

echo "**************"
echo  "  ex 2) run ex 1)"
echo  "       Solver: Full diagonalization"
echo "**************"
echo "mkdir output_ex2"
mkdir output_ex2
echo "./QS3.exe < ./input_ex2/input.dat > output_ex2/output.dat"
./QS3.exe < input_ex2/input.dat > output_ex2/output.dat
echo "cd output_ex2"
cd output_ex2
echo "../DSFLan.exe"
../DSFLan.exe
echo "cd ../"
cd ../


echo "**************"
echo "  ex 3) AFM HB model on the 10x10x10 cubic lattice" 
echo "       with u(1) and translational symmetry"
echo "       Solver: Lanczos"
echo "**************"
echo "mkdir output_3"
mkdir output_ex3
echo "./QS3.exe < ./input_ex3/input.dat > output_ex3/output.dat"
./QS3.exe < input_ex3/input.dat > output_ex3/output.dat
echo "cd output_ex3"
cd output_ex3
echo "../DSFLan.exe"
../DSFLan.exe
echo "cd ../"
cd ../

echo "**************"
echo "  ex 4) run ex1) by use of the only U(1) symmetry"
echo "       Solver: Lanczos"
echo "**************"
echo "mkdir output_ex4"
mkdir output_ex4
echo "./QS3_only_u1.exe < ./input_ex4/input.dat > output_ex4/output.dat"
./QS3_only_u1.exe < input_ex4/input.dat > output_ex4/output.dat
echo "cd output_ex4"
cd output_ex4
echo "../DSFLan.exe"
../DSFLan.exe
echo "cd ../"
cd ../

echo "**************"
echo  "  ex 5) AFM HB model on the 6x6 triangular lattice"
echo  "       with u(1) and translational symmetry"
echo  "       Solver: Lanczos"
echo "**************"
echo "mkdir output_ex5"
mkdir output_ex5
echo "./QS3.exe < input_ex5/input.dat > output_ex5/output.dat"
./QS3.exe < input_ex5/input.dat > output_ex5/output.dat

echo "**************"
