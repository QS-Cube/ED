# QS<sup>3</sup> for ED

This is a program package based on the (thick-restart) Lanczos method for analyzing a spin-1/2 XXZ-type quantum spin models on spatially uniform lattices near the fully polarized states. All of calculations including eigenvalue problems, expectation values for one/two-body spin operators, static/dynamical spin structure factors with QS<sup>3</sup> is performed in the symmetry-adapted bases specified by the number <i>N</i><sub>↓</sub> of down spin and the wave-number <i>k</i> associated with the translational symmetry and without using the bit representation for specifying the spin configuration. Since this treatment, QS<sup>3</sup> supports large-scale quantum systems containing more than 1000 sites with dilute <i>N</i><sub>↓</sub>.

# Licence

MIT License

# External routines/libraries 

BLAS, LAPACK

# Nature of problems

Physical properties (such as the total energy, the magnetic moment, the two-point spin correlation function, the dynamical structure factor)

# Solution method

Application software based on the full diagonalization method?, and the exact diagonalization using the Lanczos and thick-restart Lanczos methods for quantum spin <i>S</i> = 1/2 models such as the XXZ model.

# Restrictions

Spin <i>S</i> = 1/2 systems with the U(1) and translational symmetries.

# Requirements

This package, containing the Fortran source codes, samples, and manual, is available. For the building, Fortran compiler with BLAS/LAPACK library is prerequisite. 

# Compile

For those who have their own Git accounts, simply clone the repository on their local computers: 

$ git clone https://github.com/QS-Cube/ED.git

Otherwise, go to the web page and click the "Code" button and "Download ZIP" to get "ED-main.zip". The zip file is unpacked as 

$ unzip ED-main.zip <br>
$ cd ED-main

A simple Makefile is provided to build the executable files "QS3.exe" for systems preserving the translational symmetry and "QS3\_only\_u1.exe" for systems without the translational symmetry.
The following procedures after the cloning or downloading will give the executable file and execute sample programs

$ cd script <br>
$ ./make.sh <br>

After executing samples, all results are stored in output directories.

# Developers

Hiroshi Ueda and Tokuro Shimokawa
