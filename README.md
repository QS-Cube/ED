# QS<sup>3</sup> for ED

This is a program package, proposed in arXiv:xxxx.xxxx, based on the (thick-restart) Lanczos method for analyzing a spin-1/2 XXZ-type quantum spin models on spatially uniform lattices near the fully polarized states. All of calculations including eigenvalue problems, expectation values for one/two-body spin operators, static/dynamical spin structure factors with QS<sup>3</sup> is performed in the symmetry-adapted bases specified by the number <i>N</i><sub>↓</sub> of down spin and the wave-number <i>k</i> associated with the translational symmetry and without using the bit representation for specifying the spin configuration. Since this treatment, QS<sup>3</sup> supports large-scale quantum systems containing more than 1000 sites with dilute <i>N</i><sub>↓</sub>.

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

This package includes the Fortran source codes, samples, and manual. 
Intel Fortran compiler (ver. XX or later) + MKL libraries or gfrotran + LAPACK (ver. YY or later) + BLAS (ver. ZZ or later) are prerequisites.

# Compile

If you have your Git account, please clone this repository on your local computer: 

$ git clone https://github.com/QS-Cube/ED.git

If you do not have, please push "Code" button and "Download ZIP" on this web page and get "ED-main.zip" and unpack the zip: 

$ unzip ED-main.zip <br>
$ cd ED-main

A simple Makefile is prepared to build the executable file "QS3.exe", and the following procedures after the cloning/downloading will give the executable and execute sample programs

$ cd script <br>
$ ./make.sh <br>

After executing the sample, all results are stored in output directories.

# Developers

Hiroshi Ueda and Tokuro Shimokawa
