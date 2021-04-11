# QS<sup>3</sup> for ED

This is a program package, proposed in arXiv:xxxx.xxxx, based on the (thick-restart) Lanczos method for analyzing a spin-1/2 XXZ-type quantum spin models on spatially uniform lattices near the fully polarized states. All of calculations including eigenvalue problems, expectation values for one/two-body spin operators, static/dynamical spin structure factors with QS<sup>3</sup> is performed in the symmetry-adapted bases specified by the number <i>N</i><sub>↓</sub> of down spin and the wave-number <i>k</i> associated with the translational symmetry and without using the bit representation for specifying the spin configuration. Since this treatment, QS<sup>3</sup> supports large-scale quantum systems containing more than 1000 sites with dilute <i>N</i><sub>↓</sub>.

# Requirements

This package includes the Fortran source codes, samples, and manual. Intel Fortran compiler (ver. XX or later) + MKL libraries or gfrotran + LAPACK (ver. YY or later) + BLAS (ver. ZZ or later) are prerequisites.

# Compile

A simple Makefile is prepared to build the executable, and the following procedures after downloading the QS3.tar.gz will give the executable file, QS3.exe.

$ tar -zxvf QS3.tar.gz

$ cd QS3/source

$ make

# How to run
