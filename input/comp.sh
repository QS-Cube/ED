rm list_*.dat condition.dat
rm *~
ifort make-list-xxz-cubic.f90
./a.out
rm a.out
