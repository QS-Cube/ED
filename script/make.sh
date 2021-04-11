#!/bin/sh
cd ../source
make
cd ../script
../a.out < ../input/input.dat
