#!/bin/sh
for i in {0..62}; do
if [ $i -lt 21 ]; then
    echo 2 $i `grep E_0 output_${i}.dat`
elif [ $i -lt 42 ]; then
    echo 3 $(( i - 21 )) `grep E_0 output_${i}.dat`
else
    echo 1 $(( i - 42 )) `grep E_0 output_${i}.dat`
fi
done
