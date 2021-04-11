#!/bin/sh
i=42
for nod in {1..1}; do
while read kx ky kz; do
echo $kx $ky $kz
cat input_class.dat | \
sed -e "s/NNOD/$nod/g" |
sed -e "s/NKX/$kx/g" |
sed -e "s/NKY/$ky/g" |
sed -e "s/NKZ/$kz/g" > input_$i.dat
i=$(( i + 1 ))
done < ./path.dat
done
