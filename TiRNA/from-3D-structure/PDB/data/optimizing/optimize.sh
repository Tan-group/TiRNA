#!/bin/bash
#SBATCH -n 16
cd ..
cd folding/secondary
cp ch_0.dat stem.dat stem_kissing.dat ../../optimizing
cd ../
cp config1.dat ../optimizing
cd ..
cd optimizing
gcc -O3 -Wall -fopenmp TiRNA_optimize.c -o TiRNA_optimize -lm
./TiRNA_optimize
for((i=0;i<16;i++))
do
	cat conf_${i}.dat >> conf_all.dat
done
g++ center.c
./a.out
g++ tc.c
./a.out




