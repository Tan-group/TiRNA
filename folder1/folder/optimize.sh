#!/bin/bash
#SBATCH -n 16
rm -rf optimize
mkdir optimize
cd secondary
cp ch_0.dat stem.dat stem_kissing.dat ../optimize
cd ..
	
cp TiRNA_optimize.c center.c tc.c config1.dat optimize
cd optimize
gcc -Wall -fopenmp TiRNA_optimize.c -o TiRNA_optimize -lm
./TiRNA_optimize
for((i=0;i<16;i++))
do
	cat conf_${i}.dat >> conf_all.dat
done
g++ center.c
./a.out
g++ tc.c
./a.out
cd ..




