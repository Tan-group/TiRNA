#!/bin/bash
#SBATCH -n 16
rm -rf model
mkdir model
bash initial.sh
bash secondary.sh
cp TiRNA_remc.c config1.dat model
cd secondary
cp ch_0.dat ../model
cd ..
cd initial
cp RNA_type ../
cd ..
cd model
gcc -Wall -fopenmp TiRNA_remc.c -o TiRNA_remc -lm
./TiRNA_remc
cd ..
bash optimize.sh
bash re_scoring.sh
bash rebuild.sh
g++ tm.c
./a.out
