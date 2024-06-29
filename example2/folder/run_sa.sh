#!/bin/bash
#SBATCH -n 16
rm -rf model
mkdir model
bash initial.sh
bash secondary.sh
cp TiRNA_sa.c config1.dat model
cd secondary
cp ch_0.dat ../model
cd ..
cd model
gcc -Wall TiRNA_sa.c -o TiRNA_sa -lm
./TiRNA_sa
cd ..
bash optimize.sh
bash re_scoring.sh
bash rebuild.sh
g++ tm.c
./a.out
