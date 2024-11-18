#!/bin/bash
#SBATCH -n 1
rm -rf model
mkdir model
bash initial.sh
cp TiRNA_sa.c config1.dat model
cd model
gcc -O3 -Wall TiRNA_sa.c -o TiRNA_sa -lm
./TiRNA_sa
cd ..
cp t1.c model
cd model
g++ t1.c
./a.out
cd ..
bash scoring.sh
bash secondary.sh
bash optimize.sh
bash re_scoring.sh
bash rebuild.sh
bash tm.sh
