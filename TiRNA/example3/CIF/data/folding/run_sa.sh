#!/bin/bash
#SBATCH -n 1
cp config1.dat ../scoring-function
cp config1.dat ../optimizing
gcc -O3 -Wall TiRNA_mcsa.c -o TiRNA_mcsa -lm
./TiRNA_mcsa
g++ convert_pdb.c
./a.out
cd ..
cd scoring-function
bash scoring.sh
cd ..
cd folding
bash secondary.sh
cp RNA_type ../wham
cp RNA_type ../optimizing
cd ..
cd optimizing
bash optimize.sh
cd ..
cd scoring-function
bash re-scoring.sh
cd ..
cd rebuild-CG-to-AA
bash rebuild.sh
cd ..
cd wham
g++ tm.c
./a.out
cd ..
cd folding
