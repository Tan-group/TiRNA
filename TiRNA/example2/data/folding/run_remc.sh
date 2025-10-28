#!/bin/bash
#SBATCH -n 16
cd ..
cd initial-file
bash initial.sh
cd ..
cd folding
cp config1.dat ../scoring-function
cp RNA_type ../wham
cp RNA_type ../optimizing
bash secondary.sh
gcc -O3 -Wall -fopenmp TiRNA_remc.c -o TiRNA_remc -lm
./TiRNA_remc
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
rm conf*.dat Bp*.dat Energy*.dat sec_stru*.dat RNA_type wham*.c
rm -rf secondary tm
cd ..
cd optimizing
rm conf*.dat second_stru*.dat cconf_0.dat stem*.dat
cd ..
cd rebuild-CG-to-AA
rm cs.dat rebuild RNA_type secondary state.dat stem*.dat
