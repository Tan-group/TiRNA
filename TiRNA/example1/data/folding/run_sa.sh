#!/bin/bash
#SBATCH -n 1
cd ..
cd initial-file
bash initial.sh
cd ..
cd folding
cp config1.dat ../scoring-function
cp config1.dat ../optimizing
gcc -O3 -Wall TiRNA_mcsa.c -o TiRNA_mcsa -lm
./TiRNA_mcsa
g++ convert_pdb.c
./a.out
cp RNA_type ../wham
cd ..
cd scoring-function
bash scoring.sh
cd ..
cd folding
bash secondary.sh
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
rm conf*.dat Bp*.dat Energy*.dat sec_stru*.dat top1.pdb RNA_type wham*.c
rm -rf secondary tm
cd ..
cd optimizing
rm conf*.dat second_struc*.dat cconf_0.dat stem*.dat
cd ..
cd rebuild-CG-to-AA
rm cs.dat rebuild RNA_type secondary state.dat stem*.dat
