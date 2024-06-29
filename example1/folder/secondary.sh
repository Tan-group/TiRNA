#!/bin/bash
#SBATCH -n 16
chmod +x op
rm -rf secondary
mkdir secondary
cp secondary.c sx.c cat.c all.sh op secondary
cd model
cp top1.pdb ch.dat ../secondary
cd ..
cp ch.c sdCG_ch.c secondary
cd secondary
gcc -Wall secondary.c -o secondary -lm
./secondary
cp RNA_type ../
cp state.dat ../model
g++ ch.c
./a.out
./op
bash all.sh
cd ..
