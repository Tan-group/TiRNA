#!/bin/bash
#SBATCH -n 16
chmod +x op
rm -rf secondary
mkdir secondary
cp secondary.c sx.c cat.c all.sh op secondary
cp top1.pdb ch.dat secondary
cp ch.c sdCG_ch.c secondary
cd secondary
gcc -Wall secondary.c -o secondary -lm
./secondary
cp RNA_type ../
cp state.dat ../
g++ ch.c
./a.out
./op 
bash all.sh
