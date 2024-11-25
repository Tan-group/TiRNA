#!/bin/bash
bash initial.sh
cp TiRNA_remc.c model
cd model
gcc -Wall -fopenmp TiRNA_remc.c -o TiRNA_remc -lm
./TiRNA_remc
cd ..
cp t1.c model
cd model
g++ t1.c
./a.out
cd ..
bash scoring.sh
bash secondary.sh
bash optimize.sh
bash scoring.sh
cd model.pdb
cp top1.pdb ../result
cd ..
