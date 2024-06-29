#!/bin/bash
#SBATCH -n 8
cp s.c an_sc.c initial
cd initial
g++ an_sc.c
./a.out
g++ s.c
./a.out
g++ seq_initial.c
./a.out
cd ..
