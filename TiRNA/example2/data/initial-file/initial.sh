#!/bin/bash
#SBATCH -n 8
g++ an_sc.c
./a.out
g++ s.c
./a.out
g++ seq_initial.c
./a.out
cp RNA_type ../folding
