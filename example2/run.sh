#!/bin/bash
#SBATCH -n 16
rm -rf result
mkdir result
cd result
mkdir CG_strucure All_atom_structure Folding_trajectory Thermal_Stability Secondary_structure
cd ..
cp seq.dat folder/initial
cp config.dat folder
cd folder
python config.py
g++ main.c
./a.out
cd ..
