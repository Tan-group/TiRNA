#!/bin/bash
#SBATCH -n 16
rm -rf result
mkdir result
cd result
mkdir CG_structure All_atom_structure Folding_trajectory Thermal_Stability Secondary_structure
cd ..
cp seq.dat data/initial-file
cp config.dat data/folding
cd data/folding
python config.py
python op_algorithm.py

