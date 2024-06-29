#!/bin/bash
cd rebuild
cd ..
cd optimize
N=`ls *.pdb | wc -l`
cd ..
for((i=1;i<=N;i++))
do
 	cd optimize
 	mv CG_top${i}.pdb ../rebuild
 	cd ..
 	cd rebuild
 	mv CG*.pdb CG.pdb
 	bash run.sh
 	mv All_atom.pdb All_atom_top${i}.pdb
 	mv CG.pdb CG_top${i}.pdb
 	mv sec_struc.dat sec_struc_top${i}.dat
 	mv CG_top${i}.pdb ../../result/CG_strucure
 	mv All_atom_top${i}.pdb ../../result/All_atom_structure
 	mv sec_struc_top${i}.dat ../../result/Secondary_structure
 	cd ..
done
