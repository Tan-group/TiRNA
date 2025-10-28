#!/bin/bash
cd ..
cd optimizing
N=`ls *.pdb | wc -l`
cd ..
for((i=1;i<=N;i++))
do
 	cd optimizing
 	mv CG_top${i}.pdb ../rebuild-CG-to-AA
 	cd ..
 	cd rebuild-CG-to-AA
 	mv CG*.pdb CG.pdb
 	bash run.sh
 	mv All_atom.pdb All_atom_top${i}.pdb
 	mv CG.pdb CG_top${i}.pdb
 	mv sec_struc.dat sec_struc_top${i}.dat
 	mv CG_top${i}.pdb ../../result/CG_structure
 	mv All_atom_top${i}.pdb ../../result/All_atom_structure
 	mv sec_struc_top${i}.dat ../../result/Secondary_structure
 	cd ..
done
cd rebuild-CG-to-AA
