#!/bin/bash
cp wham.c A_state.c model
cd model
mkdir tm
cd tm
mkdir fragment
cd ..
gcc -Wall A_state.c -o A_state -lm
./A_state
for((i=0;i<15;i++))
do
   cp Energy_${i}.dat tm/fragment
   mv bp_${i}.dat tm/fragment
done
cp wham.c tm
cd tm
g++ wham.c
./a.out
cp thermal_stability.dat ../../../result/Thermal_Stability
cd ../../
cp center.c tc.c A_state1.c model
cd model
mv ch.dat ch_0.dat
for((i=0;i<15;i++))
do
	cat conf_${i}.dat >> conf_all.dat
done
g++ center.c 
./a.out
g++ tc.c
./a.out
mv cf.pdb folding_trj.pdb
rm conf_all.dat cconf_0.dat
mv folding_trj.pdb ../../result/Folding_trajectory
gcc -Wall A_state1.c -o A_state1 -lm
./A_state1
for((i=0;i<15;i++))
do
	mv secondary_stru_${i}.dat ../../result/Thermal_Stability
done
cd ..
	
	

