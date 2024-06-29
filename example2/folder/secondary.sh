#!/bin/bash
#SBATCH -n 16
chmod +x op
rm -rf secondary
mkdir secondary
cp secondary.c sx.c cat.c all.sh op config1.dat secondary
cd initial
cp ch.dat stem.dat cs1.dat stem_kissing.dat ../secondary
cd ..
cd secondary
./op
bash all.sh
cd ..
