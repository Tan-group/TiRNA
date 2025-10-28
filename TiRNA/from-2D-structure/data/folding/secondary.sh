#!/bin/bash
#SBATCH -n 16
chmod +x op
rm -rf secondary
mkdir secondary
cp secondary.c sx.c cat.c all.sh op config1.dat secondary
cd ../
cd initial-file
cp ch.dat stem.dat cs1.dat stem_kissing.dat ../folding/secondary
cd ..
cd folding/secondary
./op
bash all.sh
cp ch_0.dat ../
cd ../
