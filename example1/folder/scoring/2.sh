#!/bin/bash
cp 111.py -r data cgRNASP.c pai2.c run_batch.sh config1.dat ../optimize
cd ..
cd optimize
rm -rf xx
mkdir xx
mv cf.pdb xx
python3 111.py
gcc cgRNASP.c -lm -o cgRNASP
cp -r data cgRNASP run_batch.sh pai2.c config1.dat xx
cd xx
rm cf.pdb
bash run_batch.sh
g++ pai2.c
./a.out
bash tiao.sh
cd ../
rm cg* pai2.c 111.py run_batch.sh
rm -rf data xx
cd ..
