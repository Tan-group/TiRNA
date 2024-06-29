#!/bin/bash
cp 111.py -r data cgRNASP.c pai.c run_batch.sh ../model
cd ..
cd model
rm -rf xx
mkdir xx
mv cf.pdb xx
python3 111.py
gcc cgRNASP.c -lm -o cgRNASP
cp -r data cgRNASP run_batch.sh pai.c xx
cd xx
rm cf.pdb
bash run_batch.sh
g++ pai.c
./a.out
bash tiao.sh
cd ../
mv cf-*.pdb top1.pdb
rm cg* pai.c 111.py run_batch.sh
rm -rf data xx
cd ..
