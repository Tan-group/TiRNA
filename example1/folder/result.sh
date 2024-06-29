#!/bin/bash
folder_to_delete="secondary"
if [ -d "folder_to_delete" ]; then
rm -rf "folder_to_delete"
fi
mkdir result
for((i=1;i<=10;i++))
do
	cd optimize/${i}
	cat conf.dat >> ../../result/conf_all.dat
	cd ../..
done

cd secondary
cp ch.dat ../result
cd ..
cp tc.c center.c result

cd result
g++ center.c
./a.out
g++ tc.c
./a.out
rm conf_all.dat cconf_0.dat a.out tc.c center.c ch.dat

rm -rf secondary
