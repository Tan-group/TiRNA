#!/bin/bash
for((i=0;i<15;i++))
do
	cat conf_${i}.dat >> conf_all.dat
done
gcc -Wall sx.c -o sx -lm
./sx
g++ cat.c
./a.out
