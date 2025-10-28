#!/bin/bash
gcc -Wall secondary.c -o secondary -lm
./secondary
gcc -Wall rebuild.c -o rebuild -lm
./rebuild
