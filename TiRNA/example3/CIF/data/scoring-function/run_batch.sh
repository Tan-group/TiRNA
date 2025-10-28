#!/bin/bash

ulimit -s unlimited
 N=`ls *.pdb | wc -l`
 ./cgRNASP ./ ${N} cgRNASP.txt
