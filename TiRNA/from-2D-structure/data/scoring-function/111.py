#!/bin/python

import glob

for i in ["xx"]: 
    route = i + "/" + "cf.pdb"
    with open(route) as f1:
        lines=f1.readlines()
    
    j=0
    for l in lines:
        if l[:6]=='CRYST1':
            route = i + "/" +  "cf-" + str(j) + ".pdb"
            #route = i + "/" + i + "1-" + str(j) + ".pdb"
            j+=1
        elif l[:3] != 'TER' and l[:3] != 'END':
            with open(route,"a+") as F1:
                F1.write(l)



