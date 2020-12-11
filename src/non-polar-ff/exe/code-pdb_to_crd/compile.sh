#!/bin/bash

g++ -O2 -o pdb_to_crd pdb_to_crd.cpp ff.cpp 

cp pdb_to_crd ../../../../bin/
chmod g+rx ../../../../bin/pdb_to_crd
