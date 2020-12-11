#!/bin/bash

g++ -O3 -o clustering-phi clustering.cpp 

cp clustering-phi ../../../../bin/
chmod g+rx ../../../../bin/clustering-phi
