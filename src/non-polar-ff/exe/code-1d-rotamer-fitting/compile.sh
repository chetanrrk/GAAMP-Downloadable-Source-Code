#!/bin/bash


mpicxx -O3 -o 1d-rotamer-fitting 1D-rotamer-fitting.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp 1d-rotamer-fitting ../../../../bin
chmod g+rx ../../../../bin/1d-rotamer-fitting
