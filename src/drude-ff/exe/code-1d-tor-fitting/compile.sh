#!/bin/bash


mpicxx -O3 -o 1d-fitting-drude 1D-fitting.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp 1d-fitting-drude ../../../../bin
chmod g+rx ../../../../bin/1d-fitting-drude
