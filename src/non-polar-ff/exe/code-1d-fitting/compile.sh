#!/bin/bash


mpicxx -O3 -o 1d-fitting 1D-fitting.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp 1d-fitting ../../../../bin/
chmod g+rx ../../../../bin/1d-fitting
