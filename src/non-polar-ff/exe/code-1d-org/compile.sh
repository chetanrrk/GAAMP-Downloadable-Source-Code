#!/bin/bash


mpicxx -O2 -o 1d-org 1D-fitting.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp 1d-org ../../../../bin
chmod g+rx ../../../../bin/1d-org

