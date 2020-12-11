#!/bin/bash


mpicxx -O3 -o fitcharge-again-drude fitcharge.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp fitcharge-again-drude ../../../../bin/
chmod g+rx ../../../../bin/fitcharge-again-drude

