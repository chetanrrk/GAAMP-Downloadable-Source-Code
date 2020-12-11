#!/bin/bash


icc -O3 -o fitcharge-again fitcharge.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp fitcharge-again ../../../../bin/
chmod g+rx ../../../../bin/fitcharge-again

