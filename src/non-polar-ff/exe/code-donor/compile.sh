#!/bin/bash

icc -O3 -o donor mol_H_donor.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp donor ../../../../bin
chmod g+rx ../../../../bin/donor
