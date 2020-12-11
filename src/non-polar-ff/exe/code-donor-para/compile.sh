#!/bin/bash


icc -O3 -o donor_para mol_H_donor.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp donor_para ../../../../bin/
chmod g+rx ../../../../bin/donor_para
